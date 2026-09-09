# =====================================================================================
#  Three-phase Volt-VAr droop in a distribution OPF
#  HOST   : IVACOPF, a three-phase current-voltage AC-OPF, linearised (Soltani, Khorsand
#                      & Ma, IEEE OJIA 5, 2024, doi:10.1109/OJIA.2024.3367547), solved by
#                      successive linearisation
#  DROOP  : Lambda / SOS2 (convex combination of breakpoints, adjacency forced by binaries)
#
#  Test system : network_5_Feeder_2, a real 415/240 V LV feeder from the ENWL dataset,
#                Kron-reduced to three wires. 194 buses, 193 lines, 489 m, and 18
#                single-phase loads 
#  Horizon     : 24 h at 15-minute resolution (96 steps), residential load shape.
#  Objective   : minimise total PV curtailment over the day.
#
#  Equation numbers in the comments below are the tutorial's:
#  https://epsrlab-ub.github.io/SmartInverter-3P-DOPF-POWERUP.jl/dev/tutorial_voltvar/#Three-phase-IVACOPF:-the-near-exact-host

# =====================================================================================

using JuMP, Gurobi                    
using JSON3, Printf, Plots, LinearAlgebra

const PHASES = 1:3
const DATA   = joinpath(@__DIR__, "data")

# =====================================================================================
# Interactive tutorial configuration
# =====================================================================================
# By default this script is interactive. For unattended runs, set TP_INTERACTIVE=0 and
# use the existing environment variables (TP_NPV, TP_STEPS, TP_WARMSTART, TP_IMAXSEG).
const INTERACTIVE = get(ENV, "TP_INTERACTIVE", "1") != "0"

function ask_int(prompt::AbstractString, default::Int; allowed = nothing)
    while true
        print("$prompt [$default]: ")
        s = strip(readline())
        x = isempty(s) ? default : tryparse(Int, s)

        if x === nothing
            println("  Please enter an integer.")
            continue
        end

        if allowed !== nothing && !(x in allowed)
            println("  Choose one of: $(join(collect(allowed), ", "))")
            continue
        end

        return x
    end
end

function ask_yesno(prompt::AbstractString, default::Bool = true)
    tag = default ? "Y/n" : "y/N"
    while true
        print("$prompt [$tag]: ")
        s = lowercase(strip(readline()))
        isempty(s) && return default
        s in ("y", "yes") && return true
        s in ("n", "no") && return false
        println("  Please enter y or n.")
    end
end

function pause_tutorial(message::AbstractString = "Press ENTER to continue...")
    INTERACTIVE || return
    println()
    print(message)
    readline()
    println()
end

# Four inverter classes. Each phase receives one of each when N_PV_PER_PHASE = 4, so
# size is not confounded with phase or distance from the substation. Because q̄ = S_max,
# the four classes follow four different droop curves with common voltage breakpoints.
const PV_CLASSES = (("A — 3 kW",  3.0),
                    ("B — 5 kW",  5.0),
                    ("C — 8 kW",  8.0),
                    ("D — 12 kW", 12.0))
const S_OVER_P       = 1.1
const VLIM           = (0.95, 1.05)
const VBP            = [0.88, 0.90, 0.97, 1.00, 1.02, 1.10]
const QSHAPE         = [1.0, 1.0, 0.0, 0.0, -1.0, -1.0]
const SBASE_KVA      = 100.0
const MIP_GAP        = 1e-3
const METHOD         = "Lambda / SOS2"
const FIGSUF         = "_iva"

# Successive-linearisation controls.
const MAX_ITER = parse(Int, get(ENV, "TP_MAXITER", "15"))
const TOL      = parse(Float64, get(ENV, "TP_TOL", "1e-6"))

# Defaults can still be overridden with environment variables.
npv_default      = parse(Int, get(ENV, "TP_NPV", "4"))
steps_default    = parse(Int, get(ENV, "TP_STEPS", "96"))
warm_default     = get(ENV, "TP_WARMSTART", "sweep") == "sweep"
imaxseg_default  = parse(Int, get(ENV, "TP_IMAXSEG", "0"))

if INTERACTIVE
    println()
    println("============================================================")
    println(" THREE-PHASE SMART-INVERTER OPF — INTERACTIVE TUTORIAL")
    println("============================================================")
    println()
    println("Configure a three-phase distribution-system simulation with")
    println("IEEE 1547 Volt-VAr control encoded using Lambda / SOS2.")
    println("Press ENTER at a prompt to accept the value shown in brackets.")
    println()

    global N_PV_PER_PHASE = ask_int(
        "Number of PV inverters per phase",
        npv_default in 1:4 ? npv_default : 4;
        allowed = 1:4,
    )

    global T = ask_int(
        "Number of time steps over 24 h",
        steps_default in (24, 48, 96) ? steps_default : 96;
        allowed = (24, 48, 96),
    )

    global WARMSTART = ask_yesno(
        "Use exact three-phase power-flow warm start?",
        warm_default,
    )

    thermal_default = imaxseg_default > 0
    enforce_thermal = ask_yesno(
        "Enforce thermal line-current limits?",
        thermal_default,
    )
    global IMAX_SEG = enforce_thermal ? (imaxseg_default > 0 ? imaxseg_default : 8) : 0

    println()
    println("---------------- Selected configuration ----------------")
    println("PV sites / phase : $N_PV_PER_PHASE")
    println("Time steps       : $T")
    println("Warm start       : $(WARMSTART ? "three-phase sweep" : "flat start")")
    println("Thermal limits   : $(IMAX_SEG > 0 ? "enforced with $IMAX_SEG-sided polygon" : "reported only")")
    println("--------------------------------------------------------")
    pause_tutorial("Press ENTER to load the feeder and build the simulation...")
else
    N_PV_PER_PHASE = npv_default
    T = steps_default
    WARMSTART = warm_default
    IMAX_SEG = imaxseg_default
end

T in (24, 48, 96) || error("TP_STEPS must be one of 24, 48, or 96 for the supplied 96-point profiles.")
N_PV_PER_PHASE in 1:4 || error("TP_NPV must be between 1 and 4 for this tutorial.")
# =====================================================================================

# ============================================================== 1) read the BMOPF case ==
# Feeder and horizon may be overridden from the environment, so the identical model can
# be run on a much larger network without touching the code:
#   TP_CASE=network_17_Feeder_6 TP_STEPS=24 julia --project=. <script>
const CASE = get(ENV, "TP_CASE", "network_5_Feeder_2")
net    = JSON3.read(read(joinpath(DATA, "$CASE.bmopf.json"), String))
loadpr = JSON3.read(read(joinpath(DATA, "load_profile_residential_15min.json"), String))
solar  = JSON3.read(read(joinpath(DATA, "solar_profile.json"), String))

BUSES  = [String(k) for k in keys(net.bus)]
bus_id = Dict(b => i for (i, b) in enumerate(BUSES))
nb     = length(BUSES)

SLACK      = String(net.voltage_source.source.bus)
VBASE      = Float64(net.voltage_source.source.v_magnitude[1])   # phase-to-neutral, V
SBASE      = SBASE_KVA * 1e3                                     # VA, per phase
ZBASE      = VBASE^2 / SBASE
IBASE      = SBASE / VBASE                                       # A, per phase
const VNOM = 1.0

step   = 96 ÷ T
Pmult  = collect(Float64, loadpr.P_percent)[1:step:end][1:T] ./ 100
Qmult  = collect(Float64, loadpr.Q_percent)[1:step:end][1:T] ./ 100
G      = collect(Float64, solar.G_percent)[1:step:end][1:T]  ./ 100

# ---- lines: (from, to) and their 3×3 series impedance in p.u. -------------------------
struct Branch
    from::String
    to::String
    Z::Matrix{ComplexF64}     # 3×3, per unit
    imax::Float64             # conductor rating, p.u. current
end

function linecode_Z(lc, len_m)
    Z = zeros(ComplexF64, 3, 3)
    for p in PHASES, q in PHASES
        r = Float64(lc[Symbol("R_series_$(p)_$(q)")])      # Ω per metre
        x = Float64(lc[Symbol("X_series_$(p)_$(q)")])
        Z[p, q] = complex(r, x) * len_m / ZBASE
    end
    return Z
end

BR = Branch[]
for (_, l) in pairs(net.line)
    lc = net.linecode[Symbol(l.linecode)]
    push!(BR, Branch(String(l.bus_from), String(l.bus_to),
                     linecode_Z(lc, Float64(l.length)),
                     Float64(lc.i_max[1]) / IBASE))
end
nbr = length(BR)

# The IVACOPF line equations (35) want the impedance split into its real and imaginary
# 3×3 parts. Note that no rotation operator appears anywhere: unlike LinDist3Flow, this
# host never assumes the three phases are near-balanced, so the raw coupled impedance is
# all it needs.
Rm = [real.(b.Z) for b in BR]
Xm = [imag.(b.Z) for b in BR]

# Shunt admittances y^{φp}: the Σ y V terms of (34) and (35). The ENWL feeders are
# supplied Kron-reduced to three wires with zero shunt, so these terms vanish here; they
# are kept explicit so a dataset that carries them can be dropped in.
YSH = [zeros(ComplexF64, 3, 3) for _ in BR]
const HAS_SHUNT = any(any(!iszero, y) for y in YSH)

# ---- loads: single-phase, phase-to-neutral, in p.u. -----------------------------------
# Pload[bus, phase] is the peak; the residential shape scales it over the day.
Pload_pk = zeros(nb, 3)
Qload_pk = zeros(nb, 3)
for (_, ld) in pairs(net.load)
    b = bus_id[String(ld.bus)]
    φ = parse(Int, String(ld.terminal_map[1]))           # ("1","n") → phase 1
    Pload_pk[b, φ] += Float64(ld.p_nom[1]) / SBASE
    Qload_pk[b, φ] += Float64(ld.q_nom[1]) / SBASE
end

# ================================================== 2) topology helpers and the sweep ===
out_br = [Int[] for _ in 1:nb]
in_br  = [Int[] for _ in 1:nb]
for (k, br) in enumerate(BR)
    push!(out_br[bus_id[br.from]], k)
    push!(in_br[bus_id[br.to]], k)
end
islack = bus_id[SLACK]

# The three-phase slack reference: 1∠0°, 1∠−120°, 1∠+120°. This is also the flat start of
# the linearisation, exactly as the paper prescribes. Seeding all three phases at 1∠0°
# instead is a silent and expensive mistake: the mutual terms then add rather than
# largely cancelling.
const V0 = ComplexF64[1, exp(-2π * im / 3), exp(2π * im / 3)]

# Farthest load buses on each phase, by |Z| along the path from the source.
function distance_from_slack()
    adj = Dict(b => Tuple{String,Int}[] for b in BUSES)
    for (k, br) in enumerate(BR)
        push!(adj[br.from], (br.to, k))
        push!(adj[br.to],   (br.from, k))
    end
    dist   = Dict(SLACK => 0.0)
    parent = Dict{String,Tuple{String,Int}}()
    order  = [SLACK]
    queue  = [SLACK]
    while !isempty(queue)
        b = popfirst!(queue)
        for (nbb, k) in adj[b]
            if !haskey(dist, nbb)
                dist[nbb]   = dist[b] + abs(BR[k].Z[1, 1])
                parent[nbb] = (b, k)
                push!(order, nbb); push!(queue, nbb)
            end
        end
    end
    return dist, parent, order
end
DIST, PARENT, ORDER = distance_from_slack()

CHILD_BR = [Int[] for _ in 1:nb]
for (k, br) in enumerate(BR); push!(CHILD_BR[bus_id[br.from]], k); end

# ---- exact three-phase backward/forward sweep -----------------------------------------
#  Kron-reduced three-wire model with implicit ground and zero shunts, so the sweep is
#  exact for these injections. It does double duty here: it produces the warm-start
#  linearisation point *before* the optimisation, and it audits the solution afterwards.
#  `Pd`/`Qd` are npv×T arrays of inverter dispatch in p.u.
function sweep_state(Pd, Qd, t; iters = 60)
    Vc  = [copy(V0) for _ in 1:nb]
    Ibr = [zeros(ComplexF64, 3) for _ in 1:nbr]
    rev = reverse(ORDER)
    for _ in 1:iters
        S = [zeros(ComplexF64, 3) for _ in 1:nb]
        for b in 1:nb, φ in PHASES
            S[b][φ] = complex(Pload_pk[b, φ] * Pmult[t], Qload_pk[b, φ] * Qmult[t])
        end
        for (i, g) in enumerate(PV)
            S[g.bus][g.phase] -= complex(Pd[i, t], Qd[i, t])
        end
        Ibus = [conj.(S[b] ./ Vc[b]) for b in 1:nb]        # current *drawn* by each bus
        for bname in rev                                   # backward: currents
            b = bus_id[bname]
            b == islack && continue
            (_, k) = PARENT[bname]
            Ibr[k] = Ibus[b] + sum((Ibr[c] for c in CHILD_BR[b]), init = zeros(ComplexF64, 3))
        end
        for bname in ORDER                                 # forward: voltages
            bname == SLACK && continue
            (par, k) = PARENT[bname]
            Vc[bus_id[bname]] = Vc[bus_id[par]] - BR[k].Z * Ibr[k]
        end
    end
    return Vc, Ibr
end

# ============================================================ 3) place the PV inverters ==
PV = NamedTuple{(:bus, :phase, :cls, :Pmax, :Smax),Tuple{Int,Int,Int,Float64,Float64}}[]
for φ in PHASES
    cand = [b for b in 1:nb if Pload_pk[b, φ] > 0]
    sort!(cand, by = b -> -DIST[BUSES[b]])
    for (rank, b) in enumerate(cand[1:min(N_PV_PER_PHASE, length(cand))])
        # rotate the class order by phase, so the farthest site is not always the same size
        ci = mod(rank - 1 + (φ - 1), length(PV_CLASSES)) + 1
        p  = PV_CLASSES[ci][2] * 1e3 / SBASE
        push!(PV, (bus = b, phase = φ, cls = ci, Pmax = p, Smax = S_OVER_P * p))
    end
end
npv = length(PV)
ncls = length(PV_CLASSES)

pv_at = [Int[] for _ in 1:nb, _ in 1:3]
for (i, g) in enumerate(PV); push!(pv_at[g.bus, g.phase], i); end

# irradiance ceiling per inverter per step
Pavail = [PV[i].Pmax * G[t] for i in 1:npv, t in 1:T]

@printf("network : %d buses, %d lines, %d single-phase loads\n", nb, nbr, length(net.load))
@printf("bases   : V = %.2f V (phase-neutral), S = %.0f kVA/phase, Z = %.4f Ω, I = %.1f A\n",
        VBASE, SBASE_KVA, ZBASE, IBASE)
@printf("load    : %.2f kW peak total\n", sum(Pload_pk) * SBASE / 1e3)
totkW = sum(g.Pmax for g in PV) * SBASE / 1e3
@printf("PV      : %d smart inverters in %d classes, %.1f kW total (%.0f%% of peak load)\n",
        npv, ncls, totkW, 100 * totkW / (sum(Pload_pk) * SBASE / 1e3))
@printf("host    : IVACOPF, %s start, tol %.0e on max(MAPB, MRPB, MVM), ≤ %d passes\n",
        WARMSTART ? "swept" : "flat", TOL, MAX_ITER)

println("\ninverter classes")
println("  class        P_rated   S_max    q̄ = S_max    q̄ (p.u.)   sites")
for ci in 1:ncls
    n = count(g -> g.cls == ci, PV)
    s = S_OVER_P * PV_CLASSES[ci][2]
    @printf("  %-11s  %5.1f kW  %5.2f kVA  %6.2f kvar    %.4f      %d\n",
            PV_CLASSES[ci][1], PV_CLASSES[ci][2], s, s, s * 1e3 / SBASE, n)
end

println("\nsmart-inverter sites ($(N_PV_PER_PHASE) electrically farthest load buses per phase)")
println("  bus   phase   class        |Z| from source   local load   P_rated    q̄")
println("                                    (Ω)           (kW)       (kW)    (kvar)")
for g in sort(PV, by = g -> (g.phase, -DIST[BUSES[g.bus]]))
    @printf("  %-5s   %d     %-11s   %8.4f        %6.2f     %6.1f   %6.2f\n",
            BUSES[g.bus], g.phase, PV_CLASSES[g.cls][1], DIST[BUSES[g.bus]] * ZBASE,
            Pload_pk[g.bus, g.phase] * SBASE / 1e3,
            g.Pmax * SBASE / 1e3, g.Smax * SBASE / 1e3)
end

# ================================================ 4) the linearisation point (the "∘") ==
#  Everything marked `_pr` below is a *constant* held from the previous pass, not a
#  variable. Pass 1 needs somewhere to start:
#    · flat  : 1∠0°, 1∠−120°, 1∠+120° with all currents zero, as the paper prescribes;
#    · sweep : an exact power flow at full PV and zero VArs. It costs one sweep per time
#              step and is a physically consistent state, which typically halves the
#              number of passes. This is the three-phase analogue of the package's
#              `warm_start = :lindistflow`.
v_r_pr    = zeros(nb, 3, T)
v_im_pr   = zeros(nb, 3, T)
Ibs_r_pr  = zeros(nb, 3, T)
Ibs_im_pr = zeros(nb, 3, T)

warm_seconds = @elapsed begin
    if WARMSTART
        Qzero = zeros(npv, T)
        for t in 1:T
            Vc, Ibr = sweep_state(Pavail, Qzero, t)
            for b in 1:nb, φ in PHASES
                v_r_pr[b, φ, t]  = real(Vc[b][φ])
                v_im_pr[b, φ, t] = imag(Vc[b][φ])
                inj = sum((Ibr[k][φ] for k in out_br[b]), init = zero(ComplexF64)) -
                      sum((Ibr[k][φ] for k in in_br[b]),  init = zero(ComplexF64))
                Ibs_r_pr[b, φ, t]  = real(inj)
                Ibs_im_pr[b, φ, t] = imag(inj)
            end
        end
    else
        for b in 1:nb, φ in PHASES, t in 1:T
            v_r_pr[b, φ, t]  = real(V0[φ])
            v_im_pr[b, φ, t] = imag(V0[φ])
        end
    end
end
WARMSTART && @printf("\nwarm start: %d exact sweeps in %.1f s\n", T, warm_seconds)

if INTERACTIVE
    println()
    println("============================================================")
    println(" SOLVING THE THREE-PHASE IVACOPF")
    println("============================================================")
    println()
    println("Each successive-linearisation pass will:")
    println("  1. build a linear approximation around the current state,")
    println("  2. solve the resulting mixed-integer optimization problem,")
    println("  3. update the voltage/current linearisation point, and")
    println("  4. evaluate MAPB, MRPB, and MVM against the nonlinear relations.")
    println()
    @printf("Convergence requires max(MAPB, MRPB, MVM) < %.1e.\n", TOL)
    pause_tutorial("Press ENTER to start the optimization...")
end

# ==================================================== 5) successive-linearisation loop ==
# One Gurobi environment, reused by every pass: the licence banner is printed once and
# each rebuild is cheaper than spinning up a fresh environment.
const GRB_ENV = Gurobi.Env()

iterlog = NamedTuple{(:iter, :seconds, :objective, :MAPB, :MRPB, :MVM, :status),
                     Tuple{Int,Float64,Float64,Float64,Float64,Float64,String}}[]
model = nothing; status = nothing; secs = 0.0; total_solve = 0.0
V = Vr = Vi = nothing
Pdg_v = Qdg_v = PVC_v = nothing
Ibr_r_v = Ibr_im_v = nothing
converged = false
nvar = nbin = ncon = 0

for it in 1:MAX_ITER
    global model, status, secs, total_solve, V, Vr, Vi, Pdg_v, Qdg_v, PVC_v
    global Ibr_r_v, Ibr_im_v, converged, nvar, nbin, ncon

    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_silent(model)
    set_optimizer_attribute(model, "MIPGap", MIP_GAP)

    # ---- host variables ---------------------------------------------------------------
    # v is the voltage *magnitude*, eq. (40), bounded by (41); v_r/v_im are the rectangular
    # network equations are actually written in.
    @variable(model, VLIM[1] <= v[1:nb, PHASES, 1:T] <= VLIM[2], start = VNOM)
    @variable(model, v_r[b = 1:nb, φ = PHASES, t = 1:T],  start = v_r_pr[b, φ, t])
    @variable(model, v_im[b = 1:nb, φ = PHASES, t = 1:T], start = v_im_pr[b, φ, t])
    @variable(model, Ibr_r[1:nbr, PHASES, 1:T])       # branch current, real part
    @variable(model, Ibr_im[1:nbr, PHASES, 1:T])      # branch current, imaginary part
    @variable(model, Ibs_r[b = 1:nb, φ = PHASES, t = 1:T],  start = Ibs_r_pr[b, φ, t])
    @variable(model, Ibs_im[b = 1:nb, φ = PHASES, t = 1:T], start = Ibs_im_pr[b, φ, t])
    @variable(model, Pg[PHASES, 1:T])                 # substation injection
    @variable(model, Qg[PHASES, 1:T])
    @variable(model, 0 <= Pdg[1:npv, 1:T])
    @variable(model, Qdg[1:npv, 1:T])
    @variable(model, 0 <= PVC[1:npv, 1:T])

    # ---- slack reference: 1∠0°, 1∠−120°, 1∠+120° --------------------------------------
    @constraint(model, [φ in PHASES, t in 1:T], v_r[islack, φ, t]  == real(V0[φ]))
    @constraint(model, [φ in PHASES, t in 1:T], v_im[islack, φ, t] == imag(V0[φ]))

    # ---- line current constraints, eq. (35): exact, linear, fully phase-coupled --------
    #  V_n − V_m = Z_nm I_nm  (shunt terms of (34) omitted, zero on this dataset), split
    #  into real and imaginary parts. Every mutual term Z[φ,ψ], ψ ≠ φ, is carried: no
    #  transposition, no balance and no sequence decomposition is assumed anywhere.
    @constraint(model, [k in 1:nbr, φ in PHASES, t in 1:T],
        v_r[bus_id[BR[k].from], φ, t] - v_r[bus_id[BR[k].to], φ, t] ==
            sum(Rm[k][φ, ψ] * Ibr_r[k, ψ, t] - Xm[k][φ, ψ] * Ibr_im[k, ψ, t] for ψ in PHASES))
    @constraint(model, [k in 1:nbr, φ in PHASES, t in 1:T],
        v_im[bus_id[BR[k].from], φ, t] - v_im[bus_id[BR[k].to], φ, t] ==
            sum(Rm[k][φ, ψ] * Ibr_im[k, ψ, t] + Xm[k][φ, ψ] * Ibr_r[k, ψ, t] for ψ in PHASES))

    # ---- bus current injection, eq. (36): KCL per bus and phase, exact and linear ------
    @constraint(model, [b in 1:nb, φ in PHASES, t in 1:T],
        Ibs_r[b, φ, t] == sum(Ibr_r[k, φ, t] for k in out_br[b]; init = zero(AffExpr))
                        - sum(Ibr_r[k, φ, t] for k in in_br[b];  init = zero(AffExpr)))
    @constraint(model, [b in 1:nb, φ in PHASES, t in 1:T],
        Ibs_im[b, φ, t] == sum(Ibr_im[k, φ, t] for k in out_br[b]; init = zero(AffExpr))
                         - sum(Ibr_im[k, φ, t] for k in in_br[b];  init = zero(AffExpr)))

    # ---- power balance, eq. (37), linearised as eq. (38)–(39) --------------------------
    #  The only nonlinearity in the whole network model, and it lives at the buses rather
    #  than along the lines. Each product xy is replaced by x∘y + y∘x − x∘y∘.
    Plin(b, φ, t) = v_r_pr[b,φ,t]  * Ibs_r[b,φ,t]  + Ibs_r_pr[b,φ,t]  * v_r[b,φ,t] +
                    v_im_pr[b,φ,t] * Ibs_im[b,φ,t] + Ibs_im_pr[b,φ,t] * v_im[b,φ,t] -
                    v_r_pr[b,φ,t]  * Ibs_r_pr[b,φ,t] - v_im_pr[b,φ,t] * Ibs_im_pr[b,φ,t]
    Qlin(b, φ, t) = v_im_pr[b,φ,t] * Ibs_r[b,φ,t]  + Ibs_r_pr[b,φ,t]  * v_im[b,φ,t] -
                    v_r_pr[b,φ,t]  * Ibs_im[b,φ,t] - Ibs_im_pr[b,φ,t] * v_r[b,φ,t] -
                    v_im_pr[b,φ,t] * Ibs_r_pr[b,φ,t] + v_r_pr[b,φ,t]  * Ibs_im_pr[b,φ,t]

    @expression(model, netP[b = 1:nb, φ = PHASES, t = 1:T],
        (b == islack ? Pg[φ, t] : zero(AffExpr))
        + sum(Pdg[i, t] for i in pv_at[b, φ]; init = zero(AffExpr))
        - Pload_pk[b, φ] * Pmult[t])
    @expression(model, netQ[b = 1:nb, φ = PHASES, t = 1:T],
        (b == islack ? Qg[φ, t] : zero(AffExpr))
        + sum(Qdg[i, t] for i in pv_at[b, φ]; init = zero(AffExpr))
        - Qload_pk[b, φ] * Qmult[t])

    @constraint(model, [b in 1:nb, φ in PHASES, t in 1:T], netP[b,φ,t] == Plin(b,φ,t))
    @constraint(model, [b in 1:nb, φ in PHASES, t in 1:T], netQ[b,φ,t] == Qlin(b,φ,t))
    #                                                      ↑ where the droop meets the network

    # ---- voltage magnitude, eq. (40) --------------------------------------------------
    #  This v is the single quantity the droop module reads.
    @constraint(model, [b in 1:nb, φ in PHASES, t in 1:T],
        v[b, φ, t] == (v_r_pr[b,φ,t]  / hypot(v_r_pr[b,φ,t], v_im_pr[b,φ,t])) * v_r[b,φ,t]
                    + (v_im_pr[b,φ,t] / hypot(v_r_pr[b,φ,t], v_im_pr[b,φ,t])) * v_im[b,φ,t])

    # ---- thermal line limit, eq. (42), as a linear polygon inscribing the circle -------
    #  IVACOPF carries the line current as a decision variable, so (42) is a constraint
    #  you can simply write. LinDist3Flow has no I to write it about. Off by default:
    #  see IMAX_SEG above.
    if IMAX_SEG > 0
        for l in 1:IMAX_SEG
            θ = l * π / IMAX_SEG
            @constraint(model, [k in 1:nbr, φ in PHASES, t in 1:T],
                cos(θ) * Ibr_r[k,φ,t] + sin(θ) * Ibr_im[k,φ,t] <=  BR[k].imax)
            @constraint(model, [k in 1:nbr, φ in PHASES, t in 1:T],
                cos(θ) * Ibr_r[k,φ,t] + sin(θ) * Ibr_im[k,φ,t] >= -BR[k].imax)
        end
    end

    # ==================== DROOP BLOCK : Lambda / SOS2, unchanged from single phase ======
    #  v at the inverter's own bus and phase, and its reactive output, share one set of
    #  weights λ, so the operating point cannot leave the curve. The binaries z force the
    #  two nonzero weights to be adjacent (the SOS2 condition), which keeps the point on a
    #  segment rather than anywhere in the convex hull of the breakpoints.
    @variable(model, λ[1:6, 1:npv, 1:T] >= 0)
    @variable(model, z[1:5, 1:npv, 1:T], Bin)

    @constraint(model, [i in 1:npv, t in 1:T], sum(λ[j, i, t] for j in 1:6) == 1)
    @constraint(model, [i in 1:npv, t in 1:T], sum(z[j, i, t] for j in 1:5) == 1)
    @constraint(model, [i in 1:npv, t in 1:T], λ[1, i, t] <= z[1, i, t])
    @constraint(model, [j in 2:5, i in 1:npv, t in 1:T], λ[j, i, t] <= z[j-1, i, t] + z[j, i, t])
    @constraint(model, [i in 1:npv, t in 1:T], λ[6, i, t] <= z[5, i, t])

    @constraint(model, [i in 1:npv, t in 1:T],
        v[PV[i].bus, PV[i].phase, t] == sum(λ[j, i, t] * VBP[j] for j in 1:6))
    @constraint(model, [i in 1:npv, t in 1:T],
        Qdg[i, t] == sum(λ[j, i, t] * QSHAPE[j] * PV[i].Smax for j in 1:6))
    # ============================ END DROOP BLOCK ======================================

    # ---- inverter capability: 16-segment polygon inscribing the apparent-power circle --
    let k = 16
        for l in 1:k
            θ = l * π / k
            @constraint(model, [i in 1:npv, t in 1:T],
                cos(θ) * Pdg[i, t] + sin(θ) * Qdg[i, t] <=  PV[i].Smax)
            @constraint(model, [i in 1:npv, t in 1:T],
                cos(θ) * Pdg[i, t] + sin(θ) * Qdg[i, t] >= -PV[i].Smax)
        end
    end
    @constraint(model, [i in 1:npv, t in 1:T], Pdg[i, t] <= Pavail[i, t])

    # ---- curtailment and objective ------------------------------------------------------
    @constraint(model, [i in 1:npv, t in 1:T], PVC[i, t] == Pavail[i, t] - Pdg[i, t])
    @objective(model, Min, sum(PVC))

    nvar = num_variables(model)
    nbin = count(is_binary, all_variables(model))
    ncon = num_constraints(model; count_variable_in_set_constraints = false)
    if it == 1
        @printf("\nmodel   : %d variables (%d binary), %d constraints per pass\n\n",
                nvar, nbin, ncon)
        println(" pass    solve (s)    objective       MAPB        MRPB         MVM     status")
    end

    secs = @elapsed optimize!(model)
    total_solve += secs
    status = termination_status(model)
    status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ||
        error("pass $it terminated with status $status")

    # comprehensions rather than `value.(...)` so every result is a plain Array: the
    # PHASES axis makes JuMP hand back a DenseAxisArray, which does not broadcast into one
    bus3(x)  = [value(x[b, φ, t]) for b in 1:nb,  φ in PHASES, t in 1:T]
    br3(x)   = [value(x[k, φ, t]) for k in 1:nbr, φ in PHASES, t in 1:T]
    V       = bus3(v)
    Vr, Vi  = bus3(v_r), bus3(v_im)
    Ibs_r_v, Ibs_im_v = bus3(Ibs_r), bus3(Ibs_im)
    Ibr_r_v, Ibr_im_v = br3(Ibr_r), br3(Ibr_im)
    netP_v,  netQ_v   = bus3(netP), bus3(netQ)
    Pdg_v, Qdg_v, PVC_v = value.(Pdg), value.(Qdg), value.(PVC)

    # ---- linearisation error, measured against the TRUE nonlinear relations -----------
    #  MAPB and MRPB: the exact bilinear v·I of (37) minus the linear (38)–(39), which the
    #  constraint has forced equal to the net injection.
    #  MVM: the exact magnitude √(v_r² + v_im²) minus the linearised (40).
    #  Converging on these rather than on the model's own residual is what makes the
    #  returned point a genuine power-flow solution and not a solution of the approximation.
    MAPB = maximum(abs(Vr[b,φ,t]*Ibs_r_v[b,φ,t] + Vi[b,φ,t]*Ibs_im_v[b,φ,t] - netP_v[b,φ,t])
                   for b in 1:nb, φ in PHASES, t in 1:T)
    MRPB = maximum(abs(Vi[b,φ,t]*Ibs_r_v[b,φ,t] - Vr[b,φ,t]*Ibs_im_v[b,φ,t] - netQ_v[b,φ,t])
                   for b in 1:nb, φ in PHASES, t in 1:T)
    MVM  = maximum(abs(hypot(Vr[b,φ,t], Vi[b,φ,t]) - V[b,φ,t])
                   for b in 1:nb, φ in PHASES, t in 1:T)

    push!(iterlog, (iter = it, seconds = secs, objective = objective_value(model),
                    MAPB = MAPB, MRPB = MRPB, MVM = MVM, status = string(status)))
    @printf("  %2d   %9.2f   %10.6f   %.3e   %.3e   %.3e   %s\n",
            it, secs, objective_value(model), MAPB, MRPB, MVM, status)
    flush(stdout)

    # ---- refresh the linearisation point ----------------------------------------------
    v_r_pr .= Vr; v_im_pr .= Vi; Ibs_r_pr .= Ibs_r_v; Ibs_im_pr .= Ibs_im_v

    if max(MAPB, MRPB, MVM) < TOL
        converged = true
        break
    end
end
converged || @warn "stopped at MAX_ITER = $MAX_ITER without clearing TOL = $TOL"

# ==================================================================== 6) results ========
kWh(x) = x * SBASE / 1e3 / 4                      # p.u. summed over 15-min steps → kWh
E_avail = kWh(sum(Pavail))
E_curt  = kWh(sum(PVC_v))

# network losses, computed post hoc from the converged complex state: Σ ΔV·conj(I).
# LinDist3Flow cannot report this at all: it drops the loss term to become linear.
branch_loss(k, φ, t) =
    complex(Vr[bus_id[BR[k].from], φ, t] - Vr[bus_id[BR[k].to], φ, t],
            Vi[bus_id[BR[k].from], φ, t] - Vi[bus_id[BR[k].to], φ, t]) *
    conj(complex(Ibr_r_v[k, φ, t], Ibr_im_v[k, φ, t]))
Sloss    = sum(branch_loss(k, φ, t) for k in 1:nbr, φ in PHASES, t in 1:T)
loss_kWh = kWh(real(Sloss))

# how close the thermal limit (42) came to binding
Iutil = maximum(hypot(Ibr_r_v[k,φ,t], Ibr_im_v[k,φ,t]) / BR[k].imax
                for k in 1:nbr, φ in PHASES, t in 1:T)

println("\n================ RESULTS — IVACOPF (3-phase) + $METHOD ================")
@printf("passes              : %d %s (%.0e tolerance on MAPB / MRPB / MVM)\n",
        length(iterlog), converged ? "(converged)" : "(NOT converged)", TOL)
@printf("solve time          : %.1f s total, %.1f s last pass  (%s)\n",
        total_solve, secs, status)
WARMSTART && @printf("  + warm start      : %.1f s of exact sweeps\n", warm_seconds)
@printf("PV energy available : %.1f kWh\n", E_avail)
@printf("PV energy delivered : %.1f kWh\n", E_avail - E_curt)
@printf("PV curtailment      : %.2f kWh  (%.3f %%)\n", E_curt, 100 * E_curt / E_avail)
@printf("network losses      : %.2f kWh  (modelled, not neglected)\n", loss_kWh)
@printf("voltage range        : %.4f – %.4f p.u.\n", minimum(V), maximum(V))
for φ in PHASES
    @printf("  phase %d            : %.4f – %.4f p.u.\n", φ,
            minimum(V[:, φ, :]), maximum(V[:, φ, :]))
end
@printf("peak line loading   : %.1f %% of i_max  (thermal limit (42) %s)\n",
        100 * Iutil, IMAX_SEG > 0 ? "enforced" : "reported only")

# ---- verification 1: every dispatch point must lie ON the droop curve ------------------
droop_q(vv, qb) = vv <= VBP[2] ? qb :
                  vv <= VBP[3] ? qb * (VBP[3] - vv) / (VBP[3] - VBP[2]) :
                  vv <= VBP[4] ? 0.0 :
                  vv <= VBP[5] ? -qb * (vv - VBP[4]) / (VBP[5] - VBP[4]) : -qb
dev = maximum(abs(Qdg_v[i, t] - droop_q(V[PV[i].bus, PV[i].phase, t], PV[i].Smax))
              for i in 1:npv, t in 1:T)
@printf("\nmax |q_dispatch − q_curve| : %.3e p.u.   (exactness of the %s encoding)\n", dev, METHOD)

# ---- verification 2: IVACOPF voltages vs an exact three-phase power flow ---------------
#  Same backward/forward sweep used for the warm start, now re-run on the *solved*
#  dispatch. This is the number that separates the two hosts: it asks whether the voltage
#  the model told each inverter to read is the voltage it would really see.
sweep(t) = (Vc = sweep_state(Pdg_v, Qdg_v, t)[1]; [abs(Vc[b][φ]) for b in 1:nb, φ in PHASES])

tmax  = argmax([sum(Pdg_v[:, t]) for t in 1:T])            # busiest PV step
Vtrue = sweep(tmax)
gap   = maximum(abs.(V[:, :, tmax] .- Vtrue))
@printf("IVACOPF vs exact AC at t=%d : max |Δv| = %.3e p.u.  (true range %.4f – %.4f)\n",
        tmax, gap, minimum(Vtrue), maximum(Vtrue))
dev_true = maximum(abs(Qdg_v[i, tmax] - droop_q(Vtrue[PV[i].bus, PV[i].phase], PV[i].Smax))
                   for i in 1:npv)
@printf("droop residual at the TRUE voltages : %.3e p.u.\n", dev_true)

println("\npass-by-pass")
println("  pass   solve (s)     objective       MAPB        MRPB         MVM")
for r in iterlog
    @printf("  %3d  %9.2f   %11.6f   %.3e   %.3e   %.3e\n",
            r.iter, r.seconds, r.objective, r.MAPB, r.MRPB, r.MVM)
end

if INTERACTIVE
    println()
    println("============================================================")
    println(" RESULT VISUALIZATION")
    println("============================================================")
    println()
    println("The optimization is complete. The first figure overlays every")
    println("optimized inverter dispatch point on its corresponding Volt-VAr curve.")
    pause_tutorial("Press ENTER to display the full Volt-VAr dispatch plot...")
end

# ==================================================================== 7) figures ========
gr(size = (900, 520), legend = :topright, framestyle = :box, grid = true, gridalpha = 0.15,
   left_margin = 6Plots.mm, bottom_margin = 4Plots.mm)
hours = range(0, 24 - 24 / T, length = T)

#  One curve per inverter class, drawn in ABSOLUTE p.u. VArs so the four saturation
#  levels q̄ are visible; normalising by q̄ would collapse them onto a single line and
#  hide the point. Every "+" is one 15-min dispatch point and must sit on the curve of
#  its own class.
CLSCOL = [:seagreen, :orangered, :dodgerblue, :mediumorchid]
qmax   = maximum(g.Smax for g in PV)

p1 = plot(size = (1050, 780), grid = false, framestyle = :axes,
          title = "Dispatch vs. the IEEE 1547 droop — IVACOPF, $METHOD, four classes, three phases",
          titlefontsize = 12,
          xlabel = "Voltage (p.u.)", ylabel = "VAR Gen. (p.u.)",
          xlims = (VBP[1], VBP[6]), ylims = (-1.15qmax, 1.15qmax),
          xticks = 0.90:0.05:1.10, guidefontsize = 12, tickfontsize = 10,
          legendfontsize = 9, legend = :outertop, legend_columns = 4,
          foreground_color_legend = :black, background_color_legend = :white,
          left_margin = 6Plots.mm, right_margin = 8Plots.mm, bottom_margin = 5Plots.mm)

vspan!(p1, [VLIM[1], VLIM[2]], color = :lightblue, alpha = 0.30, lw = 0, label = false)
vline!(p1, [VLIM[1], VLIM[2]], ls = :dash, lw = 1.5, color = :gray65, label = false)
hline!(p1, [0.0], ls = :dash, lw = 1.5, color = :gray65, label = false)

for ci in 1:ncls
    qb = S_OVER_P * PV_CLASSES[ci][2] * 1e3 / SBASE
    plot!(p1, VBP, QSHAPE .* qb, lw = 3, color = CLSCOL[ci], label = false)
    scatter!(p1, VBP[2:5], (QSHAPE .* qb)[2:5], m = :circle, ms = 6,
             mc = CLSCOL[ci], msc = CLSCOL[ci], label = false)
    idx = [i for i in 1:npv if PV[i].cls == ci]
    scatter!(p1, vec([V[PV[i].bus, PV[i].phase, t] for i in idx, t in 1:T]),
             vec([Qdg_v[i, t] for i in idx, t in 1:T]),
             m = :+, ms = 7, msw = 2.5, mc = CLSCOL[ci], msc = CLSCOL[ci], label = false)
end

# legend proxies, parked outside xlims so they appear only in the key
for ci in 1:ncls
    plot!(p1, [1.5, 1.6], [0.0, 0.0], lw = 2, color = CLSCOL[ci], m = :circle, ms = 5,
          mc = CLSCOL[ci], msc = CLSCOL[ci], label = "$(PV_CLASSES[ci][1]) droop")
end
for ci in 1:ncls
    scatter!(p1, [1.5], [0.0], m = :+, ms = 7, msw = 2,
             mc = CLSCOL[ci], msc = CLSCOL[ci], label = "$(PV_CLASSES[ci][1]) Q")
end
vspan!(p1, [1.5, 1.6], color = :lightblue, alpha = 0.30, lw = 0,
       label = "Feasible Operation Region")
display(p1)
savefig(p1, joinpath(@__DIR__, "droop_dispatch_3ph$FIGSUF.png"))

# --------------------------------------------------------------------- interactive result explorer
if INTERACTIVE
    pause_tutorial("Press ENTER to explore one inverter class...")

    class_options = sort(unique(g.cls for g in PV))
    println("Available inverter classes in this simulation:")
    for (menu_idx, ci) in enumerate(class_options)
        nsites = count(g -> g.cls == ci, PV)
        println("  $menu_idx. $(PV_CLASSES[ci][1])  ($nsites sites)")
    end

    class_choice = ask_int(
        "Select an inverter class",
        1;
        allowed = 1:length(class_options),
    )
    ci = class_options[class_choice]
    idx = [i for i in 1:npv if PV[i].cls == ci]
    qb = S_OVER_P * PV_CLASSES[ci][2] * 1e3 / SBASE

    pclass = plot(
        VBP,
        QSHAPE .* qb,
        lw = 3,
        label = "$(PV_CLASSES[ci][1]) Volt-VAr curve",
        xlabel = "Voltage at inverter terminal (p.u.)",
        ylabel = "VAr output (p.u.)",
        title = "Dispatch of $(PV_CLASSES[ci][1]) inverters",
        xlims = (VBP[1], VBP[6]),
        ylims = (-1.15 * qb, 1.15 * qb),
        grid = false,
        framestyle = :axes,
        legend = :topright,
    )
    vspan!(pclass, [VLIM[1], VLIM[2]], alpha = 0.20, label = "Allowed voltage range")
    hline!(pclass, [0.0], ls = :dash, label = false)
    scatter!(
        pclass,
        vec([V[PV[i].bus, PV[i].phase, t] for i in idx, t in 1:T]),
        vec([Qdg_v[i, t] for i in idx, t in 1:T]),
        marker = :+,
        ms = 8,
        msw = 2,
        label = "Optimal dispatch",
    )
    display(pclass)
    savefig(pclass, joinpath(@__DIR__, "droop_selected_class_3ph$FIGSUF.png"))

    pause_tutorial("Press ENTER to inspect one inverter and one hour...")

    println("Inverters in $(PV_CLASSES[ci][1]):")
    for (menu_idx, i) in enumerate(idx)
        g = PV[i]
        @printf("  %d. bus %-6s  phase %d\n", menu_idx, BUSES[g.bus], g.phase)
    end

    inv_choice = ask_int(
        "Select an inverter",
        1;
        allowed = 1:length(idx),
    )
    i = idx[inv_choice]
    g = PV[i]

    hour = ask_int(
        "Choose an hour of day to inspect",
        12;
        allowed = 0:23,
    )
    tinspect = min(hour * T ÷ 24 + 1, T)
    vv = V[g.bus, g.phase, tinspect]
    qq = Qdg_v[i, tinspect]
    qcurve = droop_q(vv, g.Smax)

    println()
    println("Selected operating point")
    println("------------------------")
    @printf("Bus              : %s\n", BUSES[g.bus])
    @printf("Phase            : %d\n", g.phase)
    @printf("Inverter class   : %s\n", PV_CLASSES[g.cls][1])
    @printf("Displayed time   : %.2f h\n", hours[tinspect])
    @printf("Voltage          : %.6f p.u.\n", vv)
    @printf("Reactive output  : %.6f p.u.\n", qq)
    @printf("Droop prediction : %.6f p.u.\n", qcurve)
    @printf("|q - q_curve|    : %.3e p.u.\n", abs(qq - qcurve))

    pinspect = plot(
        VBP,
        QSHAPE .* g.Smax,
        lw = 3,
        xlabel = "Voltage at inverter terminal (p.u.)",
        ylabel = "VAr output (p.u.)",
        title = "$(PV_CLASSES[g.cls][1]) — bus $(BUSES[g.bus]), phase $(g.phase), $(round(hours[tinspect]; digits=2)) h",
        label = "Volt-VAr curve",
        xlims = (VBP[1], VBP[6]),
        ylims = (-1.15 * g.Smax, 1.15 * g.Smax),
        grid = false,
        framestyle = :axes,
        legend = :topright,
    )
    vspan!(pinspect, [VLIM[1], VLIM[2]], alpha = 0.20, label = "Allowed voltage range")
    hline!(pinspect, [0.0], ls = :dash, label = false)
    scatter!(pinspect, [vv], [qq], ms = 9, marker = :circle, label = "Selected operating point")
    display(pinspect)
    savefig(pinspect, joinpath(@__DIR__, "droop_selected_point_3ph$FIGSUF.png"))

    pause_tutorial("Press ENTER to display the feeder voltage envelope...")
end

p2 = plot(xlabel = "hour of day", ylabel = "voltage (p.u.)", xticks = 0:3:24, xlims = (0, 24),
          title = "Feeder voltage envelope by phase — IVACOPF")
for (φ, c) in zip(PHASES, (:seagreen, :orangered, :dodgerblue))
    plot!(p2, hours, [maximum(V[:, φ, t]) for t in 1:T], lw = 2, color = c, label = "phase $φ max")
    plot!(p2, hours, [minimum(V[:, φ, t]) for t in 1:T], lw = 2, ls = :dash, color = c,
          label = "phase $φ min")
end
hline!(p2, [VLIM[1], VLIM[2]], ls = :dot, lw = 1.5, color = :red, label = "limits")
display(p2)
savefig(p2, joinpath(@__DIR__, "voltage_envelope_3ph$FIGSUF.png"))

INTERACTIVE && pause_tutorial("Press ENTER to display available versus delivered PV power...")

p3 = plot(hours, [sum(Pavail[:, t]) * SBASE / 1e3 for t in 1:T], lw = 2, ls = :dash,
          color = :grey45, label = "available",
          xlabel = "hour of day", ylabel = "kW", xticks = 0:3:24, xlims = (0, 24),
          title = "Fleet PV: available vs delivered — IVACOPF")
plot!(p3, hours, [sum(Pdg_v[:, t]) * SBASE / 1e3 for t in 1:T], lw = 2, color = :darkorange2,
      fillrange = 0, fillalpha = 0.15, label = "delivered")
display(p3)
savefig(p3, joinpath(@__DIR__, "pv_dispatch_3ph$FIGSUF.png"))

println("\nwrote droop_dispatch_3ph$FIGSUF.png, voltage_envelope_3ph$FIGSUF.png, " *
        "pv_dispatch_3ph$FIGSUF.png")
if INTERACTIVE
    println("also wrote droop_selected_class_3ph$FIGSUF.png and droop_selected_point_3ph$FIGSUF.png")
    println()
    println("Tutorial complete.")
end
