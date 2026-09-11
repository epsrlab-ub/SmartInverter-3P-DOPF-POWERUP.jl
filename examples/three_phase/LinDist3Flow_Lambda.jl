# =====================================================================================
#  Three-phase Volt-VAr droop in a distribution OPF
#  HOST   : LinDist3Flow, the three-phase (multiphase) LinDistFlow linearisation
#                          of Sankur, Dobbe, Stewart, Callaway & Arnold,
#                          arXiv:1606.04492 (2016), their eqs. (20)–(23)
#  DROOP  : Lambda / SOS2 (convex combination of breakpoints, adjacency forced by binaries)
#
#  Test system : network_5_Feeder_2, a real 415/240 V LV feeder from the ENWL dataset,
#                Kron-reduced to three wires. 194 buses, 193 lines, 489 m, and 18
#                single-phase loads spread 4 / 5 / 9 across phases 1 / 2 / 3, so the
#                network is genuinely unbalanced, which is the point of modelling it
#                three-phase at all.
#  Horizon     : 24 h at 15-minute resolution (96 steps), residential load shape,
#                clear-sky irradiance.
#  Objective   : minimise total PV curtailment over the day.
#
#  This is the three-phase counterpart of the single-phase 33-bus / LinDistFlow / Lambda
#  case. The droop block is unchanged from the single-phase version: it only ever needs
#  a voltage variable per inverter, so everything new here is in the network model.
#
#  Equation numbers in the comments below are the tutorial's:
#  https://epsrlab-ub.github.io/SmartInverter-3P-DOPF.jl/dev/tutorial_voltvar/#LinDist3Flow:-the-linear-host
#
#  Run:   julia --project=. LinDist3Flow_BigM.jl
# =====================================================================================

using JuMP, Gurobi                    # MILP solver (Gurobi is free for academic use)
using JSON3, Printf, Plots

const PHASES = 1:3
const DATA   = joinpath(@__DIR__, "data")

# =====================================================================================
# Interactive tutorial configuration
# =====================================================================================
# By default this script is interactive. For unattended runs, set TP_INTERACTIVE=0 and
# use the existing environment variables TP_NPV and TP_STEPS.
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

function pause_tutorial(message::AbstractString = "Press ENTER to continue...")
    INTERACTIVE || return
    println()
    print(message)
    readline()
    println()
end

# Four inverter classes. Each phase receives one of each when N_PV_PER_PHASE = 4, so
# size is not confounded with phase or distance from the substation.
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
const FIGSUF         = "_lambda"

npv_default   = parse(Int, get(ENV, "TP_NPV", "4"))
steps_default = parse(Int, get(ENV, "TP_STEPS", "96"))

if INTERACTIVE
    println()
    println("============================================================")
    println(" THREE-PHASE SMART-INVERTER OPF — INTERACTIVE TUTORIAL")
    println("============================================================")
    println()
    println("Configure a LinDist3Flow simulation with IEEE 1547")
    println("Volt-VAr control encoded using the Lambda / SOS2 formulation.")
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

    println()
    println("---------------- Selected configuration ----------------")
    println("PV sites / phase : $N_PV_PER_PHASE")
    println("Time steps       : $T")
    println("Network model    : LinDist3Flow")
    println("Droop encoding   : $METHOD")
    println("--------------------------------------------------------")
    pause_tutorial("Press ENTER to load the feeder and build the simulation...")
else
    N_PV_PER_PHASE = npv_default
    T = steps_default
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
    push!(BR, Branch(String(l.bus_from), String(l.bus_to), linecode_Z(lc, Float64(l.length))))
end
nbr = length(BR)

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

# =============================================== 2) LinDist3Flow voltage-drop matrices ==
#  This whole block is eq. (31)–(32) of the tutorial, i.e. eqs. (21)–(23) of Sankur et al.
#
#  For a line with 3×3 phase impedance Z, and voltages that are near-balanced (the
#  standing assumption of every LinDistFlow variant), the drop in *squared* magnitude is
#
#       w_j = w_i − ( aR · P_ij + aX · Q_ij )
#
#  where, with the rotation operator α = e^{−j2π/3},
#
#       aR[φ,ψ] =  2·Re( α^{φ−ψ} · conj(Z[φ,ψ]) )
#       aX[φ,ψ] = −2·Im( α^{φ−ψ} · conj(Z[φ,ψ]) )
#
#  The ±√3 cross terms that appear in the published form of these matrices are just this
#  rotation written out. Two checks worth keeping in mind:
#    · single phase (α^0 = 1): aR = 2r, aX = 2x, recovering w_j = w_i − 2(rP + xQ);
#    · Z diagonal: aR, aX are diagonal and the phases decouple into three LinDistFlows.
#
#  We then work in magnitude rather than squared magnitude, exactly as the single-phase
#  code does: w = v² and v ≈ 1 give w_j − w_i ≈ 2·VNOM·(v_j − v_i), hence the 2·VNOM
#  divisor below.
#
#  Reference: M. D. Sankur, R. Dobbe, E. Stewart, D. S. Callaway, D. B. Arnold, "A
#  linearized power flow model for optimization in unbalanced distribution systems",
#  arXiv:1606.04492 (2016).  aR and aX below are MINUS that paper's M^P and M^Q, its
#  eqs. (22)-(23); the drop itself is its eq. (21) and the lossless power balance used
#  further down is its eq. (20). Both follow from its exact Dist3Flow eqs. (14)-(17)
#  under assumptions A1 (constant inter-phase voltage ratios) and A2 (constant losses).
#
#  Term-by-term check of the (a,b) entry, paper against this code:
#      M^P(a,b) =  r_ab - sqrt(3) x_ab   ->   aR[1,2] = -r_ab + sqrt(3) x_ab
#      M^Q(a,b) =  x_ab + sqrt(3) r_ab   ->   aX[1,2] = -x_ab - sqrt(3) r_ab
const ALPHA = exp(-2π * im / 3)

# aR, aX of eq. (32): minus the M^P, M^Q of Sankur et al. (22)–(23). The trailing
# 2·VNOM divisor converts to the magnitude form the drop constraint of eq. (33) uses.
function drop_matrices(Z::Matrix{ComplexF64})
    aR = zeros(3, 3)
    aX = zeros(3, 3)
    for p in PHASES, q in PHASES
        c = ALPHA^(p - q) * conj(Z[p, q])
        aR[p, q] =  2 * real(c)
        aX[p, q] = -2 * imag(c)
    end
    return aR ./ (2 * VNOM), aX ./ (2 * VNOM)
end

AR = [drop_matrices(b.Z)[1] for b in BR]
AX = [drop_matrices(b.Z)[2] for b in BR]

# ============================================================ 3) place the PV inverters ==
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

# irradiance ceiling per inverter per step
Pavail = [PV[i].Pmax * G[t] for i in 1:npv, t in 1:T]

@printf("network : %d buses, %d lines, %d single-phase loads\n", nb, nbr, length(net.load))
@printf("bases   : V = %.2f V (phase-neutral), S = %.0f kVA/phase, Z = %.4f Ω\n",
        VBASE, SBASE_KVA, ZBASE)
@printf("load    : %.2f kW peak total\n", sum(Pload_pk) * SBASE / 1e3)
totkW = sum(g.Pmax for g in PV) * SBASE / 1e3
@printf("PV      : %d smart inverters in %d classes, %.1f kW total (%.0f%% of peak load)\n",
        npv, ncls, totkW, 100 * totkW / (sum(Pload_pk) * SBASE / 1e3))

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

# ================================================================= 4) build and solve ===
if INTERACTIVE
    println()
    println("============================================================")
    println(" SOLVING THE THREE-PHASE LinDist3Flow OPF")
    println("============================================================")
    println()
    println("The optimization will:")
    println("  1. enforce the multiphase LinDist3Flow voltage-drop equations,")
    println("  2. enforce phase-by-phase active and reactive power balance,")
    println("  3. embed the IEEE 1547 Volt-VAr curve using the Lambda / SOS2 formulation, and")
    println("  4. minimize total PV curtailment over the selected 24-hour horizon.")
    println()
    pause_tutorial("Press ENTER to start the optimization...")
end
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "MIPGap", MIP_GAP)
set_optimizer_attribute(model, "OutputFlag", 0)

@variable(model, VLIM[1] <= v[1:nb, PHASES, 1:T] <= VLIM[2], start = VNOM)  # limits: eq. (33)
@variable(model, P[1:nbr, PHASES, 1:T])
@variable(model, Q[1:nbr, PHASES, 1:T])
@variable(model, Pg[PHASES, 1:T])
@variable(model, Qg[PHASES, 1:T])
@variable(model, 0 <= Pdg[1:npv, 1:T])
@variable(model, Qdg[1:npv, 1:T])
@variable(model, 0 <= PVC[1:npv, 1:T])

islack = bus_id[SLACK]
# slack reference: first line of eq. (33)
@constraint(model, [φ in PHASES, t in 1:T], v[islack, φ, t] == VNOM)

# ---- voltage drop: eq. (31)–(32) here, eq. (21) of Sankur et al. ---------------------
@constraint(model, [k in 1:nbr, φ in PHASES, t in 1:T],
    v[bus_id[BR[k].to], φ, t] == v[bus_id[BR[k].from], φ, t]
        - sum(AR[k][φ, ψ] * P[k, ψ, t] + AX[k][φ, ψ] * Q[k, ψ, t] for ψ in PHASES))

# ---- power balance per bus and phase: eq. (33) here, eq. (20) of Sankur et al.;
#      losses neglected, which is the whole of the linearisation ------------------------
out_br = [Int[] for _ in 1:nb]
in_br  = [Int[] for _ in 1:nb]
for (k, br) in enumerate(BR)
    push!(out_br[bus_id[br.from]], k)
    push!(in_br[bus_id[br.to]], k)
end
pv_at = [Int[] for _ in 1:nb, _ in 1:3]
for (i, g) in enumerate(PV); push!(pv_at[g.bus, g.phase], i); end

@expression(model, netP[b = 1:nb, φ = PHASES, t = 1:T],
    (b == islack ? Pg[φ, t] : zero(AffExpr))
    + sum(Pdg[i, t] for i in pv_at[b, φ]; init = zero(AffExpr))
    - Pload_pk[b, φ] * Pmult[t])
@expression(model, netQ[b = 1:nb, φ = PHASES, t = 1:T],
    (b == islack ? Qg[φ, t] : zero(AffExpr))
    + sum(Qdg[i, t] for i in pv_at[b, φ]; init = zero(AffExpr))
    - Qload_pk[b, φ] * Qmult[t])

@constraint(model, [b in 1:nb, φ in PHASES, t in 1:T],
    netP[b, φ, t] == sum(P[k, φ, t] for k in out_br[b]; init = zero(AffExpr))
                   - sum(P[k, φ, t] for k in in_br[b];  init = zero(AffExpr)))
@constraint(model, [b in 1:nb, φ in PHASES, t in 1:T],
    netQ[b, φ, t] == sum(Q[k, φ, t] for k in out_br[b]; init = zero(AffExpr))
                   - sum(Q[k, φ, t] for k in in_br[b];  init = zero(AffExpr)))

# ======================== DROOP BLOCK : Lambda / SOS2, unchanged from single phase ======
#  The IEEE 1547 curve of eq. (1), entering the host through eq. (30).
#  v at the inverter's own bus and phase, and its reactive output, share one set of
#  weights λ, so the operating point cannot leave the curve. The binaries z force the two
#  nonzero weights to be adjacent (the SOS2 condition), which keeps the point on a
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
# ================================ END DROOP BLOCK ======================================

# ---- inverter capability, eq. (30): 16-segment polygon inscribing the S-circle ---------
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

# ---- curtailment and objective, eq. (30) -----------------------------------------------
@constraint(model, [i in 1:npv, t in 1:T], PVC[i, t] == Pavail[i, t] - Pdg[i, t])
@objective(model, Min, sum(PVC))

@printf("model   : %d variables (%d binary), %d constraints\n",
        num_variables(model), count(is_binary, all_variables(model)),
        num_constraints(model; count_variable_in_set_constraints = false))

secs = @elapsed optimize!(model)
status = termination_status(model)
status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || error("solver returned $status")

V   = value.(v)
Pdg_v, Qdg_v, PVC_v = value.(Pdg), value.(Qdg), value.(PVC)

# ==================================================================== 5) results ========
Δt_hours = 24 / T
kWh(x) = x * SBASE / 1e3 * Δt_hours              # p.u. summed over selected time steps → kWh
E_avail = kWh(sum(Pavail))
E_curt  = kWh(sum(PVC_v))
println("\n================ RESULTS — LinDist3Flow + $METHOD ================")
@printf("solve time          : %.1f s  (%s)\n", secs, status)
@printf("PV energy available : %.1f kWh\n", E_avail)
@printf("PV energy delivered : %.1f kWh\n", E_avail - E_curt)
@printf("PV curtailment      : %.2f kWh  (%.3f %%)\n", E_curt, 100 * E_curt / E_avail)
@printf("voltage range        : %.4f – %.4f p.u.\n", minimum(V), maximum(V))
for φ in PHASES
    @printf("  phase %d            : %.4f – %.4f p.u.\n", φ,
            minimum(V[:, φ, :]), maximum(V[:, φ, :]))
end

# ---- verification 1: the droop-deviation metric Δ of eq. (29) --------------------------
#      every dispatch point must lie ON the droop curve of its own inverter
droop_q(vv, qb) = vv <= VBP[2] ? qb :
                  vv <= VBP[3] ? qb * (VBP[3] - vv) / (VBP[3] - VBP[2]) :
                  vv <= VBP[4] ? 0.0 :
                  vv <= VBP[5] ? -qb * (vv - VBP[4]) / (VBP[5] - VBP[4]) : -qb
dev = maximum(abs(Qdg_v[i, t] - droop_q(V[PV[i].bus, PV[i].phase, t], PV[i].Smax))
              for i in 1:npv, t in 1:T)

# ==================================================================== 6) figures ========
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

gr(size = (900, 520), legend = :topright, framestyle = :box, grid = true, gridalpha = 0.15,
   left_margin = 6Plots.mm, bottom_margin = 4Plots.mm)
hours = range(0, 24 - 24 / T, length = T)

CLSCOL = [:seagreen, :orangered, :dodgerblue, :mediumorchid]
qmax   = maximum(g.Smax for g in PV)

# --------------------------------------------------------------------- full fleet dispatch
p1 = plot(size = (1050, 780), grid = false, framestyle = :axes,
          title = "Dispatch vs. the IEEE 1547 droop — LinDist3Flow, $METHOD, four classes, three phases",
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

gui(p1)
println("\nPress ENTER when you are finished viewing the plot.")
readline()

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
        title = "Dispatch of $(PV_CLASSES[ci][1]) inverters — LinDist3Flow",
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

    gui(pclass)
    println("\nPress ENTER when you are finished viewing the plot.")
    readline()

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

    println()
    println("Selected operating point")
    println("------------------------")
    @printf("Bus              : %s\n", BUSES[g.bus])
    @printf("Phase            : %d\n", g.phase)
    @printf("Inverter class   : %s\n", PV_CLASSES[g.cls][1])
    @printf("Displayed time   : %.2f h\n", hours[tinspect])
    @printf("Voltage          : %.6f p.u.\n", vv)
    @printf("Reactive output  : %.6f p.u.\n", qq)

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

    gui(pinspect)
    println("\nPress ENTER when you are finished viewing the plot.")
    readline()

    pause_tutorial("Press ENTER to display the feeder voltage envelope...")
end

# --------------------------------------------------------------------- feeder voltage envelope
p2 = plot(xlabel = "hour of day", ylabel = "voltage (p.u.)", xticks = 0:3:24, xlims = (0, 24),
          title = "Feeder voltage envelope by phase — LinDist3Flow")
for (φ, c) in zip(PHASES, (:seagreen, :orangered, :dodgerblue))
    plot!(p2, hours, [maximum(V[:, φ, t]) for t in 1:T], lw = 2, color = c, label = "phase $φ max")
    plot!(p2, hours, [minimum(V[:, φ, t]) for t in 1:T], lw = 2, ls = :dash, color = c,
          label = "phase $φ min")
end
hline!(p2, [VLIM[1], VLIM[2]], ls = :dot, lw = 1.5, color = :red, label = "limits")

gui(p2)
println("\nPress ENTER when you are finished viewing the plot.")
readline()

INTERACTIVE && pause_tutorial("Press ENTER to display available versus delivered PV power...")

# --------------------------------------------------------------------- available vs delivered PV
p3 = plot(hours, [sum(Pavail[:, t]) * SBASE / 1e3 for t in 1:T], lw = 2, ls = :dash,
          color = :grey45, label = "available",
          xlabel = "hour of day", ylabel = "kW", xticks = 0:3:24, xlims = (0, 24),
          title = "Fleet PV: available vs delivered — LinDist3Flow")
plot!(p3, hours, [sum(Pdg_v[:, t]) * SBASE / 1e3 for t in 1:T], lw = 2, color = :darkorange2,
      fillrange = 0, fillalpha = 0.15, label = "delivered")

gui(p3)
println("\nPress ENTER when you are finished viewing the plot.")
readline()


# ---- interactive verification: LinDist3Flow vs exact three-phase power flow ------------

if INTERACTIVE
    println()
    println("============================================================")
    println(" LINDIST3FLOW ACCURACY CHECK")
    println("============================================================")
    println()
    println("LinDist3Flow is a linear approximation of the three-phase feeder.")
    println("We can compare its voltage prediction against an exact")
    println("three-phase backward/forward sweep using the same optimized")
    println("PV active- and reactive-power dispatch.")
    println()

    # Automatically identify the time of maximum fleet PV output
    tmax_auto = argmax([sum(Pdg_v[:, t]) for t in 1:T])
    hour_auto = (tmax_auto - 1) * 24 / T

    @printf("The highest-PV operating point occurs at approximately %.2f h.\n", hour_auto)
    println()
    println("You may inspect that operating point or choose another hour.")

    hour_check = ask_int(
        "Choose an hour of day for the exact AC comparison",
        round(Int, hour_auto);
        allowed = 0:23,
    )

    tcheck = min(hour_check * T ÷ 24 + 1, T)

else
    # For unattended runs, retain the original behavior:
    # validate at the time of maximum PV production.
    tcheck = argmax([sum(Pdg_v[:, t]) for t in 1:T])
end



# ---- verification 2: LinDist3Flow voltages vs an exact three-phase power flow ----------
#  Backward/forward sweep on the same dispatch. Kron-reduced three-wire model with
#  implicit ground and zero shunts, so the sweep is exact for these injections.
function sweep(t)
    # Start from a properly rotated three-phase set: 1∠0°, 1∠−120°, 1∠+120°. Seeding all
    # three phases at 1∠0° instead is a silent and expensive mistake: the mutual terms
    # then add rather than largely cancelling, and the sweep reports roughly twice the
    # true voltage deviation.
    V0  = ComplexF64[1, exp(-2π * im / 3), exp(2π * im / 3)]
    Vc  = [copy(V0) for _ in 1:nb]
    Ibr = [zeros(ComplexF64, 3) for _ in 1:nbr]
    child_of = [Int[] for _ in 1:nb]
    for (k, br) in enumerate(BR); push!(child_of[bus_id[br.from]], k); end
    rev = reverse(ORDER)
    for _ in 1:60
        S = [zeros(ComplexF64, 3) for _ in 1:nb]
        for b in 1:nb, φ in PHASES
            S[b][φ] = complex(Pload_pk[b, φ] * Pmult[t], Qload_pk[b, φ] * Qmult[t])
        end
        for (i, g) in enumerate(PV)
            S[g.bus][g.phase] -= complex(Pdg_v[i, t], Qdg_v[i, t])
        end
        Ibus = [conj.(S[b] ./ Vc[b]) for b in 1:nb]
        for bname in rev                                   # backward: currents
            b = bus_id[bname]
            b == islack && continue
            (_, k) = PARENT[bname]
            Ibr[k] = Ibus[b] + sum((Ibr[c] for c in child_of[b]), init = zeros(ComplexF64, 3))
        end
        for bname in ORDER                                 # forward: voltages
            bname == SLACK && continue
            (par, k) = PARENT[bname]
            Vc[bus_id[bname]] = Vc[bus_id[par]] - BR[k].Z * Ibr[k]
        end
    end
    return [abs(Vc[b][φ]) for b in 1:nb, φ in PHASES]
end


# Run exact three-phase backward/forward sweep
Vtrue = sweep(tcheck)

# Compare LinDist3Flow voltage magnitudes against exact AC voltages
Vapprox = V[:, :, tcheck]
Verror  = abs.(Vapprox .- Vtrue)
gap     = maximum(Verror)



println()
println("---------------- Voltage-model comparison ----------------")
@printf("Displayed time          : %.2f h\n", (tcheck - 1) * 24 / T)
@printf("LinDist3Flow range       : %.4f – %.4f p.u.\n",
        minimum(Vapprox), maximum(Vapprox))
@printf("Exact AC range           : %.4f – %.4f p.u.\n",
        minimum(Vtrue), maximum(Vtrue))
@printf("Maximum voltage error    : %.3e p.u.\n", gap)
println("----------------------------------------------------------")


if INTERACTIVE
    println()
    println("Tutorial complete.")
    pause_tutorial("Press ENTER when you are finished viewing the plots and want to exit...")
else
    println("Press ENTER when you are finished viewing the plots and want to exit...")
    readline()
end
