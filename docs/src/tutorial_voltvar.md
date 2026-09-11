# Modeling Smart Inverters in Three-Phase Distribution Optimal Power Flow

A smart inverter follows a Volt-VAr curve where it measures its own terminal voltage and autonomously decides how much reactive power to inject or absorb. A distribution optimal power flow (DOPF) problem that ignores the curve will return a reactive dispatch that could be unsuitable at the local inverter controller level.

This tutorial presents three ways to embed the Volt-VAr curve into a three-phase DOPF, so that every dispatch point the solver returns is one the inverter would actually produce. 

To keep the separation between *encoding* and *network model* measurable rather than merely asserted, the three encodings are run against **two** three-phase hosts, a linear one and a near-exact one. Both feeders used below are real Electricity North West low-voltage networks from the
*Low Voltage Network Solutions* project, Kron-reduced to three wires: `network_5_Feeder_2`
[[14]](#ref-14) for the case study and `network_17_Feeder_6` [[15]](#ref-15) for the
scalability check. The Kron reduction follows [[16]](#ref-16), and the conductor impedances
are those of [[17]](#ref-17).

```@setup tut
using JSON3, Plots, Printf, Markdown
gr(fmt = :svg, size = (760, 420), legendfontsize = 8, titlefontsize = 10,
   guidefontsize = 9, tickfontsize = 8, framestyle = :box, grid = true,
   gridalpha = 0.15, dpi = 150)

# Everything on this page is drawn from the committed three-phase results, so the
# documentation builds without an optimisation solver.
TPR  = joinpath("assets", "results", "threephase")
tpc  = JSON3.read(read(joinpath(TPR, "case.json"), String))
tpr  = Dict(m => JSON3.read(read(joinpath(TPR, "$m.json"), String))
            for m in ("bigm", "lambda", "heaviside"))
tpi  = Dict(m => JSON3.read(read(joinpath(TPR, "iva_$m.json"), String))
            for m in ("bigm", "lambda", "heaviside"))
tpsc = JSON3.read(read(joinpath(TPR, "scalability.json"), String))

const NAMES = Dict("bigm" => "Big-M", "lambda" => "Lambda / SOS2", "heaviside" => "Heaviside")
const ORDER = ["bigm", "lambda", "heaviside"]
const TPHOSTS = [("LinDist3Flow", tpr), ("IVACOPF", tpi)]

hours = range(0, 24 - 24/tpc.n_steps, length = tpc.n_steps)   # quarter-hourly steps
fmt(x, n) = @sprintf("%.*f", n, x)
sci(x)    = @sprintf("%.2e", x)
md(rows...) = Markdown.parse(join(rows, "\n"))

Vbp, qshape = collect(Float64, tpc.Vbp), collect(Float64, tpc.qshape)

# Table builders. These live here rather than in the visible blocks so that the page
# shows tables, not the string-mangling that produces them, while still deriving every
# number from the committed results rather than hard-coding it.
breakpoint_table() = md(
    "| | " * join(["``V^{\\text{bp}}_$i``" for i in 1:6], " | ") * " |",
    "|:--|" * repeat("--:|", 6),
    "| voltage (p.u.) | " * join(fmt.(Vbp, 2), " | ") * " |",
    "| ``q/\\bar q`` | " * join(fmt.(qshape, 0), " | ") * " |")

tp_class_table() = md(
    "| class | ``P`` rated | ``S_{\\max}`` | ``\\bar q`` (p.u.) | sites | buses |",
    "|:--|--:|--:|--:|--:|:--|",
    join([let c = tpc.classes[k]
              bs = [s for s in tpc.sites if s.class_idx == k]
              "| $(c.name) | $(fmt(c.P_kW, 0)) kW | $(fmt(c.S_kVA, 2)) kVA | " *
              "$(fmt(c.qbar_pu, 4)) | $(length(bs)) | " *
              join(["$(s.bus) (φ$(s.phase))" for s in bs], ", ") * " |"
          end for k in eachindex(tpc.classes)], "\n"))

tp_scale_table() = md(
    "| encoding | feeder | time steps | variables | binaries | solve (s) | max droop deviation |",
    "|:--|:--|--:|--:|--:|:--|--:|",
    join([let r = row
              solve = r.ok ? fmt(r.solve_seconds, 1) : "**did not solve**"
              dev   = r.max_droop_deviation === nothing ? "n/a" : sci(r.max_droop_deviation)
              "| $(NAMES[r.encoding]) | $(r.feeder) | $(r.steps) | " *
              "$(r.nvar) | $(r.nbin) | $solve | $dev |"
          end for row in tpsc.runs], "\n"))

# Every encoding on every host, side by side. `passes` is 1 for the linear host, which
# has no outer loop at all.
tp_host_table() = md(
    "| host | encoding | class | solver | variables | binaries | passes | solve (s) | " *
    "curtailed (kWh) | curtailed (%) | losses (kWh) | voltage range (p.u.) |",
    "|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|:--|:--|",
    join([let r = res[m]
              "| $hname | $(NAMES[m]) | $(r.model_class) | `$(r.solver)` | $(r.nvar) | " *
              "$(r.nbin) | $(get(r, :n_passes, 1)) | $(fmt(r.solve_seconds, 1)) | " *
              "$(fmt(r.E_curt_kWh, 2)) | $(fmt(r.curt_percent, 3)) | " *
              "$(haskey(r, :loss_kWh) ? fmt(r.loss_kWh, 2) : "not modelled") | " *
              "$(fmt(r.Vmin, 4)) – $(fmt(r.Vmax, 4)) |"
          end for (hname, res) in TPHOSTS for m in ORDER], "\n"))

# Exactness of the *encoding* inside each host: does the returned dispatch lie on the
# curve the model itself reports? This is a different question from the audit below.
tp_exact_table() = md(
    "| encoding | LinDist3Flow | IVACOPF |",
    "|:--|--:|--:|",
    join(["| $(NAMES[m]) | $(sci(tpr[m].max_droop_deviation)) | " *
          "$(sci(tpi[m].max_droop_deviation)) |" for m in ORDER], "\n"))

# The audit that separates the hosts: take each dispatch, solve the EXACT three-phase AC
# power flow for those injections, and ask what the inverters would really have seen.
tp_audit_table() = md(
    "| host | its own ``v`` vs the true AC ``v`` | droop residual at the **true** voltage | " *
    "bus-steps outside ``[0.95, 1.05]`` | true voltage range (p.u.) |",
    "|:--|--:|--:|--:|:--|",
    join([let r = res["lambda"]
              "| $hname | $(sci(r.audit.v_gap)) p.u. | $(sci(r.audit.droop_residual_true_v)) p.u. | " *
              "$(r.audit.n_limit_violations) | $(fmt(r.audit.true_Vmin, 4)) – " *
              "$(fmt(r.audit.true_Vmax, 4)) |"
          end for (hname, res) in TPHOSTS], "\n"))

# The successive-linearisation loop, pass by pass, measured against the error
# metrics MAPB / MRPB / MVM rather than against the model's internal residual.
tp_pass_table(m = "lambda") = md(
    "| pass | solve (s) | objective (p.u. curtailed) | MAPB | MRPB | MVM | solver status |",
    "|--:|--:|--:|--:|--:|--:|:--|",
    join(["| $(r.iter) | $(fmt(r.seconds, 1)) | $(fmt(r.objective, 6)) | $(sci(r.MAPB)) | " *
          "$(sci(r.MRPB)) | $(sci(r.MVM)) | `$(r.status)` |" for r in tpi[m].iterations], "\n"))

const TPCOL = [:seagreen, :orangered, :dodgerblue, :mediumorchid]

# The PV resource the fleet is working against, and how much of it survives the droop.
function tp_pv_figure(m = "lambda"; res = tpi)
    r = res[m]
    p = plot(hours, collect(Float64, r.P_avail_kW), lw = 2, ls = :dash, color = :grey45,
             label = "available", xlabel = "hour of day", ylabel = "kW",
             title = "PV across the twelve inverters: available and delivered",
             xticks = 0:3:24, xlims = (0, 24), legend = :topleft)
    plot!(p, hours, collect(Float64, r.P_disp_kW), lw = 2.4, color = :darkorange2,
          fillrange = 0, fillalpha = 0.18, label = "delivered")
    p
end

function tp_droop_figure(m = "lambda"; res = tpr, host = "LinDist3Flow")
    r    = res[m]
    qmax = maximum(c.qbar_pu for c in tpc.classes)
    p = plot(size = (860, 620), grid = false, framestyle = :axes,
             title = "Three-phase dispatch vs. the droop — $host, $(r.method)",
             titlefontsize = 11,
             xlabel = "voltage at the inverter terminal (p.u.)", ylabel = "VAr output (p.u.)",
             xlims = (Vbp[1], Vbp[6]), ylims = (-1.15qmax, 1.15qmax),
             xticks = 0.90:0.05:1.10, legend = :outertop, legend_columns = 4,
             legendfontsize = 8, foreground_color_legend = :black,
             background_color_legend = :white, left_margin = 4Plots.mm)
    vspan!(p, [tpc.Vmin_limit, tpc.Vmax_limit], color = :lightblue, alpha = 0.30,
           lw = 0, label = false)
    hline!(p, [0.0], ls = :dash, lw = 1.2, color = :gray65, label = false)
    for k in eachindex(tpc.classes)
        qb  = tpc.classes[k].qbar_pu
        idx = [i for i in eachindex(tpc.sites) if tpc.sites[i].class_idx == k]
        plot!(p, Vbp, qshape .* qb, lw = 2.5, color = TPCOL[k], label = false)
        scatter!(p, vcat([collect(Float64, r.Vdg_series[i]) for i in idx]...),
                 vcat([collect(Float64, r.Qdg_series[i]) for i in idx]...),
                 m = :+, ms = 5, msw = 2, mc = TPCOL[k], msc = TPCOL[k], label = false)
        plot!(p, [1.5, 1.6], [0.0, 0.0], lw = 2, color = TPCOL[k], m = :circle, ms = 4,
              mc = TPCOL[k], msc = TPCOL[k], label = tpc.classes[k].name)
    end
    p
end

function tp_envelope_figure(m = "lambda"; res = tpr, host = "LinDist3Flow")
    r = res[m]
    p = plot(xlabel = "hour of day", ylabel = "voltage (p.u.)", xticks = 0:3:24,
             xlims = (0, 24), legend = :topright,
             title = "Voltage envelope by phase, $host — the three phases do not coincide")
    for (φ, c) in zip(1:3, (:seagreen, :orangered, :dodgerblue))
        plot!(p, hours, collect(Float64, r.Vmax_t[φ]), lw = 2, color = c, label = "phase $φ max")
        plot!(p, hours, collect(Float64, r.Vmin_t[φ]), lw = 2, ls = :dash, color = c,
              label = "phase $φ min")
    end
    hline!(p, [tpc.Vmin_limit, tpc.Vmax_limit], ls = :dot, lw = 1.5, color = :red,
           label = "limits")
    p
end

# The voltage envelopes of the two hosts on one axis: same feeder, same dispatch problem, and
# a visible offset that is entirely the network model's doing.
function tp_host_envelope_figure(m = "lambda")
    p = plot(xlabel = "hour of day", ylabel = "voltage (p.u.)", xticks = 0:3:24,
             xlims = (0, 24), legend = :bottomleft, legend_columns = 2,
             title = "Feeder voltage envelope: LinDist3Flow vs IVACOPF")
    for ((hname, res), c) in zip(TPHOSTS, (:orangered, :navy))
        r  = res[m]
        hi = [maximum(collect(Float64, r.Vmax_t[φ])[t] for φ in 1:3) for t in eachindex(hours)]
        lo = [minimum(collect(Float64, r.Vmin_t[φ])[t] for φ in 1:3) for t in eachindex(hours)]
        plot!(p, hours, hi, lw = 2, color = c, label = "$hname max")
        plot!(p, hours, lo, lw = 2, ls = :dash, color = c, label = "$hname min")
    end
    hline!(p, [tpc.Vmin_limit, tpc.Vmax_limit], ls = :dot, lw = 1.5, color = :red,
           label = "limits")
    p
end
```

## Prerequisites

### Getting Julia

This package requires **Julia 1.10 or newer**, and is tested on 1.10 and on the current
release. If you do not have Julia yet, install it with
[juliaup](https://github.com/JuliaLang/juliaup), or via the
[Microsoft Store](https://apps.microsoft.com/detail/9NJNWW8PVKMN) (or
`winget install julia -s msstore`) on Windows, or with
`curl -fsSL https://install.julialang.org | sh` on macOS and Linux.

Alternatively, take an installer from
[julialang.org/downloads](https://julialang.org/downloads/).

### Get the code

The six example scripts, both feeders, the load and irradiance profiles and the committed
results all live in the repository, so the first step is to clone it. Every command on
this page is run from the directory this creates:

```bash
git clone https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl
cd SmartInverter-3P-DOPF-POWERUP.jl
```

### Choosing an environment

Julia installs packages into an *environment*, and a fresh environment starts out empty.
Plain `julia` uses the shared default environment; `julia --project=.` uses the one
described by the `Project.toml` in the current directory.

The three-phase example scripts are **standalone**: they carry their own `Project.toml`,
they do not depend on this repository being installed as a package, and one command
installs everything they need:

```bash
julia --project=examples/three_phase -e "using Pkg; Pkg.instantiate()"
```


### Solvers

Two of the three encodings produce a mixed-integer linear program (MILP) and need an MILP
solver; the third produces a nonlinear program (NLP) and needs an NLP solver. Both
solvers used here are in the General registry:

```julia
julia --% --project=. -e "using Pkg; Pkg.add([""Gurobi"", ""Ipopt""])"
```

| encoding | model class | solver used here |
|:--|:--|:--|
| Big-M, Lambda / SOS2 | MILP | Gurobi |
| Heaviside | NLP | Ipopt |

Gurobi is commercial and needs a licence, and is
[free for academic users](https://www.gurobi.com/academia/academic-program-and-licenses/).
Ipopt is open source and needs no licence, so the Heaviside route runs with no
commercial software at all.

#### Installing Gurobi

`Pkg.add("Gurobi")` installs the wrapper and, with it, the Gurobi binaries from
[`Gurobi_jll`](https://github.com/jump-dev/Gurobi_jll.jl), so there is no separate solver
download to do. What it does not install is a **licence**, and the size-limited trial
licence that ships with those binaries is nowhere near enough for the models here: the
case study builds a few hundred thousand variables and the scalability check 3.3 million.

To get one, register at [gurobi.com](https://www.gurobi.com) and request a licence,
which is [free for academics](https://www.gurobi.com/academia/academic-program-and-licenses/).
What you do next depends on the licence type.

A **Web License Service (WLS)** licence is a file named `gurobi.lic`, holding your
`WLSACCESSID`, `WLSSECRET` and `LICENSEID`. Save it in your home directory and nothing
further is needed:

| | home directory | the file goes at |
|:--|:--|:--|
| Windows | `C:\Users\<you>`, that is `%USERPROFILE%` | `C:\Users\<you>\gurobi.lic` |
| macOS | `/Users/<you>` | `~/gurobi.lic` |
| Linux | `/home/<you>` | `~/gurobi.lic` |

To keep it somewhere else, set the `GRB_LICENSE_FILE` environment variable to the file's
full path and Gurobi will read it from there instead.

A **named-user** licence is fetched with `grbgetkey`:

```julia
using Pkg
Pkg.add("Gurobi_jll")
import Gurobi_jll
key = "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"     # your own key
run(`$(Gurobi_jll.grbgetkey()) $key`)
```

If you already run a full Gurobi installation of your own and want Julia to use that
instead of the bundled binaries, point `GUROBI_HOME` at it, opt out, and rebuild:

```julia
ENV["GUROBI_HOME"] = "C:\\gurobi1200\\win64"    # or /Library/gurobi1200/macos_universal2
ENV["GUROBI_JL_USE_GUROBI_JLL"] = "false"
Pkg.add("Gurobi")
Pkg.build("Gurobi")
```

To confirm the licence is live, solve something trivial:

```julia
using JuMP, Gurobi
model = Model(Gurobi.Optimizer)
@variable(model, x >= 0)
@objective(model, Min, x)
optimize!(model)        # prints the licence banner, then reports OPTIMAL
```

Ipopt needs none of this. `Pkg.add("Ipopt")` is the whole installation.


## Why an "if-else" cannot go straight into a solver

The Volt-VAr law is a definition by cases, and that is exactly what a solver cannot read.
Two distinct obstacles follow.

**Conditional logic.** Which of the five expressions applies depends on ``v_i``, which is
itself a decision variable. Branching on the value of an unknown is not an algebraic
constraint, and an algebraic constraint is the only thing a solver accepts.

**Non-differentiability.** Even setting the branching aside, the slope jumps at every
breakpoint. Newton and interior-point methods build their steps from derivatives, and at
a kink the derivative does not exist.

## Try it yourself

Before moving to solver-compatible formulations, try the two minimal examples:

1. [`ifelse_numericworks.jl`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl/blob/main/examples/minimal/ifelse_numericworks.jl)  
   See that the ordinary `if-else` Volt-VAr function works when voltage is a known numerical value.

2. [`ifelse_variablefails.jl`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl/blob/main/examples/minimal/ifelse_variablefails.jl)  
   Then make voltage a JuMP decision variable and observe what happens when the same function is used in a constraint.


There are two ways out, and they define the rest of this tutorial:

- **Introduce integer variables** to encode the logic exactly. The model becomes an
  MILP. This is the Big-M and Lambda / special-ordered-set-of-type-2 (SOS2) route.
- **Write the logic in closed algebraic form** using step functions. The model stays
  integer-free but becomes non-smooth, so it needs an NLP solver. This is the Heaviside
  route.

### The one scalar the droop needs from the network

The droop is a self-contained module. Whatever DOPF you use, it exposes a voltage
magnitude at each inverter terminal; the droop module adds the relationship tying that
inverter's reactive output to that voltage:

```
        ┌────────────────────────────┐
        │       three-phase DOPF     │
        │  network model + limits    │
        └───────┬────────────▲───────┘
  exposes v_b^φ │            │ q_i^G sets
        ┌───────▼────────────┴───────┐
        │      Q-V droop module      │
        │   the IEEE 1547 curve      │
        └────────────────────────────┘
```

On a three-phase network that interface needs one sentence of care, and it is the only
place phases enter the droop at all. Take the bus set ``\Upsilon``, the phase set
``\Psi = \{a,b,c\}`` and the inverter fleet ``\mathcal{G}``. A rooftop inverter is a
single-phase device connected line-to-neutral, so inverter ``i \in \mathcal{G}`` has a
bus ``b(i) \in \Upsilon`` **and a phase** ``\varphi(i) \in \Psi``, and the voltage it
senses is that bus on that phase alone:

```math
v_i \;:=\; v_{b(i)}^{\varphi(i)},
\qquad i \in \mathcal{G}
```

Every constraint in the three sections below is written in terms of this one scalar per
inverter per time step, together with that inverter's reactive output ``q_i^{G}`` and its
reactive capability ``\bar q_i``. The time index ``t \in \{1,\dots,T\}`` is suppressed
throughout; on the case study ``T = 96``, a full day at 15-minute resolution.



## Method A — Big-M

*Following Savasci, Inaolaji and Paudyal [[5]](#ref-5), where this formulation was introduced for
a second-order-cone DOPF; also Chapter 4 of Inaolaji's dissertation [[9]](#ref-9).*

**The idea in one sentence.** Give every segment its own on/off switch, and write
constraints that are *switched off*, made trivially true, whenever their segment is
not the active one.

That switching-off is what "big-M" means. Take any constraint you want to enforce only
when a binary ``\delta`` equals 1, and add ``M(1-\delta)`` to its right-hand side. If
``\delta = 1`` the added term vanishes and the constraint bites. If ``\delta = 0`` the
right-hand side becomes so large that the constraint cannot possibly be violated; it is
still *present* in the model, but it no longer restricts anything. One constant, ``M``,
buys you an if-statement.

Everything below is written per inverter ``i \in \mathcal{G}``, and the voltage it reasons
about is the one its own phase at its own bus, ``v_i = v_{b(i)}^{\varphi(i)}``. A fleet of
twelve single-phase inverters spread over three phases therefore carries twelve
independent copies of this system at every time step.

**Step 1: exactly one segment is active.** Introduce a binary ``\delta_{i,b}`` for each of
the five segments of inverter ``i`` and require

```math
\sum_{b=1}^{5}\delta_{i,b}=1,
\qquad \forall i \in \mathcal{G} . \tag{2}
```

**Step 2: each switch owns a voltage window.** If segment ``b`` is the active one, then
the sensed voltage must lie in that segment's range ``[V^{\text{bp}}_{b},
V^{\text{bp}}_{b+1}]``. Written once with the phase visible, so there is no doubt which
voltage is meant,

```math
V^{\text{bp}}_{b} - M(1-\delta_{i,b}) \;\le\; v_{b(i)}^{\varphi(i)}
   \;\le\; V^{\text{bp}}_{b+1} + M(1-\delta_{i,b}),
```

and in big-M form that is one two-sided inequality per segment. Writing all five out, with
``v_i`` for ``v_{b(i)}^{\varphi(i)}`` from here on, gives the complete window system:

```math
\begin{aligned}
-(1-\delta_{i,1})M + V^{l}\;\; &\le v_i \le\;\; V^{\text{bp}}_{2} + (1-\delta_{i,1})M\\
-(1-\delta_{i,2})M + V^{\text{bp}}_{2} &\le v_i \le\;\; V^{\text{bp}}_{3} + (1-\delta_{i,2})M\\
-(1-\delta_{i,3})M + V^{\text{bp}}_{3} &\le v_i \le\;\; V^{\text{bp}}_{4} + (1-\delta_{i,3})M\\
-(1-\delta_{i,4})M + V^{\text{bp}}_{4} &\le v_i \le\;\; V^{\text{bp}}_{5} + (1-\delta_{i,4})M\\
-(1-\delta_{i,5})M + V^{\text{bp}}_{5} &\le v_i \le\;\; V^{u} + (1-\delta_{i,5})M
\end{aligned} \tag{3}
```

Each row is vacuous when its ``\delta_{i,b} = 0`` and binding when ``\delta_{i,b} = 1``, so
together with Step 1 the solver is forced to pick the segment that genuinely contains
``v_i``. Note the two outer rows: the first segment is bounded below by the voltage
variable's own lower bound ``V^{l}`` and the last above by ``V^{u}``, rather than by
``V^{\text{bp}}_{1}`` and ``V^{\text{bp}}_{6}``. That keeps the model feasible if ``v_i``
ever sits outside the range the curve was drawn over; the saturated laws simply continue
to apply.

**Step 3: assemble the droop law, and watch it turn nonlinear.** With the switches in
place, ``q_i^G`` is just the sum of the five segment laws, each weighted by its own
binary. Segments 1, 3 and 5 contribute constants (``\bar q_i``, ``0``, ``-\bar q_i``);
the two sloped segments contribute their affine laws, written in slope-intercept form:

```math
q_i^G \;=\; \delta_{i,1}\,\bar q_i
\;+\; \delta_{i,2}\!\left(\alpha_{i,1} v_i + \frac{\bar q_i V^{\text{bp}}_3}{V^{\text{bp}}_3 - V^{\text{bp}}_2}\right)
\;+\; \delta_{i,3}\cdot 0
\;+\; \delta_{i,4}\!\left(\alpha_{i,2} v_i + \frac{\bar q_i V^{\text{bp}}_4}{V^{\text{bp}}_5 - V^{\text{bp}}_4}\right)
\;+\; \delta_{i,5}\left(-\bar q_i\right) \tag{4}
```

with slopes ``\alpha_{i,1} = -\bar q_i/(V^{\text{bp}}_3-V^{\text{bp}}_2)`` and
``\alpha_{i,2} = -\bar q_i/(V^{\text{bp}}_5-V^{\text{bp}}_4)``. Both are inverter-specific,
because ``\bar q_i`` is: on the case study below the fleet carries four size classes, so
four different pairs of slopes appear in one model.

This is a correct statement of the curve: exactly one ``\delta_{i,b}`` equals 1, so exactly
one bracket survives and ``q_i^G`` takes that segment's value. But it is **not linear**.
Multiply the two sloped brackets out and the offending terms appear:

```math
\underbrace{\delta_{i,2}\,\alpha_{i,1} v_i}_{\text{bilinear}} \qquad\text{and}\qquad
\underbrace{\delta_{i,4}\,\alpha_{i,2} v_i}_{\text{bilinear}}
```

Each is a **product of two decision variables**: one binary, one continuous. Everything
else in the expression is a variable times a constant. So the whole difficulty of the
Big-M formulation reduces to these two products, and if they can be removed the model
becomes a plain MILP.

**Step 4: remove the two products, exactly.** The saving grace is that ``\delta_{i,b}`` is
binary rather than merely continuous, and ``v_i`` is bounded. Under those two conditions
each product can be replaced by a new continuous variable ``W_{i,b} := \delta_{i,b} v_i``
and four linear inequalities, with **no approximation whatsoever**:

```math
-M(1-\delta_{i,b}) \;\le\; v_i - W_{i,b} \;\le\; M(1-\delta_{i,b}), \qquad
V^{\text{bp}}_{b}\,\delta_{i,b} \;\le\; W_{i,b} \;\le\; V^{\text{bp}}_{b+1}\,\delta_{i,b} . \tag{5}
```

Check the two cases and the exactness is immediate. If ``\delta_{i,b} = 1``, the left pair
forces ``W_{i,b} = v_i`` and the right pair confines ``v_i`` to the segment. If
``\delta_{i,b} = 0``, the right pair forces ``W_{i,b} = 0`` (both bounds collapse to zero)
while the left pair goes slack. Either way ``W_{i,b}`` equals ``\delta_{i,b} v_i`` exactly:
this is a reformulation, not a relaxation.

Only segments 2 and 4 need this treatment, and for those the ``W_{i,b}`` bounds already pin
``v_i`` into the segment, so their Step-2 window rows are replaced rather than added to.
The complete constraint system for the Big-M droop is therefore:

```math
\begin{aligned}
-(1-\delta_{i,1})M + V^{l}\;\; &\le v_i \le\; V^{\text{bp}}_{2} + (1-\delta_{i,1})M\\[2pt]
-M(1-\delta_{i,2}) \;&\le\; v_i - W_{i,2} \;\le\; (1-\delta_{i,2})M\\
V^{\text{bp}}_{2}\,\delta_{i,2} \;&\le\; W_{i,2} \;\le\; V^{\text{bp}}_{3}\,\delta_{i,2}\\[2pt]
-(1-\delta_{i,3})M + V^{\text{bp}}_{3} &\le v_i \le\; V^{\text{bp}}_{4} + (1-\delta_{i,3})M\\[2pt]
-M(1-\delta_{i,4}) \;&\le\; v_i - W_{i,4} \;\le\; (1-\delta_{i,4})M\\
V^{\text{bp}}_{4}\,\delta_{i,4} \;&\le\; W_{i,4} \;\le\; V^{\text{bp}}_{5}\,\delta_{i,4}\\[2pt]
-(1-\delta_{i,5})M + V^{\text{bp}}_{5} &\le v_i \le\; V^{u} + (1-\delta_{i,5})M
\end{aligned} \tag{6}
```

Read alongside the Step-2 system, the change is visible: rows 2 and 4, the sloped
segments, have each become a ``W`` definition plus a ``W`` range, while the three flat
segments keep their original windows unchanged.

Now substitute ``\delta_{i,2} v_i \to W_{i,2}`` and ``\delta_{i,4} v_i \to W_{i,4}`` in the
Step 3 expression. Nothing else changes, and the droop law becomes a single **linear**
equation in which every coefficient is a constant:

```math
q_i^G = \delta_{i,1}\bar q_i
      + \alpha_{i,1} W_{i,2} + \delta_{i,2}\frac{\bar q_i V^{\text{bp}}_3}{V^{\text{bp}}_3 - V^{\text{bp}}_2}
      + \alpha_{i,2} W_{i,4} + \delta_{i,4}\frac{\bar q_i V^{\text{bp}}_4}{V^{\text{bp}}_5 - V^{\text{bp}}_4}
      - \delta_{i,5}\bar q_i,
\qquad \forall i \in \mathcal{G} \tag{7}
```

Compare it with the Step 3 version: the two bracketed sloped terms have simply been split
into a ``W`` term and a ``\delta`` term. That substitution is the entire content of the
Big-M droop model.


!!! tip "Choose M as tightly as you can justify"
    ``M`` only has to dominate the largest possible violation of a deactivated
    constraint, which here is set by the voltage bounds. A needlessly large ``M`` leaves
    the LP relaxation loose, the branch-and-bound tree deep, and the solve slow. The
    value used here is `1.1`.


## Method B — Lambda / SOS2

*Following Inaolaji, Savasci and Paudyal [[6]](#ref-6) and its three-phase extension [[7]](#ref-7), which
apply the classical lambda method to Volt-VAr and Volt-Watt droops on a LinDistFlow
host; see also Chapter 5 of [[9]](#ref-9).*

**The idea in one sentence.** Instead of asking *which segment am I on*, describe the
operating point directly as a blend of two neighbouring breakpoints.

Big-M starts from the case distinction and works to make it linear. Lambda never forms
the case distinction at all. It uses a fact about piecewise-linear curves: **every point
on the curve is a weighted average of two adjacent breakpoints**, and nothing else is.

So attach a weight ``\lambda_{i,b} \ge 0`` to each of the six breakpoints of inverter
``i``, make the weights sum to one, and build *both* coordinates from the same weights.
Written once with the phase visible,

```math
v_{b(i)}^{\varphi(i)} = \sum_{b=1}^{6}\lambda_{i,b} V^{\text{bp}}_b,
```

and then, with ``v_i`` for that same quantity,

```math
v_i = \sum_{b=1}^{6}\lambda_{i,b} V^{\text{bp}}_b, \qquad
q_i^G = \sum_{b=1}^{6}\lambda_{i,b}\, q^{\text{bp}}_{i,b}, \qquad
\sum_{b=1}^{6}\lambda_{i,b} = 1, \qquad \lambda_{i,b} \ge 0 . \tag{8}
```

The single shared ``\lambda`` is the whole trick. Because one set of weights generates
the voltage *and* the reactive power, the pair ``(v_i, q_i^G)`` cannot drift off the
curve: move the weights and both coordinates slide together along it. Both
``V^{\text{bp}}_b`` and the ordinates ``q^{\text{bp}}_{i,b} = \bar q_i\,(q/\bar q)_b`` are
constants, so these are ordinary linear constraints, and no ``M`` needs choosing anywhere.
Note that only the ordinates carry the inverter index: the six breakpoint *voltages* are
the utility's setting and are shared by the whole fleet, while the six reactive values are
scaled by each inverter's own capability.

**The catch.** As written, the weights describe the *convex hull* of the six
breakpoints, not the curve. Nothing yet stops the solver putting weight on
``\lambda_{i,1}`` and ``\lambda_{i,5}`` simultaneously, which lands the operating point
somewhere in the interior of that hull, a ``(v, q)`` pair the inverter would never
produce. Since interior points give the optimiser more reactive power at a given voltage
than the real device offers, it will happily take them.

**The fix** is the classical **SOS2** condition: at most two weights may be nonzero, and
they must be *adjacent*. That is exactly the "blend of two neighbouring breakpoints"
statement, imposed rather than hoped for. Introduce one binary ``z_{i,b}`` per segment,
five of them for six breakpoints, and write, in full:

```math
\begin{aligned}
\lambda_{i,1} &\le z_{i,1}\\
\lambda_{i,2} &\le z_{i,1} + z_{i,2}\\
\lambda_{i,3} &\le z_{i,2} + z_{i,3}\\
\lambda_{i,4} &\le z_{i,3} + z_{i,4}\\
\lambda_{i,5} &\le z_{i,4} + z_{i,5}\\
\lambda_{i,6} &\le z_{i,5}\\[2pt]
\sum_{b=1}^{5} z_{i,b} &= 1, \qquad z_{i,b} \in \{0,1\}
\end{aligned} \tag{9}
```

Read it as: ``z_{i,b} = 1`` names the active segment; a weight ``\lambda_{i,b}`` is
allowed to be nonzero only if breakpoint ``b`` is an endpoint of that segment. Since
exactly one ``z_{i,b}`` is 1, precisely two adjacent weights survive and every other weight
is forced to zero. The blend is back on the curve.

Trace one case to see it work. Suppose ``z_{i,3} = 1`` and every other ``z_{i,b} = 0``.
Rows 1, 2 and 6 then force ``\lambda_{i,1} = \lambda_{i,2} = \lambda_{i,6} = 0``; row 5
forces ``\lambda_{i,5} = 0``; and only ``\lambda_{i,3} \le 1`` and ``\lambda_{i,4} \le 1``
survive. With ``\sum_b \lambda_{i,b} = 1`` the operating point is a blend of breakpoints 3
and 4 alone, that is, a point on segment 3, the dead-band.

Collecting everything, the complete Lambda droop model is:

```math
\begin{aligned}
v_i &= \sum_{b=1}^{6}\lambda_{i,b} V^{\text{bp}}_b\\
q_i^G &= \sum_{b=1}^{6}\lambda_{i,b}\, q^{\text{bp}}_{i,b}\\
\sum_{b=1}^{6}\lambda_{i,b} &= 1, \qquad \lambda_{i,b} \ge 0\\
\lambda_{i,1} \le z_{i,1}, \quad \lambda_{i,b} &\le z_{i,b-1} + z_{i,b} \;\;(b=2,\dots,5), \quad \lambda_{i,6} \le z_{i,5}\\
\sum_{b=1}^{5} z_{i,b} &= 1, \qquad z_{i,b} \in \{0,1\}
\end{aligned}
\qquad \forall i \in \mathcal{G} \tag{10}
```


## Method C — Heaviside

*Following Inaolaji, Savasci and Paudyal [[10]](#ref-10), which introduced this encoding precisely
to remove the integer variables from the two formulations above, on a
current-voltage DOPF host of the same family used here; see also Chapter 6 of [[9]](#ref-9).*

Both previous methods spend integer variables to answer "which segment?". Integers are
what make a model combinatorial: the count grows with inverters × time steps, and
branch-and-bound has to search over them. On an unbalanced LV feeder that product is the
whole problem, because the fleet is made of single-phase devices and there can be one at
every service connection. The motivation in [[10]](#ref-10) is to get rid of the integers
altogether, which also makes the model a candidate for real-time use.

The observation is that an "if" is just an on/off switch, and the unit step *is* an
on/off switch written as a function:

```math
H(x) = \begin{cases} 1, & x \ge 0\\ 0, & x < 0\end{cases} \tag{12}
```

Shift it to flip at a breakpoint and subtract two of them, and you get a **window** that
equals 1 on one segment and 0 everywhere else:

```math
\mathcal{W}_{i,b} \;=\; H\!\left(v_i - V^{\text{bp}}_{b}\right) - H\!\left(v_i - V^{\text{bp}}_{b+1}\right),
\qquad v_i = v_{b(i)}^{\varphi(i)} \tag{13}
```

which is precisely the condition ``V^{\text{bp}}_b \le v_i \le V^{\text{bp}}_{b+1}``: the
if-else of segment ``b``, written without logic and without binaries. Multiply each
segment's law by its own window and add them up. The windows are disjoint, so at any
voltage all but one vanish and the sum collapses to the single active law.

Written out with every window expanded, and with all five segments present so the
structure is visible:

```math
\begin{aligned}
q_i^G \;=\; &\;\;\;\;\bar q_i \big[\,H(v_i - V^{\text{bp}}_1) - H(v_i - V^{\text{bp}}_2)\,\big] \;+\\
&\;\alpha_{i,1}\!\left(v_i - V^{\text{bp}}_3\right)\big[\,H(v_i - V^{\text{bp}}_2) - H(v_i - V^{\text{bp}}_3)\,\big] \;+\\
&\;\;\;\;0\,\big[\,H(v_i - V^{\text{bp}}_3) - H(v_i - V^{\text{bp}}_4)\,\big] \;+\\
&\;\alpha_{i,2}\!\left(v_i - V^{\text{bp}}_4\right)\big[\,H(v_i - V^{\text{bp}}_4) - H(v_i - V^{\text{bp}}_5)\,\big] \;-\\
&\;\;\;\;\bar q_i \big[\,H(v_i - V^{\text{bp}}_5) - H(v_i - V^{\text{bp}}_6)\,\big]
\end{aligned}
\qquad \forall i \in \mathcal{G} \tag{14}
```

with the same slopes as before,
``\alpha_{i,1} = -\bar q_i/(V^{\text{bp}}_3-V^{\text{bp}}_2)`` and
``\alpha_{i,2} = -\bar q_i/(V^{\text{bp}}_5-V^{\text{bp}}_4)``.



`op_ifelse` and `op_greater_than_or_equal_to` are JuMP's nonlinear operators
(JuMP ≥ 1.15); they build the expression correctly outside a macro.

No extra variables at all, just one algebraic expression per inverter per time step. The
price is paid in solver behaviour. ``H(\cdot)`` is discontinuous, so the derivative is
undefined at every breakpoint and the problem is non-convex. Two consequences follow: the
model needs an NLP solver rather than an MILP one, and the non-smoothness is expensive to
differentiate, which makes this the slowest of the three encodings on the case study and
the first to break down as the network grows. [Does it scale?](@ref) puts numbers on both.

## The three-phase hosts

The droop needs a host, and this tutorial provides two, deliberately, because the pair
makes the separation between *encoding* and *host* measurable rather than merely asserted:

**Table 3.** The two three-phase hosts, and the script family implementing each.

| script family | model | class | solve |
|:--|:--|:--|:--|
| `LinDist3Flow_*.jl` | **LinDist3Flow**: multiphase linearised branch flow [[12]](#ref-12) | linear approximation | one pass |
| `IVACOPF3Ph_*.jl` | **IVACOPF**: three-phase current-voltage AC-OPF [[4]](#ref-4) | near-exact AC | successive linearisation, iterated |

Both are set out in full below, briefly, because the subject of this page is the droop
rather than the network model, and both carry the three droop blocks just written, the
same feeder, the same fleet and the same objective.


### LinDist3Flow: the linear host

The multiphase form of the LinDistFlow linearisation [[3]](#ref-3) of the Baran and Wu
branch-flow model [[2]](#ref-2), from Sankur, Dobbe, Stewart, Callaway and Arnold
[[12]](#ref-12). Each line carries a 3×3 phase impedance ``Z`` rather than a scalar and
the phases couple, so a scalar ``rP + xQ`` drop no longer suffices.

That paper writes KVL and KCL in three-phase vector form and derives the exact
**Dist3Flow** equations, Eqs. (14)–(17) of [[12]](#ref-12). Two things in them are
nonlinear: the loss terms, and the *ratio of voltages between phases* at a node,
``\gamma_n^{\varphi\psi} = V_n^{\varphi}/V_n^{\psi}``, which scales and rotates the
off-diagonal impedances. **LinDist3Flow is what follows from holding both constant**, under the
paper's assumptions **A1** (``\gamma`` constant) and **A2** (loss terms constant). Fixing
``\gamma`` at its nominal value, ``1\angle{\pm}120^{\circ}``, and dropping losses
altogether gives Eqs. (20)–(23) of [[12]](#ref-12), the model used here. Per phase
``\varphi``:

```math
w_j^{\varphi} = w_i^{\varphi} - \sum_{\psi} \Big( a^R_{\varphi\psi} P_{ij}^{\psi}
                                               + a^X_{\varphi\psi} Q_{ij}^{\psi} \Big),
\qquad
\begin{aligned}
a^R_{\varphi\psi} &= 2\,\mathrm{Re}\!\left(\alpha^{\psi-\varphi} Z_{\varphi\psi}\right)\\
a^X_{\varphi\psi} &= 2\,\mathrm{Im}\!\left(\alpha^{\psi-\varphi} Z_{\varphi\psi}\right)
\end{aligned} \tag{16}
```

with ``w = \lvert V\rvert^2`` the squared voltage magnitude and ``\alpha = e^{-j2\pi/3}``
the 120° rotation. This is Eq. (21) of [[12]](#ref-12). Writing the rotation out term by
term recovers its coefficient matrices, Eqs. (22)–(23), exactly, with
``a^R = -\mathbb{M}^P`` and ``a^X = -\mathbb{M}^Q``:

```math
a^R_{ij} =
\begin{bmatrix}
 2r^{aa} & -r^{ab}+\sqrt{3}\,x^{ab} & -r^{ac}-\sqrt{3}\,x^{ac}\\
-r^{ba}-\sqrt{3}\,x^{ba} &  2r^{bb} & -r^{bc}+\sqrt{3}\,x^{bc}\\
-r^{ca}+\sqrt{3}\,x^{ca} & -r^{cb}-\sqrt{3}\,x^{cb} &  2r^{cc}
\end{bmatrix} \tag{17}
```

and ``a^X`` identically, with ``r`` and ``x`` exchanged and the sign of every
``\sqrt{3}`` term flipped. Those ``\pm\sqrt{3}`` cross-terms are the ``120^{\circ}``
rotation written out, and they are what makes this a *three-phase* model rather than three
single-phase ones running side by side. Two checks are worth carrying: for a single phase
``\alpha^0 = 1`` gives ``a^R = 2r`` and ``a^X = 2x``, recovering
``w_j = w_i - 2(rP + xQ)``; and for diagonal ``Z`` the matrices are diagonal and the
phases decouple into three independent LinDistFlows.

The implementation works in magnitude rather than squared magnitude
(``w_j - w_i \approx 2 V^{\mathrm{nom}}(v_j - v_i)`` near nominal) so that the droop
breakpoints stay in ordinary p.u. voltage. That done, the complete host is four
equations: the lossless power balance per bus *and* per phase, which is Eq. (20) of
[[12]](#ref-12), the coupled drop above, the slack, and the voltage limits its DOPF,
Eq. (24), imposes:

```math
\begin{aligned}
v_0^{\varphi} &= V^{\mathrm{nom}} & &\forall \varphi \in \Psi\\
p_j^{G,\varphi} - p_j^{L,\varphi} &= \sum_{k:(j,k)\in\mathcal{L}} P_{jk}^{\varphi}
                                   - \sum_{i:(i,j)\in\mathcal{L}} P_{ij}^{\varphi}
                                   & &\forall j \in \mathcal{B},\ \varphi \in \Psi\\
q_j^{G,\varphi} - q_j^{L,\varphi} &= \sum_{k:(j,k)\in\mathcal{L}} Q_{jk}^{\varphi}
                                   - \sum_{i:(i,j)\in\mathcal{L}} Q_{ij}^{\varphi}
                                   & &\forall j \in \mathcal{B},\ \varphi \in \Psi\\
v_j^{\varphi} &= v_i^{\varphi} - \sum_{\psi \in \Psi} \Big( \tilde a^R_{\varphi\psi} P_{ij}^{\psi}
                                       + \tilde a^X_{\varphi\psi} Q_{ij}^{\psi} \Big)
                                   & &\forall (i,j) \in \mathcal{L},\ \varphi \in \Psi\\
V^{\min} &\le v_j^{\varphi} \le V^{\max} & &\forall j \in \mathcal{B},\ \varphi \in \Psi
\end{aligned} \tag{18}
```

with ``\Psi = \{a,b,c\}`` the phase set, ``\varphi`` and ``\psi`` phases within it, and
``\tilde a = a / (2V^{\mathrm{nom}})`` the magnitude-form coefficients. Note what is
*absent*: there is no current variable and no loss term. That is exactly what buys the
linearity, and exactly what it costs.



### Three-phase IVACOPF: the near-exact host

The **current-voltage AC optimal power flow** of Soltani, Khorsand and Ma [[4]](#ref-4),
in its native three-phase unbalanced form, the setting [[4]](#ref-4) was written for. Its
appeal here is structural. Write the network in rectangular current and voltage
coordinates and the *line* equations become exactly linear, mutual coupling and all; the
only nonlinearity left is the ``v \cdot I`` power balance and the voltage magnitude, and
both live at the **buses**. In a distribution feeder the buses are the endpoints and the
lines are everything else, so this confines the nonlinearity to a small, well-behaved part
of the model instead of spreading it along every branch, as a power-voltage formulation
does.

The formulation is that of Soltani, Khorsand and Ma [[4]](#ref-4); the equations are
numbered here in this tutorial's own sequence. Bus set ``\Upsilon``, phase set
``\Psi = \{a,b,c\}`` with phases indexed ``\varphi`` and ``p``; the time index
``t \in \{1,\dots,96\}`` is suppressed throughout.

**Line current constraints.** For the line from ``n`` to ``m``, with 3×3
impedance ``Z_{nm}^{\varphi p} = R_{nm}^{\varphi p} + jX_{nm}^{\varphi p}`` and shunt
admittance ``y_{nm}^{p,k}``:

```math
V_n^{\varphi} - V_m^{\varphi} = \sum_{p\in\Psi} Z_{nm}^{\varphi p} I_{nm}^{p}
   \;-\; \tfrac{1}{2}\sum_{p\in\Psi} Z_{nm}^{\varphi p}
          \Big( \sum_{k\in\Psi} y_{nm}^{p,k} V_n^{k} \Big),
\qquad \forall \varphi \in \Psi \tag{19}
```

Three physical contributions, in two sums: the current in the same phase (the ``p = \varphi``
term), the currents in the *other* phases reaching this one through the mutual impedances,
and the shunt current. Splitting (19) into real and imaginary parts gives, for the
Kron-reduced three-wire feeders used here, where ``y = 0``,

```math
\begin{aligned}
v_n^{\mathrm{re},\varphi} - v_m^{\mathrm{re},\varphi}
  &= \sum_{p\in\Psi}\Big( R_{nm}^{\varphi p} I_{nm}^{\mathrm{re},p}
                        - X_{nm}^{\varphi p} I_{nm}^{\mathrm{im},p} \Big)\\
v_n^{\mathrm{im},\varphi} - v_m^{\mathrm{im},\varphi}
  &= \sum_{p\in\Psi}\Big( R_{nm}^{\varphi p} I_{nm}^{\mathrm{im},p}
                        + X_{nm}^{\varphi p} I_{nm}^{\mathrm{re},p} \Big)
\end{aligned}
\qquad \forall (n,m) \in \mathcal{L},\ \varphi \in \Psi \tag{20}
```

These are **exact and linear**. No rotation operator appears, nothing is transposed, and
no near-balance is assumed anywhere; compare the ``\alpha^{\psi-\varphi}`` of
LinDist3Flow, which is precisely where that host's balanced-voltage assumption enters.

**Bus current injection**, KCL per bus and phase, also exact and linear:

```math
I_n^{\mathrm{re},\varphi} = \sum_{m:(n,m)\in\mathcal{L}} I_{nm}^{\mathrm{re},\varphi}
                          - \sum_{k:(k,n)\in\mathcal{L}} I_{kn}^{\mathrm{re},\varphi},
\qquad
I_n^{\mathrm{im},\varphi} = \sum_{m:(n,m)\in\mathcal{L}} I_{nm}^{\mathrm{im},\varphi}
                          - \sum_{k:(k,n)\in\mathcal{L}} I_{kn}^{\mathrm{im},\varphi} \tag{21}
```

**Power balance**, the first of the two nonlinear relations:

```math
\begin{aligned}
p_n^{G,\varphi} - p_n^{L,\varphi}
   &= v_n^{\mathrm{re},\varphi} I_n^{\mathrm{re},\varphi}
    + v_n^{\mathrm{im},\varphi} I_n^{\mathrm{im},\varphi}\\
q_n^{G,\varphi} - q_n^{L,\varphi}
   &= v_n^{\mathrm{im},\varphi} I_n^{\mathrm{re},\varphi}
    - v_n^{\mathrm{re},\varphi} I_n^{\mathrm{im},\varphi}
\end{aligned}
\qquad \forall n \in \Upsilon,\ \varphi \in \Psi \tag{22}
```

**Linearised power balance.** Each product ``xy`` in (22) is replaced by its first-order
Taylor expansion about the previous iterate, ``xy \approx x^{\circ}y + y^{\circ}x -
x^{\circ}y^{\circ}``, where ``\circ`` marks a value **fixed from the previous pass**, a
constant, not a variable:

```math
\begin{aligned}
\mathcal{P}_n^{\varphi} &:= v_n^{\mathrm{re},\varphi\circ} I_n^{\mathrm{re},\varphi}
   + I_n^{\mathrm{re},\varphi\circ} v_n^{\mathrm{re},\varphi}
   + v_n^{\mathrm{im},\varphi\circ} I_n^{\mathrm{im},\varphi}
   + I_n^{\mathrm{im},\varphi\circ} v_n^{\mathrm{im},\varphi}
   - v_n^{\mathrm{re},\varphi\circ} I_n^{\mathrm{re},\varphi\circ}
   - v_n^{\mathrm{im},\varphi\circ} I_n^{\mathrm{im},\varphi\circ}\\
\mathcal{Q}_n^{\varphi} &:= v_n^{\mathrm{im},\varphi\circ} I_n^{\mathrm{re},\varphi}
   + I_n^{\mathrm{re},\varphi\circ} v_n^{\mathrm{im},\varphi}
   - v_n^{\mathrm{re},\varphi\circ} I_n^{\mathrm{im},\varphi}
   - I_n^{\mathrm{im},\varphi\circ} v_n^{\mathrm{re},\varphi}
   - v_n^{\mathrm{im},\varphi\circ} I_n^{\mathrm{re},\varphi\circ}
   + v_n^{\mathrm{re},\varphi\circ} I_n^{\mathrm{im},\varphi\circ}
\end{aligned} \tag{23}
```

which are then set equal to the net injection at each class of bus:

```math
\begin{aligned}
p_0^{\mathrm{grid},\varphi} &= \mathcal{P}_0^{\varphi}, &
q_0^{\mathrm{grid},\varphi} &= \mathcal{Q}_0^{\varphi} & &\text{substation}\\
-p_n^{L,\varphi} &= \mathcal{P}_n^{\varphi}, &
-q_n^{L,\varphi} &= \mathcal{Q}_n^{\varphi} & &\text{load-only bus and phase}\\
p_i^{G} - p_n^{L,\varphi} &= \mathcal{P}_n^{\varphi}, &
q_i^{G} - q_n^{L,\varphi} &= \mathcal{Q}_n^{\varphi} & &\text{inverter } i \text{ at } (n,\varphi)
\end{aligned} \tag{24}
```

The last line is where the droop enters the network: ``q_i^{G}`` is exactly the variable
the three encodings constrain.

**Voltage magnitude, and its linearisation**, the second nonlinear relation and the
single quantity the droop module reads:

```math
v_n^{\varphi} = \sqrt{\big(v_n^{\mathrm{re},\varphi}\big)^2 + \big(v_n^{\mathrm{im},\varphi}\big)^2}
\;\;\longrightarrow\;\;
v_n^{\varphi} = \frac{v_n^{\mathrm{re},\varphi\circ}}
   {\sqrt{\big(v_n^{\mathrm{re},\varphi\circ}\big)^2 + \big(v_n^{\mathrm{im},\varphi\circ}\big)^2}}\, v_n^{\mathrm{re},\varphi}
 + \frac{v_n^{\mathrm{im},\varphi\circ}}
   {\sqrt{\big(v_n^{\mathrm{re},\varphi\circ}\big)^2 + \big(v_n^{\mathrm{im},\varphi\circ}\big)^2}}\, v_n^{\mathrm{im},\varphi} \tag{25}
```

**Voltage limits.** The band every bus and phase must stay inside:

```math
V^{\min} \le v_n^{\varphi} \le V^{\max},
\qquad \forall n \in \Upsilon,\ \varphi \in \Psi \tag{26}
```

**Thermal line limits.** The conductor rating, per line and phase:

```math
\big(I_{nm}^{\mathrm{re},\varphi}\big)^2 + \big(I_{nm}^{\mathrm{im},\varphi}\big)^2
   \le \big(I_{nm}^{\max,\varphi}\big)^2,
\qquad \forall (n,m) \in \mathcal{L},\ \varphi \in \Psi \tag{27}
```

Constraint (27) is worth pausing on: IVACOPF carries the line current as a decision
variable, so a thermal limit is something you simply *write*. LinDist3Flow has no ``I`` to
write it about. It is quadratic, so the scripts offer it as a polygon inscribing the circle,
which keeps the model an MILP, and leave it off by default, because on these Electricity
North West (ENWL) feeders [[14]](#ref-14), [[15]](#ref-15) the peak flow is about a fifth
of the conductor rating; the loading is reported either way.

**Slack reference.** The three-phase substation, which is also the flat start prescribed
by Soltani, Khorsand and Ma [[4]](#ref-4):

```math
v_0^{\mathrm{re},\varphi} = \cos\theta_{\varphi},\quad
v_0^{\mathrm{im},\varphi} = \sin\theta_{\varphi},
\qquad \theta = (0°,\, -120°,\, +120°) \tag{28}
```

Seeding all three phases at ``1\angle 0°`` instead is a silent and expensive mistake: the
mutual terms then add rather than largely cancelling.

**Convergence.** After each pass, the linearisation error is measured against the *true*
nonlinear relations, not against the model's own residual. Following [[4]](#ref-4), three
metrics are used: the **maximum absolute active power balance** error (MAPB), the
**maximum absolute reactive power balance** error (MRPB), and the **maximum
voltage-magnitude** error (MVM). Writing (22) minus (23) and (25) exact minus linearised,

```math
\begin{aligned}
\text{MAPB} &= \max_{n\in\Upsilon,\,\varphi\in\Psi}
   \Big| \big(v_n^{\mathrm{re},\varphi} I_n^{\mathrm{re},\varphi}
            + v_n^{\mathrm{im},\varphi} I_n^{\mathrm{im},\varphi}\big)
         - \mathcal{P}_n^{\varphi} \Big|\\
\text{MRPB} &= \max_{n\in\Upsilon,\,\varphi\in\Psi}
   \Big| \big(v_n^{\mathrm{im},\varphi} I_n^{\mathrm{re},\varphi}
            - v_n^{\mathrm{re},\varphi} I_n^{\mathrm{im},\varphi}\big)
         - \mathcal{Q}_n^{\varphi} \Big|\\
\text{MVM} &= \max_{n\in\Upsilon,\,\varphi\in\Psi}
   \Big| \sqrt{\big(v_n^{\mathrm{re},\varphi}\big)^2 + \big(v_n^{\mathrm{im},\varphi}\big)^2}
         - v_n^{\varphi} \Big|
\end{aligned} \tag{29}
```

the loop repeats with a refreshed ``\circ`` point until
``\max(\text{MAPB}, \text{MRPB}, \text{MVM}) < \epsilon``, here ``10^{-6}``. Checking
against the true relations is what makes the converged point a genuine power-flow solution
rather than a solution of the approximation, and the audit further down confirms it
independently.



#### Where the linearisation loop starts

A flat start (``1\angle0°,\,1\angle{-120°},\,1\angle{+120°}`` with all currents zero) is
what Soltani, Khorsand and Ma [[4]](#ref-4) prescribe, and what the scripts fall back to
with `TP_WARMSTART=flat`. By
default they do something cheaper to converge from: one exact three-phase
backward/forward sweep per time step, at full PV and zero VArs, which costs a few seconds
and hands the first linearisation a physically consistent state instead of a guess. The
principle is worth stating on its own: a cheap, physically consistent starting point buys
passes off the outer loop, and the same sweep is reused afterwards to audit the answer.

## The three-phase case study

The case study puts twelve inverters on `network_5_Feeder_2` [[14]](#ref-14), a real
unbalanced low-voltage (LV) feeder with 194 buses and eighteen single-phase loads split four,
five and nine across the phases, in **four size classes**. Because ``\bar q_i = S_i^{\max}``,
the four classes follow four *different* droop curves: same breakpoint voltages, four
saturation levels. Each phase carries one inverter of each class, 84 kW of array in total,
over a full day at 15-minute resolution, ``T = 96`` time steps. Bus voltages are limited to
``[0.95, 1.05]`` p.u. on every phase, and the objective is to **minimise total PV
curtailment** over the day.


Every host × encoding pair has its own standalone script in
[`examples/three_phase/`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl/tree/main/examples/three_phase).
All six share their skeleton verbatim (data, PV placement, verification, figures); a
`diff` between any two shows only the droop block, or only the network model:

**Table 11.** The six three-phase example scripts, one per host and encoding.

| | Big-M | Lambda / SOS2 | Heaviside |
|:--|:--|:--|:--|
| **LinDist3Flow** | [`LinDist3Flow_BigM.jl`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl/blob/main/examples/three_phase/LinDist3Flow_BigM.jl) | [`LinDist3Flow_Lambda.jl`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl/blob/main/examples/three_phase/LinDist3Flow_Lambda.jl) | [`LinDist3Flow_Heaviside.jl`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl/blob/main/examples/three_phase/LinDist3Flow_Heaviside.jl) |
| **IVACOPF** | [`IVACOPF3Ph_BigM.jl`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl/blob/main/examples/three_phase/IVACOPF3Ph_BigM.jl) | [`IVACOPF3Ph_Lambda.jl`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl/blob/main/examples/three_phase/IVACOPF3Ph_Lambda.jl) | [`IVACOPF3Ph_Heaviside.jl`](https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl/blob/main/examples/three_phase/IVACOPF3Ph_Heaviside.jl) |

The Big-M and Lambda scripts need an MILP solver (Gurobi); the Heaviside ones need only
Ipopt.

Before running the Gurobi-based examples, make sure that Gurobi is installed
and that a valid license is available. If Gurobi does not automatically locate
your license, set the `GRB_LICENSE_FILE` environment variable to the location
of your `gurobi.lic` file.

For example, 
```bash
$env:GRB_LICENSE_FILE="C:\gurobi1300\gurobi.lic"
```

```bash
julia --project=examples/three_phase examples/three_phase/IVACOPF3Ph_Lambda.jl
```

**Table 12.** Environment overrides accepted by every three-phase script.

| variable | default | meaning |
|:--|:--|:--|
| `TP_CASE` | `network_5_Feeder_2` | ENWL feeder to load: `network_5_Feeder_2` [[14]](#ref-14) or `network_17_Feeder_6` [[15]](#ref-15) |
| `TP_STEPS` | `96` | time steps in the day |
| `TP_NPV` | `4` | smart inverters per phase |
| `TP_WARMSTART` | `sweep` | IVACOPF only: `flat` for the flat start of Soltani, Khorsand and Ma [[4]](#ref-4) |
| `TP_TOL` | `1e-6` | IVACOPF only: stop tolerance on ``\max(\text{MAPB}, \text{MRPB}, \text{MVM})``, eq. (29) |
| `TP_MAXITER` | `15` | IVACOPF only: pass limit |
| `TP_IMAXSEG` | `0` | IVACOPF only: sides of the polygon enforcing (27); 0 disables it |

```bash
TP_CASE=network_17_Feeder_6 TP_STEPS=24 julia --project=examples/three_phase examples/three_phase/IVACOPF3Ph_Lambda.jl
```



**Table 4.** The four inverter size classes of the case study on `network_5_Feeder_2` [[14]](#ref-14). Because ``\bar q_i = S_i^{\max}``, each class follows a different droop curve.

```@example tut
tp_class_table()   # hide
```

## The two hosts, side by side

Six runs: three encodings on each of two hosts, everything else held fixed.

**Table 5.** Three encodings on each of two three-phase hosts, everything else held fixed. Passes is 1 for the linear host, which has no outer loop.

```@example tut
tp_host_table()   # hide
```

Read Table 5 in two directions. **Down each host block**, the three encodings agree on
curtailment, on losses and on voltage range, which is the empirical statement of the claim
that these are three encodings of one curve. **Across the two blocks**,
the hosts do not agree, and that difference is the network model's alone.

The gap is about 3.6 kWh, some 8 % of the curtailed energy, with LinDist3Flow curtailing
more. Which host curtails more is not a general rule: it depends on the feeder, and the
direction can reverse on another one. What does not depend on the feeder is
the mechanism. LinDist3Flow drops losses from the balance entirely; IVACOPF measures them
at 14.6 kWh over the day, about 3 % of the available PV energy, and having them in the
model changes which dispatch clears the voltage band. A host that cannot represent losses
cannot be expected to agree with one that can, in either direction.

The droop is reproduced exactly in every one of the six:

**Table 6.** Exactness of the encoding within each host: the largest gap between dispatched reactive power and the curve at the voltage that host reports.

```@example tut
tp_exact_table()   # hide
```

That is the separation this section exists to make. **Exactness of the encoding is a
property of the encoding; accuracy is a property of the host.** Every cell above is at
round-off (the inverters sit on their curves *within whatever model they are placed in*),
and Table 6 says nothing whatever about whether that model is right.


## Scalability

Now, we explore whether the encodings survive a network worth
calling realistic. The three LinDist3Flow scripts were run unchanged on a second real feeder
from the same ENWL family, `network_17_Feeder_6` [[15]](#ref-15), with **3856 buses, 3855 lines,
223 single-phase loads**, twenty times `network_5_Feeder_2` [[14]](#ref-14), by setting an
environment variable:

```bash
TP_CASE=network_17_Feeder_6 julia --project=examples/three_phase     examples/three_phase/LinDist3Flow_Lambda.jl
```


`TP_HOSTS=ivacopf` or `TP_HOSTS=lindist3flow` regenerates just one family. The scalability
table has its own sweep, which shells out to the scripts one per process:

```bash
julia --project=examples/three_phase examples/three_phase/scalability.jl
```

**Table 10.** Scalability of the three encodings on the LinDist3Flow host, across the 194-bus `network_5_Feeder_2` [[14]](#ref-14) and the 3856-bus `network_17_Feeder_6` [[15]](#ref-15). “Max droop deviation” is ``\Delta`` of (30), evaluated at each host's own voltages.

```@example tut
tp_scale_table()   # hide
```

Two things to take from this.

**The mixed-integer encodings scale.** Big-M and Lambda both carry a 3.3-million-variable
model over the full 96-step day and solve it in about a minute, roughly twelve times the
small feeder's solve for eighteen times the network, and the droop is still reproduced to
solver tolerance. The binary count does not move at all between the two feeders, because
it depends on inverters × time steps and not on network size. That is the useful property:
**enlarging the network grows the linear part of the problem, not the combinatorial part.**

**The integer-free encoding does not.** Heaviside is the cheapest of the three by variable
count (it adds nothing to the model) and it is comfortably the most expensive to solve.
On the small feeder it costs several times Lambda. On the large one at the full horizon
Ipopt gives up with `ERROR`; shortening the day to twelve steps brings it back to a
model an eighth the size, which then solves, in minutes rather than the seconds the
mixed-integer encodings need, but it solves, and the last row of Table 10 records what
comes back. The non-smoothness that costs nothing to write costs a great deal to
differentiate, and it is what limits this encoding long before the network does.

None of this changes which encoding is *correct*: all three reproduce the curve exactly,
here as before. It changes which one you would reach for on a feeder with an inverter at
every service connection.

The sweep is run on LinDist3Flow, because it is the host that isolates the *encodings’*
scaling: one solve each, no outer loop, so what Table 10 measures is the cost of the droop
block and nothing else. IVACOPF multiplies every row by its pass count, three passes on the
case study, on top of a larger model per pass, but the binary counts, which are the thing
at issue here, are identical in both hosts.

The sweep runs one process per row, separately from the case study above, so its timings
will not match that table to the second. Read both as orders of magnitude, as the warning
further up says.


## References

**Standard and host models**

```@raw html
<a id="ref-1"></a>
```
**[1]** IEEE Std 1547-2018, *IEEE Standard for Interconnection and Interoperability of
Distributed Energy Resources with Associated Electric Power Systems Interfaces*.
[doi:10.1109/IEEESTD.2018.8332112](https://doi.org/10.1109/IEEESTD.2018.8332112)

```@raw html
<a id="ref-2"></a>
```
**[2]** M. E. Baran and F. F. Wu, "Network reconfiguration in distribution systems for loss
reduction and load balancing," *IEEE Transactions on Power Delivery*, vol. 4, no. 2,
pp. 1401–1407, 1989. [doi:10.1109/61.25627](https://doi.org/10.1109/61.25627).
The branch-flow (DistFlow) model that LinDistFlow linearises.

```@raw html
<a id="ref-3"></a>
```
**[3]** K. Turitsyn, P. Šulc, S. Backhaus, and M. Chertkov, "Local control of reactive power
by distributed photovoltaic generators," *2010 First IEEE International Conference on
Smart Grid Communications (SmartGridComm)*, pp. 79–84, 2010.
[doi:10.1109/SMARTGRID.2010.5622021](https://doi.org/10.1109/SMARTGRID.2010.5622021).
**LinDistFlow**, the single-phase linearisation that LinDist3Flow [[12]](#ref-12)
generalises to unbalanced multiphase networks.

```@raw html
<a id="ref-4"></a>
```
**[4]** Z. Soltani, M. Khorsand, and S. Ma, "Current–Voltage Unbalanced Distribution AC
Optimal Power Flow for Advanced Distribution Management System Applications,"
*IEEE Open Journal of Industry Applications*, vol. 5, 2024.
[doi:10.1109/OJIA.2024.3367547](https://doi.org/10.1109/OJIA.2024.3367547).
**IVACOPF**, the origin of the current-voltage host; the successive-linearisation
scheme built on it here is developed further in [[11]](#ref-11).

**Embedding the Volt-VAr droop curve in a DOPF**

```@raw html
<a id="ref-5"></a>
```
**[5]** A. Savasci, A. Inaolaji, and S. Paudyal, "Distribution Grid Optimal Power Flow
Integrating Volt-Var Droop of Smart Inverters," *2021 IEEE Green Technologies
Conference (GreenTech)*, pp. 54–59, 2021.
[doi:10.1109/GreenTech48523.2021.00020](https://doi.org/10.1109/GreenTech48523.2021.00020).
**Big-M**, on a second-order-cone DOPF.

```@raw html
<a id="ref-6"></a>
```
**[6]** A. Inaolaji, A. Savasci, and S. Paudyal, "Distribution Grid Optimal Power Flow with
Volt-VAr and Volt-Watt Settings of Smart Inverters," *2021 IEEE Industry Applications
Society Annual Meeting (IAS)*, 2021.
[doi:10.1109/IAS48185.2021.9715792](https://doi.org/10.1109/IAS48185.2021.9715792).
**Lambda / SOS2**, on a LinDistFlow host; also the source of the breakpoints and the
16-segment capability polygon used here.

```@raw html
<a id="ref-7"></a>
```
**[7]** A. Inaolaji, A. Savasci, and S. Paudyal, "Distribution Grid Optimal Power Flow in
Unbalanced Multiphase Networks with Volt-VAr and Volt-Watt Droop Settings of Smart
Inverters," *IEEE Transactions on Industry Applications*, vol. 58, no. 5, 2022.
[doi:10.1109/TIA.2022.3181110](https://doi.org/10.1109/TIA.2022.3181110).
Lambda, extended to three-phase unbalanced networks.

```@raw html
<a id="ref-8"></a>
```
**[8]** A. Savasci, A. Inaolaji, and S. Paudyal, "Distribution Grid Optimal Power Flow with
Adaptive Volt-VAr Droop of Smart Inverters," *2021 IEEE Industry Applications Society
Annual Meeting (IAS)*, 2021. [doi:10.1109/IAS48185.2021.9677119](https://doi.org/10.1109/IAS48185.2021.9677119).
Big-M with an adaptive ``Q(\Delta V)`` droop responding to temporal voltage deviation.

```@raw html
<a id="ref-9"></a>
```
**[9]** A. Inaolaji, *Accurate and Efficient Optimal Power Flow Methods with Control of Smart
Inverters*, PhD dissertation, Florida International University, 2023. A book-length
treatment covering all three encodings and the host models they sit in.

**Optimising the Volt-VAr droop curve itself**

```@raw html
<a id="ref-10"></a>
```
**[10]** A. Inaolaji, A. Savasci, and S. Paudyal, "Optimal Droop Settings of Smart Inverters,"
*2021 IEEE 48th Photovoltaic Specialists Conference (PVSC)*, pp. 2584–2589, 2021.
[doi:10.1109/PVSC43889.2021.9518650](https://doi.org/10.1109/PVSC43889.2021.9518650).
The source of the **Heaviside** encoding used here: integer-free, on a current–voltage
DOPF solved with Ipopt/JuMP. The breakpoint voltages are themselves decision variables of the DOPF rather than fixed settings, so the curve is
optimised, not merely respected.

```@raw html
<a id="ref-11"></a>
```
**[11]** R. Emami Mirak and A. Inaolaji, "Adaptive and fair optimization of smart inverter
    droop curves in distribution grids," *Electric Power Systems Research*, vol. 262,
    2027, Art. no. 113613.
    [doi:10.1016/j.epsr.2026.113613](https://doi.org/10.1016/j.epsr.2026.113613).
    Lambda / SOS2 with the breakpoints promoted to decision variables.

**Three-phase network model**

```@raw html
<a id="ref-12"></a>
```
**[12]** M. D. Sankur, R. Dobbe, E. Stewart, D. S. Callaway, and D. B. Arnold, "A
linearized power flow model for optimization in unbalanced distribution systems,"
*arXiv:1606.04492*, 2016.
[arXiv:1606.04492](https://arxiv.org/abs/1606.04492).
**LinDist3Flow**, the multiphase linearisation used for the three-phase case.

```@raw html
<a id="ref-13"></a>
```
**[13]** D. Shirmohammadi, H. W. Hong, A. Semlyen, and G. X. Luo, "A compensation-based
power flow method for weakly meshed distribution and transmission networks," *IEEE
Transactions on Power Systems*, vol. 3, no. 2, pp. 753–762, 1988.
[doi:10.1109/59.192932](https://doi.org/10.1109/59.192932).
The **backward/forward sweep** used here as the exact AC reference.

**Test feeders**

Both three-phase feeders are real Electricity North West low-voltage networks from the
*Low Voltage Network Solutions* project, Kron-reduced to three wires. They reach this
tutorial through two independent open repositories, and carry the same lineage and the
same CC BY 4.0 licence.

```@raw html
<a id="ref-14"></a>
```
**[14]** F. Geth, *BMOPFDraftData*, draft benchmark datasets for the IEEE PES Task Force
on Benchmarking Multiconductor OPF.
[github.com/frederikgeth/BMOPFDraftData](https://github.com/frederikgeth/BMOPFDraftData).
Source of `network_5_Feeder_2`
([`output/ENWLvariants/Three-wire-Kron-reduced/`](https://github.com/frederikgeth/BMOPFDraftData/tree/main/output/ENWLvariants/Three-wire-Kron-reduced)),
derived from the CSIRO four-wire LV dataset,
[doi:10.25919/jaae-vc35](https://doi.org/10.25919/jaae-vc35)

```@raw html
<a id="ref-15"></a>
```
**[15]** R. Heidari, *PMDlab.jl*, test networks and functionality built on
PowerModelsDistribution.jl.
[github.com/hei06j/PMDlab.jl](https://github.com/hei06j/PMDlab.jl).
Source of `network_17_Feeder_6`
([`data/three-wire/network_17/Feeder_6`](https://github.com/hei06j/PMDlab.jl/tree/main/data/three-wire/network_17/Feeder_6)),
used here for the scalability check

```@raw html
<a id="ref-16"></a>
```
**[16]** F. Geth, R. Heidari, and A. Koirala, "Computational analysis of impedance
transformations for four-wire power networks with sparse neutral grounding," *Proceedings
of the Thirteenth ACM International Conference on Future Energy Systems (e-Energy '22)*,
pp. 105–113, 2022.
[doi:10.1145/3538637.3538844](https://doi.org/10.1145/3538637.3538844).
The impedance transformation behind the three-wire Kron reduction of both feeders.

```@raw html
<a id="ref-17"></a>
```
**[17]** A. J. Urquhart and M. Thomson, "Cable impedance data," figshare, 2019.
[hdl:2134/15544](https://hdl.handle.net/2134/15544).
The length-normalised conductor impedances the feeders were rebuilt with.
