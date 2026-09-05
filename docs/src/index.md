# SmartInverterDOPF.jl

*Modeling Smart Inverters in Distribution Optimal Power Flow*

A smart inverter does not take a reactive-power set-point. It follows a Volt-VAr curve
based on its own terminal voltage. A distribution optimal power flow (DOPF) that ignores
that curve returns a dispatch the inverter will never deliver. This package puts the curve inside the optimisation, three
different ways, on either of two DOPF **host** models, and shows that the
encodings agree.

```julia
using SmartInverterDOPF, Gurobi

case = load_case()
res  = solve_dopf(case, Gurobi.Optimizer; method = :lambda)   # host = :ivacopf (default)

kWh(case, sum(res.PVC))     # PV energy curtailed over the day
extrema(res.V)              # voltage range across the feeder
```

The `host` keyword picks the network model; the droop encoding is unchanged by it:

```julia
solve_dopf(case, Gurobi.Optimizer; method = :lambda, host = :lindistflow)
```

## The three encodings

The Volt-VAr law is a five-segment piecewise-linear function, a definition by cases,
which is exactly what a solver cannot read. We discuss three ways to rewrite it
as constraints a solver accepts:

| `method` | class | idea | needs |
|:--|:--|:--|:--|
| `:bigm` | MILP | one binary per segment activates that segment's voltage window and affine law | MILP solver |
| `:lambda` | MILP | the operating point is a convex combination of breakpoints, with SOS2 forcing adjacency | MILP solver |
| `:heaviside` | NLP | segment masks built from unit steps, summed into one closed-form expression | NLP solver |

All three are exact, reproducing the curve rather than approximating it, and on the
bundled case study they return the same dispatch to within solver tolerance. They differ
in the solver technology they demand and in how they scale with inverters × time steps.

The [Tutorial](@ref "Modeling Smart Inverters in Distribution Optimal Power Flow") derives all
three, verifies that every optimised operating point lands on the curve, and compares
them side by side.

## The two host models

The droop needs a network model to sit inside. Either of these can be selected with
`host`, and the droop constraints are identical in both:

| `host` | model | accuracy | solve |
|:--|:--|:--|:--|
| `:ivacopf` *(default)* | **IVACOPF**: current-voltage AC-OPF | very accurate, near-exact AC | iterative: re-linearised until the residual of the exact power-flow identity clears a tolerance |
| `:lindistflow` | **LinDistFlow**: linearised branch flow | approximate: losses dropped, small voltage deviations assumed | run once, much faster and far cheaper |

!!! warning "Use `:ivacopf` for quantitative work"
    Audited against an exact AC power flow on the bundled case, the IVACOPF dispatch
    reproduces the true solution to ~10⁻⁹ p.u. and sits on the droop curve to ~10⁻⁷. The
    LinDistFlow dispatch is off the droop by 6 % of inverter rating and puts 17 of the 96
    time steps below the lower voltage limit, so it is not deliverable. `:lindistflow` is
    for a fast first look and for warm-starting, not for reporting.

The cheap host also makes the accurate one cheaper. `warm_start = :lindistflow` seeds the
IVACOPF linearisation from the LinDistFlow solution instead of a flat profile, which on
the bundled case halves the passes and is ~1.7× faster overall, for the same answer:

```julia
solve_dopf(case, Gurobi.Optimizer; method = :lambda, warm_start = :lindistflow)
```

## What's in the box

- **Two host models**, IVACOPF and LinDistFlow, interchangeable behind one keyword.
- The **33-bus** feeder over 24 h at 15-minute resolution, with per-class load
  shapes and a clear-sky PV profile.
- Three smart inverters, an inverter capability polygon, and a curtailment-minimising
  objective.
- A **backward/forward sweep** power flow for the no-inverter reference case.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl")
Pkg.add(["Gurobi", "Ipopt"])                    # solvers
Pkg.add(["JuMP", "JSON3", "Plots", "Printf"])   # modelling, case files, figures, tables
```

You will also need a solver. `:bigm` and `:lambda` need an MILP solver; `:heaviside`
needs an NLP solver such as [Ipopt](https://github.com/jump-dev/Ipopt.jl). Every other
package a `using` line names has to be added as well, because a fresh project
environment starts empty; the tutorial's [Prerequisites](@ref) section lists the full
set.

!!! note "Use Gurobi"
    The results throughout this documentation were produced with **Gurobi**, which is
    [free for academic users](https://www.gurobi.com/academia/academic-program-and-licenses/).

    We also tried the open-source MILP solvers HiGHS and GLPK on this model; neither
    worked out. One returned an infeasible status inside the successive-linearisation
    loop, the other was too slow to finish.

    The `:heaviside` encoding needs no MILP solver at all, only Ipopt, which is open
    source, and reaches the same answer.

## Citing

If this material is useful in your work, please cite this repository:

> *SmartInverterDOPF.jl: Modeling Smart Inverters in Distribution Optimal Power Flow.*
> <https://github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl>

and, alongside it, the papers it builds on, listed in the tutorial's
[References](@ref) section.
