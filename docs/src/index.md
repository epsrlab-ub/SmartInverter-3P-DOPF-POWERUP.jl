# SmartInverter-3P-DOPF.jl

*Modeling Smart Inverters in Three-Phase Distribution Optimal Power Flow*

A smart inverter does not take a reactive-power set-point. It follows a Volt-VAr curve
based on its own terminal voltage. A distribution optimal power flow (DOPF) that ignores
that curve returns a dispatch the inverter will never deliver.

Real low-voltage feeders are unbalanced, and rooftop inverters are single-phase devices,
so the phase an inverter sits on decides the voltage it reads. This repository puts the
curve inside a **three-phase** optimisation, three different ways, on either of two
unbalanced three-phase DOPF **host** models, and shows that the encodings agree while the
hosts do not.

Every host and encoding pair is a standalone script:

```bash
git clone https://github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl
cd SmartInverter-3P-DOPF.jl
julia --project=examples/three_phase -e 'using Pkg; Pkg.instantiate()'
julia --project=examples/three_phase examples/three_phase/IVACOPF3Ph_Lambda.jl
```

It loads the feeder, places the fleet, builds and solves the model, and prints the model
size, the curtailed energy, the per-phase voltage range and the droop-deviation check.

## The three encodings

The Volt-VAr law is a five-segment piecewise-linear function, a definition by cases,
which is exactly what a solver cannot read. We discuss three ways to rewrite it
as constraints a solver accepts:

| encoding | class | idea | needs |
|:--|:--|:--|:--|
| Big-M | MILP | one binary per segment activates that segment's voltage window and affine law | MILP solver |
| Lambda / SOS2 | MILP | the operating point is a convex combination of breakpoints, with SOS2 forcing adjacency | MILP solver |
| Heaviside | NLP | segment masks built from unit steps, summed into one closed-form expression | NLP solver |

All three are exact, reproducing the curve rather than approximating it, and on the
bundled case study they return the same dispatch to within solver tolerance. They differ
in the solver technology they demand and in how they scale with inverters × time steps.

The three-phase network enters each of them in exactly one place. An inverter connected
line-to-neutral senses the voltage of *its own phase at its own bus*, so every encoding
constrains one scalar reactive output against one scalar voltage:

```julia
vpv(i, t) = v[PV[i].bus, PV[i].phase, t]     # the voltage inverter i actually senses
```

The [Tutorial](@ref "Modeling Smart Inverters in Three-Phase Distribution Optimal Power Flow")
states all three in full three-phase form, verifies that every optimised operating point
lands on the curve, and compares them side by side.

## The two host models

The droop needs a network model to sit inside. The droop constraints are identical in
both, and each has its own script family:

| script family | model | accuracy | solve |
|:--|:--|:--|:--|
| `IVACOPF3Ph_*.jl` | **IVACOPF**: three-phase current-voltage AC-OPF | near-exact AC: exact line equations with full mutual coupling, losses modelled | iterative: re-linearised until MAPB/MRPB/MVM clear a tolerance |
| `LinDist3Flow_*.jl` | **LinDist3Flow**: multiphase linearised branch flow | approximate: losses dropped, near-balanced voltages assumed in the drop coefficients | run once, much faster and far cheaper |

!!! warning "Use IVACOPF for quantitative work"
    Audited against an exact three-phase AC power flow on the bundled case, the IVACOPF
    dispatch reproduces the true solution to ~10⁻¹¹ p.u. and sits on the droop curve to
    the same order. The LinDist3Flow dispatch is off the droop by 6 % of inverter rating,
    so the inverters would not produce the VArs it dispatched. LinDist3Flow is for a fast
    first look, for screening, and for the scaling study, not for reporting.

    On this feeder the voltage band is not binding, so that failure shows up **only** in
    the droop residual and nowhere in the constraint report. It does not have to announce
    itself, which is why the audit is worth running every time.

## The case study

Twelve inverters in four size classes on `network_5_Feeder_2`, a real Electricity North
West low-voltage feeder with 194 buses and eighteen single-phase loads split four, five
and nine across the phases, over 24 h at 15-minute resolution, minimising PV curtailment
with voltages held inside `[0.95, 1.05]` p.u. on every phase.

| host | curtailed | losses | droop residual at the **true** AC voltage |
|:--|--:|--:|--:|
| IVACOPF | 42.69 kWh | 14.61 kWh | 2.8e-11 |
| LinDist3Flow | 46.32 kWh | not modelled | 8.0e-03, **not deliverable** |

The last column is the test that separates the hosts: take each dispatch, solve the exact
three-phase AC power flow for it, and ask whether the inverters would really have produced
those VArs at the voltages they would really have seen.

## What's in the box

- **Two three-phase host models**, IVACOPF and LinDist3Flow, sharing one droop block, one
  fleet and one objective.
- Two real unbalanced ENWL low-voltage feeders, Kron-reduced to three wires:
  `network_5_Feeder_2` (194 buses) for the case study and `network_17_Feeder_6`
  (3856 buses, 223 single-phase loads) for the scalability check.
- Twelve single-phase inverters in four size classes, an inverter capability polygon, and
  a curtailment-minimising objective.
- An exact **three-phase backward/forward sweep** power flow, used both to warm-start
  IVACOPF and to audit every solved dispatch.
- Feeder, horizon and fleet driven from the environment, so the same model runs on a
  different network without editing anything.

## Installation

The three-phase scripts are standalone and every package they use is in the General
registry:

```julia
using Pkg
Pkg.add(["Gurobi", "Ipopt"])                  # solvers
Pkg.add(["JuMP", "JSON3", "Plots", "Printf"]) # modelling, case files, figures, tables
Pkg.add("LinearAlgebra")                      # the 3×3 phase impedances
```

Big-M and Lambda need an MILP solver; Heaviside needs an NLP solver such as
[Ipopt](https://github.com/jump-dev/Ipopt.jl). Every package a `using` line names has to
be added, because a fresh project environment starts empty; the tutorial's
[Prerequisites](@ref) section lists the full set and covers the Gurobi licence.

!!! note "Use Gurobi"
    The results throughout this documentation were produced with **Gurobi**, which is
    [free for academic users](https://www.gurobi.com/academia/academic-program-and-licenses/).

    We also tried the open-source MILP solvers HiGHS and GLPK on this model; neither
    worked out. One returned an infeasible status inside the successive-linearisation
    loop, the other was too slow to finish.

    The Heaviside encoding needs no MILP solver at all, only Ipopt, which is open
    source, and reaches the same answer.

## The Julia package

Alongside the three-phase scripts, `src/` carries a Julia package, still named
`SmartInverterDOPF`, that implements the same three droop encodings on two *single-phase*
hosts over a 33-bus feeder. It is not in the General registry and installs from its Git
URL:

```julia
Pkg.add(url = "https://github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl")
```

Its exported interface is documented under [API reference](@ref). This page and the
tutorial are about the three-phase model.

## Citing

If this material is useful in your work, please cite this repository:

> *SmartInverter-3P-DOPF.jl: Modeling Smart Inverters in Three-Phase Distribution Optimal Power Flow.*
> <https://github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl>

and, alongside it, the papers it builds on, listed in the tutorial's
[References](@ref) section.
