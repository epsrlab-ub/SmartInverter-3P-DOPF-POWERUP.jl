# SmartInverter-3P-DOPF.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://epsrlab-ub.github.io/SmartInverter-3P-DOPF.jl/dev/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Modeling Smart Inverters in Three-Phase Distribution Optimal Power Flow

A smart inverter does not accept a reactive-power set-point. It follows a Volt-VAr
curve, deciding from its own terminal voltage how much reactive power to inject or
absorb. A distribution optimal power flow (DOPF) that ignores that curve returns a
dispatch the inverter is never going to deliver.

Real low-voltage feeders are unbalanced, and rooftop inverters are single-phase devices:
the phase an inverter sits on decides the voltage it reads, so two identically rated
inverters can sit on different segments of the same curve at the same instant, one idle in
the dead-band and the other absorbing hard. This repository puts the curve inside a
**three-phase** optimisation, three different ways, and runs each of them against two
unbalanced three-phase DOPF host models, one linear and one near-exact. The encodings
agree; the hosts do not, and an exact three-phase AC power flow decides between them.

**[Read the tutorial →](https://epsrlab-ub.github.io/SmartInverter-3P-DOPF.jl/dev/tutorial_voltvar/)**

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

### Choosing an environment

Julia installs packages into an *environment*, and a fresh environment starts out empty.
Plain `julia` uses the shared default environment; `julia --project=.` uses the one
described by the `Project.toml` in the current directory.

The three-phase example scripts are **standalone**: they carry their own `Project.toml`,
they do not depend on this repository being installed as a package, and one command
installs everything they need. From the repository root:

```bash
julia --project=examples/three_phase -e 'using Pkg; Pkg.instantiate()'
```

### Packages

Every package the three-phase scripts use is in the General registry, so it can be added
by name:

```julia
using Pkg
Pkg.add(["JuMP", "JSON3", "Plots"])              # modelling, data files, figures
Pkg.add(["Printf", "LinearAlgebra"])             # standard library
```

| package | what it is for |
|:--|:--|
| `JuMP` | the modelling layer every formulation is written in |
| `JSON3` | reading the feeder, load and irradiance files, and the committed results |
| `Plots` | every figure |
| `Printf` | formatting the printed output and the tables |
| `LinearAlgebra` | the 3×3 phase impedances that make the network model three-phase |

The last two ship with Julia, but a project environment still has to add them before
`using` will find them.

If you also want the Julia package in `src/` (still named `SmartInverterDOPF`) rather than
only the standalone scripts, it is not in the General registry and installs from its Git
URL:

```julia
Pkg.add(url = "https://github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl")
```

### Solvers

Two of the three encodings produce a mixed-integer linear program (MILP) and need an MILP
solver; the third produces a nonlinear program (NLP) and needs an NLP solver. Both
solvers used here are in the General registry:

```julia
Pkg.add(["Gurobi", "Ipopt"])
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
which is [free for academics](https://www.gurobi.com/academia/academic-program-and-licenses/),
then follow Gurobi's
[retrieval and setup instructions](https://support.gurobi.com/hc/en-us/articles/12872879801105-How-do-I-retrieve-and-set-up-a-Gurobi-license).
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

> **On open-source MILP solvers.** We also tried HiGHS and GLPK on this model and neither
> worked out: one returned an infeasible status inside the successive-linearisation loop,
> the other was too slow to finish. The committed results use Gurobi. If no MILP licence
> is available, the Heaviside encoding needs only Ipopt and reaches the same answer.

The tutorial's
[Prerequisites](https://epsrlab-ub.github.io/SmartInverter-3P-DOPF.jl/dev/tutorial_voltvar/#Prerequisites)
section covers the same ground.

## Quick start

Every host × encoding pair is a standalone script. Run one:

```bash
julia --project=examples/three_phase examples/three_phase/IVACOPF3Ph_Lambda.jl
```

It loads the feeder, places the fleet, builds and solves the model, prints the model size,
the curtailed energy, the per-phase voltage range and the droop-deviation check, and
writes its figures alongside.

## The three encodings

The Volt-VAr law is a five-segment piecewise-linear function, a definition by cases,
which is exactly what a solver cannot read. Three standard rewrites make it tractable,
and each is stated in full three-phase form in the tutorial:

| encoding | class | idea | needs |
|:--|:--|:--|:--|
| Big-M | MILP | one binary per segment activates that segment's voltage window and affine law | MILP solver |
| Lambda / SOS2 | MILP | operating point as a convex combination of breakpoints, SOS2 forcing adjacency | MILP solver |
| Heaviside | NLP | segment masks from unit steps, summed into one closed-form expression | NLP solver |

All three are exact and, on the bundled case study, return the same dispatch to within
solver tolerance. They differ in the solver technology they demand and in how they scale
with inverters × time steps.

The three-phase network enters each of them in exactly one place. An inverter connected
line-to-neutral senses the voltage of *its own phase at its own bus*, so every encoding
constrains one scalar reactive output against one scalar voltage:

```julia
vpv(i, t) = v[PV[i].bus, PV[i].phase, t]     # the voltage inverter i actually senses
```

## The three-phase host models

The encoding picks how the droop curve is written; the host picks the network model it
sits inside. They are independent: the droop constraints are identical in both hosts.

| script family | model | accuracy | solve |
|:--|:--|:--|:--|
| `LinDist3Flow_*.jl` | multiphase linearised branch flow ([arXiv:1606.04492](https://arxiv.org/abs/1606.04492)) | approximate: losses dropped, near-balanced voltages assumed in the drop coefficients | run **once**, much faster and far lower computational effort |
| `IVACOPF3Ph_*.jl` | three-phase current-voltage AC-OPF ([doi:10.1109/OJIA.2024.3367547](https://doi.org/10.1109/OJIA.2024.3367547), extended in [doi:10.1016/j.epsr.2026.113613](https://doi.org/10.1016/j.epsr.2026.113613)) | near-exact AC: exact line equations with full mutual coupling, losses modelled | **iterative**: re-linearised until MAPB/MRPB/MVM clear a tolerance |

Six scripts, one per host and encoding:

| | Big-M | Lambda / SOS2 | Heaviside |
|:--|:--|:--|:--|
| **LinDist3Flow** | `LinDist3Flow_BigM.jl` | `LinDist3Flow_Lambda.jl` | `LinDist3Flow_Heaviside.jl` |
| **IVACOPF** | `IVACOPF3Ph_BigM.jl` | `IVACOPF3Ph_Lambda.jl` | `IVACOPF3Ph_Heaviside.jl` |

> **Use IVACOPF for quantitative work.** Audited against an exact three-phase AC power
> flow on the bundled case, the IVACOPF dispatch reproduces the true solution to ~2e-11
> p.u. and sits on the droop curve to the same order. The LinDist3Flow dispatch is off the
> droop by 6% of inverter rating, so it is not deliverable. LinDist3Flow is for a fast
> first look, for screening, and for the scaling study, where its single pass isolates the
> cost of the droop block itself.

## Case studies

Two real Electricity North West low-voltage feeders, Kron-reduced to three wires, both
over 24 h at 15-minute resolution (96 steps), both minimising PV curtailment, both with
voltages held inside `[0.95, 1.05]` p.u. on every phase.

| feeder | size | inverters | used for |
|:--|:--|:--|:--|
| `network_5_Feeder_2` | 194 buses, 18 single-phase loads split 4/5/9 | 12 in four size classes, 84 kW | the case study |
| `network_17_Feeder_6` | 3856 buses, 223 single-phase loads | 12 | the scalability check |

Sources and licence are in
[`examples/three_phase/README.md`](examples/three_phase/README.md).

Headline results, all three encodings agreeing to solver tolerance within each host:

| host | curtailed | losses | droop residual at the **true** AC voltage |
|:--|--:|--:|--:|
| IVACOPF | 42.69 kWh | 14.61 kWh | 2.8e-11 |
| LinDist3Flow | 46.32 kWh | not modelled | 8.0e-03, **not deliverable** |

The last column is the test that separates the hosts: take each dispatch, solve the exact
three-phase AC power flow for it, and ask whether the inverters would really have produced
those VArs at the voltages they would really have seen.

Feeder, horizon and fleet come from the environment (`TP_CASE`, `TP_STEPS`, `TP_NPV`, and
for IVACOPF `TP_WARMSTART`, `TP_TOL`, `TP_MAXITER`, `TP_IMAXSEG`), so the same model runs
on a different network without editing anything.

### Does it scale?

`scalability.jl` reruns the LinDist3Flow scripts on the 3856-bus feeder. The mixed-integer
encodings carry a 3.3-million-variable model over the full day in about a minute, and the
binary count does not move between feeders: it depends on inverters × time steps, not on
network size. The integer-free encoding is the one that breaks: Ipopt gives up at the full
horizon and needs a shortened day to finish. All three stay exact wherever they finish.

## Repository layout

```
examples/three_phase/   6 standalone scripts: {LinDist3Flow,IVACOPF3Ph} × {BigM,Lambda,Heaviside}
                        plus generate_results.jl, scalability.jl, plot_network.jl, and the feeders
examples/minimal/       why an if-else cannot go straight into a solver, in 40 lines
docs/                   Documenter site; builds without any solver
src/                    the SmartInverterDOPF package: case data, droop encodings, DOPF hosts
data/                   33-bus feeder, load profiles, solar profile (JSON)
scripts/                generate_results.jl for the package's own results
examples/single_phase/  6 single-phase scripts kept for reference; not covered by the tutorial
test/                   test suite
```

Every standalone script is self-contained and shares its skeleton verbatim with its
siblings, so a `diff` between any two shows only the droop block, or only the network
model.

## Regenerating the documentation results

The documentation is built from precomputed results committed under
`docs/src/assets/results/threephase/`, so building the docs needs no optimisation solver.
To recompute them, which runs all six scripts across both hosts:

```bash
julia --project=examples/three_phase examples/three_phase/generate_results.jl
```

`TP_HOSTS=ivacopf` or `TP_HOSTS=lindist3flow` regenerates one family only. The scalability
table has its own sweep:

```bash
julia --project=examples/three_phase examples/three_phase/scalability.jl
```

## Building the documentation

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## License

MIT. See [LICENSE](LICENSE).
