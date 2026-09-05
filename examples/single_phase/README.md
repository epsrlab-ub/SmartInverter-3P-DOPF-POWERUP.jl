# Single-phase examples — every host with every droop encoding

`method` picks how the Volt-VAr curve is encoded; `host` picks the network model it sits
inside. The two are independent, so there are six combinations, and each has a standalone
script here.

```bash
julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path="../..")); Pkg.instantiate()'
julia --project=. ivacopf_lambda.jl
```

The first command installs everything this folder's `Project.toml` lists: the solvers,
`Printf`, and the package itself. Every package a script names on a `using` line has to
be installed that way before it can be used.

| | Big-M | Lambda / SOS2 | Heaviside |
|:--|:--|:--|:--|
| **IVACOPF** | [`ivacopf_bigm.jl`](ivacopf_bigm.jl) | [`ivacopf_lambda.jl`](ivacopf_lambda.jl) | [`ivacopf_heaviside.jl`](ivacopf_heaviside.jl) |
| **LinDistFlow** | [`lindistflow_bigm.jl`](lindistflow_bigm.jl) | [`lindistflow_lambda.jl`](lindistflow_lambda.jl) | [`lindistflow_heaviside.jl`](lindistflow_heaviside.jl) |

Big-M and Lambda need an MILP solver (Gurobi, [free for academic
use](https://www.gurobi.com/academia/academic-program-and-licenses/)); Heaviside needs an
NLP solver (Ipopt, open source). The IVACOPF scripts pass `warm_start = :lindistflow`,
which seeds the linearisation from the linear host and roughly halves the passes.

Each script prints the model size, the solve, the curtailed energy, the voltage range, and
the check that matters: that every optimised operating point lies on the droop curve.

The model itself lives in the package; these are thin wrappers around `solve_dopf`. See the
[tutorial](https://epsrlab-ub.github.io/SmartInverter-3P-DOPF.jl/dev/tutorial_voltvar/) for the
derivations, and [`../three_phase`](../three_phase) for the LinDist3Flow version.
