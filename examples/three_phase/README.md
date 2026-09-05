# Three-phase Volt-VAr droop in a distribution OPF

The three-phase counterpart of the single-phase 33-bus case: a real unbalanced LV feeder,
**all three droop encodings**, and **two network models**, a linear one and a near-exact
one, so that what the encoding is responsible for can be told apart from what the network
model is.

| host | scripts | model | solve |
|:--|:--|:--|:--|
| **LinDist3Flow** | `LinDist3Flow_*.jl` | multiphase linearised branch flow (Sankur et al.) | one pass |
| **IVACOPF** | `IVACOPF3Ph_*.jl` | three-phase current-voltage AC-OPF (Soltani, Khorsand & Ma) | successive linearisation, iterated |

Run from this directory. The first line installs everything the scripts use, JuMP, JSON3,
Plots and the solvers, from this folder's `Project.toml`; a package cannot be `using`-ed
before it has been added:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"

julia --project=. LinDist3Flow_Lambda.jl        # MILP, Gurobi
julia --project=. LinDist3Flow_BigM.jl          # MILP, Gurobi
julia --project=. LinDist3Flow_Heaviside.jl     # NLP,  Ipopt

julia --project=. IVACOPF3Ph_Lambda.jl          # MILP per pass, Gurobi
julia --project=. IVACOPF3Ph_BigM.jl            # MILP per pass, Gurobi
julia --project=. IVACOPF3Ph_Heaviside.jl       # NLP  per pass, Ipopt

julia --project=. plot_network.jl               # feeder schematic
julia --project=. generate_results.jl           # refresh the documentation's data
julia --project=. scalability.jl                # the large-feeder sweep
```

`generate_results.jl` runs all six scripts and writes
`docs/src/assets/results/threephase/`, which is what the
[three-phase section of the tutorial](https://epsrlab-ub.github.io/SmartInverter-3P-DOPF.jl/dev/tutorial_voltvar/#Three-phases)
is built from. `TP_HOSTS=ivacopf` or `TP_HOSTS=lindist3flow` regenerates one family only.

The six scripts share their skeleton **verbatim**: data, PV placement, inverter
constraints, objective, verification, results and figures are byte-identical. Within a
host they differ only in the fenced `DROOP BLOCK` (and Heaviside swaps Gurobi for Ipopt);
across hosts they differ only in the network model. Diff any two to see exactly what one
choice costs.

Gurobi is [free for academic use](https://www.gurobi.com/academia/academic-program-and-licenses/);
Ipopt is open source, so the Heaviside scripts need no licence at all.

## Six runs, one curve

Same feeder, same inverters, same objective:

| host | encoding | class | solver | variables | binaries | constraints | solve | curtailed | max \|q − q_curve\| |
|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|
| LinDist3Flow | Lambda / SOS2 | MILP | Gurobi | 183,744 | 5,760 | 218,304 | 5.1 s | 46.32 kWh | 8.4e-07 |
| LinDist3Flow | Big-M | MILP | Gurobi | 179,136 | 5,760 | 225,216 | 5.4 s | 46.32 kWh | 6.1e-15 |
| LinDist3Flow | Heaviside | NLP | Ipopt | 171,072 | 0 | 207,936 | 33.6 s | 46.32 kWh | 1.1e-15 |
| IVACOPF | Lambda / SOS2 | MILP | Gurobi | 407,232 | 5,760 | 441,792 | 288.5 s (3 passes) | 42.69 kWh | 1.4e-17 |
| IVACOPF | Big-M | MILP | Gurobi | 402,624 | 5,760 | 448,704 | 136.1 s (3 passes) | 42.69 kWh | 4.4e-16 |
| IVACOPF | Heaviside | NLP | Ipopt | 394,560 | 0 | 431,424 | 142.8 s (3 passes) | 42.69 kWh | 1.2e-15 |

Within each host the three encodings agree to the digit on every physical quantity, and
every dispatch point lies on the droop curve of its own class. **That is a property of the
encodings.** Across the two hosts they do not agree, and that is a property of the network
model, which the audit below settles.

## Which host to believe

Both hosts put the inverters exactly on their curves *within their own model*. Only an
exact power flow can say whether that model is the right one. Each script re-solves its own
dispatch with a full three-phase backward/forward sweep and reports the gap:

| host | its own v vs the true AC v | droop residual at the **true** voltage | steps outside [0.95, 1.05] | true voltage range |
|:--|--:|--:|--:|:--|
| LinDist3Flow | 1.3e-03 p.u. | 8.0e-03 p.u. | 0 | 0.9922 – 1.0195 p.u. |
| IVACOPF | 1.9e-11 p.u. | 2.8e-11 p.u. | 0 | 0.9922 – 1.0204 p.u. |

LinDist3Flow's voltages are off by about a milli-per-unit, and because the droop is
steep, a full swing of q̄ across 0.02 p.u. on the upper segment, that turns into a
reactive-power error of a few per cent of inverter rating. The inverters would not deliver
the VArs it dispatched. IVACOPF reproduces the exact AC solution to round-off, so its
dispatch survives the transfer intact.

Neither host violates a voltage limit on this feeder, so the disagreement here is confined
to the droop residual, unlike the single-phase 33-bus case, where LinDistFlow additionally
pushes 17 time steps below the lower limit.

**Use LinDist3Flow for a fast first look; use IVACOPF for anything quantitative.**

## What the linearisation loop does

IVACOPF's two nonlinear relations, the `v·I` power balance and the voltage magnitude,
are expanded about the previous iterate and refreshed until the paper's own error metrics
clear a tolerance. Those metrics are measured against the **true** nonlinear relations, not
against the model's internal residual, which is what makes the converged point a genuine
power-flow solution:

| pass | solve | objective (p.u. curtailed) | MAPB | MRPB | MVM |
|--:|--:|--:|--:|--:|--:|
| 1 | 98.1 s | 1.800492 | 5.1e-03 | 1.6e-03 | 9.3e-04 |
| 2 | 88.3 s | 1.707866 | 4.7e-06 | 1.6e-06 | 2.8e-07 |
| 3 | 102.0 s | 1.707783 | 1.0e-09 | 3.3e-12 | 7.2e-13 |

- **MAPB / MRPB**: max absolute active / reactive power-balance error, eqs (18)–(19)
- **MVM**: max voltage-magnitude error, eq (20)

Three passes, roughly three orders of magnitude per pass. By default the loop starts from
an exact sweep at full PV and zero VArs rather than the paper's flat start
(`TP_WARMSTART=flat` for the latter), which costs a few seconds and saves a pass or two.

## What it does

| | |
|:--|:--|
| feeder | `network_5_Feeder_2`, a real ENWL LV feeder, Kron-reduced to three wires |
| size | 194 buses, 193 lines, 489 m, 415/240 V |
| loads | 18 single-phase, split **4 / 5 / 9** across phases (3.7 / 6.0 / 7.7 kW), genuinely unbalanced |
| PV | **12 smart inverters in 4 size classes**, 84 kW total, at the electrically farthest load buses, 4 per phase |
| horizon | 24 h at 15-minute resolution (96 steps), residential load shape, clear-sky irradiance |
| hosts | LinDist3Flow (phase-coupled linear branch flow) and IVACOPF (three-phase current-voltage AC-OPF) |
| droop | Lambda / SOS2, Big-M, or Heaviside on the IEEE 1547 curve |
| objective | minimise total PV curtailment |

The point of the feeder choice is the unbalance. Loads are single-phase and unevenly
distributed, so the three phases genuinely diverge, which is the only reason to model
three phases at all.

## The twelve smart inverters

Four classes. Since `q̄ = S_max`, each class follows a **different droop curve**: same
breakpoint voltages, four saturation levels, which is what makes the dispatch figure
worth looking at.

| class | P rated | S max | q̄ | q̄ (p.u.) | sites |
|:--|--:|--:|--:|--:|--:|
| A | 3 kW | 3.30 kVA | 3.30 kvar | 0.0330 | 3 |
| B | 5 kW | 5.50 kVA | 5.50 kvar | 0.0550 | 3 |
| C | 8 kW | 8.80 kVA | 8.80 kvar | 0.0880 | 3 |
| D | 12 kW | 13.20 kVA | 13.20 kvar | 0.1320 | 3 |

Each phase gets one inverter of each class, and the class order is rotated by phase, so
size is confounded with neither phase nor distance from the substation:

| phase | bus 1 | bus 2 | bus 3 | bus 4 |
|:-:|:--|:--|:--|:--|
| 1 | 184 (A) | 74 (B) | 73 (C) | 45 (D) |
| 2 | 188 (B) | 157 (C) | 153 (D) | 145 (A) |
| 3 | 193 (C) | 179 (D) | 142 (A) | 149 (B) |

## Results

```
                          LinDist3Flow        IVACOPF
PV energy available          476.9 kWh      476.9 kWh
PV curtailment               46.32 kWh      42.69 kWh
network losses             not modelled      14.61 kWh
voltage range          0.9923 – 1.0207   0.9922 – 1.0204 p.u.
peak line loading          no I variable      20.4 % of i_max
```

Each script writes three figures, suffixed by host and encoding: for example
`droop_dispatch_3ph.png` (LinDist3Flow + Lambda), `droop_dispatch_3ph_bigm.png`,
`droop_dispatch_3ph_iva.png` (IVACOPF + Lambda) and `droop_dispatch_3ph_iva_heaviside.png`,
plus `network_schematic.png` from `plot_network.jl`.

## The formulations

### LinDist3Flow

For a line with 3×3 phase impedance `Z`, the phase-coupled voltage drop is

```
w_j = w_i − (aR·P_ij + aX·Q_ij),     aR[φ,ψ] =  2·Re(α^{ψ−φ}·Z[φ,ψ])
                                      aX[φ,ψ] =  2·Im(α^{ψ−φ}·Z[φ,ψ]),   α = e^{−j2π/3}
```

derived from `|V_j|² = |V_i − Z·I|²` with the quadratic term dropped and voltages assumed
near-balanced. The ±√3 cross terms in the published form of these matrices are just that
rotation written out. The code works in magnitude rather than squared magnitude
(`w_j − w_i ≈ 2·(v_j − v_i)` near nominal), matching the single-phase version, so the
droop breakpoints stay in ordinary p.u. voltage.

Two sanity checks are worth keeping in mind: for a single phase `α⁰ = 1` gives `aR = 2r`,
`aX = 2x`, recovering `w_j = w_i − 2(rP + xQ)`; and for diagonal `Z` the phases decouple
into three independent LinDistFlows.

`aR` and `aX` are minus the `M^P`, `M^Q` of eqs. (22)-(23) of the reference; the drop is
its eq. (21) and the lossless power balance is its eq. (20). Both follow from its exact
Dist3Flow eqs. (14)-(17) under assumptions A1 (constant inter-phase voltage ratios) and
A2 (constant loss terms).

Reference: M. D. Sankur, R. Dobbe, E. Stewart, D. S. Callaway, and D. B. Arnold,
*A linearized power flow model for optimization in unbalanced distribution systems*,
arXiv:1606.04492, 2016, [arXiv:1606.04492](https://arxiv.org/abs/1606.04492).

### IVACOPF

The network in rectangular current-voltage coordinates. Ohm's law and KCL are then
**exactly linear**, mutual coupling and all, with no rotation operator, no transposition
assumption and no near-balance:

```
(4)–(5)   v_n^{re,φ} − v_m^{re,φ} = Σ_p ( R[φ,p]·I^{re,p} − X[φ,p]·I^{im,p} )
          v_n^{im,φ} − v_m^{im,φ} = Σ_p ( R[φ,p]·I^{im,p} + X[φ,p]·I^{re,p} )

(7)–(8)   I_n^{re,φ} = Σ_out I^{re,φ} − Σ_in I^{re,φ}          (and likewise imaginary)
```

Only two relations are nonlinear, and both sit at the **buses** rather than along the
lines, which is the structural point of the formulation:

```
(9)–(10)  p_n^φ = v^{re,φ}·I^{re,φ} + v^{im,φ}·I^{im,φ}
          q_n^φ = v^{im,φ}·I^{re,φ} − v^{re,φ}·I^{im,φ}
(12)      v_n^φ = sqrt( (v^{re,φ})² + (v^{im,φ})² )
```

Each is replaced by its first-order Taylor expansion about the previous iterate, (15)–(16)
and (17), and the loop repeats until MAPB, MRPB and MVM clear `TP_TOL`. Because the line
current is a decision variable, the thermal limit (13) is writable directly; it is offered
as a polygon inscribing the circle so the model stays an MILP, and left off by default
(`TP_IMAXSEG=0`) because peak flow here is about a fifth of the conductor rating.

Equation numbers follow Soltani, Khorsand & Ma, *Current–Voltage Unbalanced Distribution AC
Optimal Power Flow for Advanced Distribution Management System Applications*, IEEE Open
Journal of Industry Applications 5, 2024,
[doi:10.1109/OJIA.2024.3367547](https://doi.org/10.1109/OJIA.2024.3367547).

## The droop block is unchanged

Every droop block here is the same as in the single-phase code, and the same in both hosts,
because each only ever needs one thing: a voltage variable at the inverter's own bus and
phase.

```julia
vpv(i, t) = v[PV[i].bus, PV[i].phase, t]   # the whole three-phase interface
```

Nothing in any of the three encodings knows how many phases the network has, that the fleet
is a mix of four sizes, or which host it is sitting in. In LinDist3Flow `v` is a decision
variable directly; in IVACOPF it is the linearised magnitude (17). That is why six scripts
can share one skeleton, and why the results agree within each host to the digit.

## Things you may want to change

At the top of any script:

- `N_PV_PER_PHASE`: sites per phase. With four classes, four per phase gives each phase
  one of each; other counts still cycle through the classes.
- `PV_CLASSES`: the class names and ratings. Add or remove entries freely; the placement,
  the droop curves and the figure all follow. Setting every class to the same rating
  collapses the four curves back to one.
- `VBP`, `QSHAPE`: the droop curve. `QSHAPE` is a free vector, so asymmetric curves (for
  example AS/NZS 4777.2's −0.6 / +0.44) work without code changes.
- `VLIM`, `SBASE_KVA`, `MIP_GAP`.
- IVACOPF only: `MAX_ITER`, `TOL`, `WARMSTART`, `IMAX_SEG`.

All of these also have environment overrides, so the same model runs on a different feeder
or horizon without editing anything:

| variable | default | meaning |
|:--|:--|:--|
| `TP_CASE` | `network_5_Feeder_2` | ENWL feeder to load |
| `TP_STEPS` | `96` | time steps in the day |
| `TP_NPV` | `4` | smart inverters per phase |
| `TP_WARMSTART` | `sweep` | IVACOPF only: `flat` for the paper's flat start |
| `TP_TOL` | `1e-6` | IVACOPF only: stop tolerance on max(MAPB, MRPB, MVM) |
| `TP_MAXITER` | `15` | IVACOPF only: pass limit |
| `TP_IMAXSEG` | `0` | IVACOPF only: sides of the polygon enforcing (13); 0 disables it |

## A trap worth flagging

The validation sweep, and IVACOPF's slack reference and flat start, must use a properly
rotated three-phase set (1∠0°, 1∠−120°, 1∠+120°). Seeding all three phases at 1∠0° is
silent and expensive: the mutual impedance terms then add instead of largely cancelling,
and the sweep reports roughly **twice** the true voltage deviation, which looks exactly
like the linear host being badly wrong. It cost me a detour; the comments in `sweep_state()`
and at `V0` mark the spot.

## Data and licence

| file | what |
|:--|:--|
| `data/network_5_Feeder_2.bmopf.json` | the feeder, unmodified, in BMOPF JSON format |
| `data/network_5_Feeder_2_report.md` | its network summary, as published |
| `data/load_profile_residential_15min.json` | residential load shape, 96 steps |
| `data/solar_profile.json` | clear-sky irradiance, 96 steps |
| `data/network_17_Feeder_6.bmopf.json` | the 3856-bus feeder used for the scalability check |

## Scalability

`scalability.jl` runs the three LinDist3Flow encodings on both feeders at the full 96-step
horizon and, for any run that does not finish, again at a reduced horizon, then writes
`docs/src/assets/results/threephase/scalability.json`. It sweeps the linear host because
that is what isolates the *encodings'* scaling: one solve each, no outer loop. IVACOPF
multiplies every row by its pass count on top of a larger model per pass; its binary counts
are identical.

Any script can be pointed at the larger feeder directly:

```bash
TP_CASE=network_17_Feeder_6 julia --project=. LinDist3Flow_Lambda.jl
TP_CASE=network_17_Feeder_6 TP_STEPS=24 julia --project=. IVACOPF3Ph_Lambda.jl
```

`network_17_Feeder_6` comes from [PMDlab.jl](https://github.com/hei06j/PMDlab.jl)
(`data/three-wire/network_17/Feeder_6`), converted from OpenDSS to the same JSON schema as
the other feeder. Same ENWL LVNS origin, same CC BY 4.0 licence, same Kron reduction.

The feeder comes from
[BMOPFDraftData](https://github.com/frederikgeth/BMOPFDraftData) (Frederik Geth),
`output/ENWLvariants/Three-wire-Kron-reduced/`, derived from the CSIRO four-wire LV
dataset ([10.25919/jaae-vc35](https://doi.org/10.25919/jaae-vc35)), itself derived from
Electricity North West's LV Network Solutions data.

**Licence: CC BY 4.0**, commercial use permitted, attribution required. Note that the
per-case `meta.license` stamp present on the curated `benchmarks/` cases is *absent* on
these `output/` files; the licence comes from the repository's README table.
