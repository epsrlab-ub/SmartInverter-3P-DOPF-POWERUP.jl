# Why "if-else" cannot go straight into a solver

The IEEE 1547 Volt-VAr curve, written as a Julia `if`-`else` and handed to a solver.
That is the whole example: no network, no OPF, no data. It needs only JuMP — no solver,
because it never gets as far as solving anything.

## Run it

```bash
git clone https://github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl
cd SmartInverter-3P-DOPF.jl/examples/minimal

julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. ifelse_fails.jl
```

## The code

```julia
const VBP  = [0.88, 0.90, 0.97, 1.00, 1.02, 1.10]   # IEEE 1547 breakpoints, p.u.
const QBAR = 1.0                                    # reactive capability, p.u.
const A1   = -QBAR / (VBP[3] - VBP[2])
const A2   = -QBAR / (VBP[5] - VBP[4])

function q_droop(v)
    if     v <= VBP[2];  return  QBAR               # saturated, full injection
    elseif v <= VBP[3];  return  A1 * (v - VBP[3])  # sloping down to zero
    elseif v <= VBP[4];  return  0.0                # dead-band
    elseif v <= VBP[5];  return  A2 * (v - VBP[4])  # sloping down again
    else                 return -QBAR               # saturated, full absorption
    end
end

model = Model()
@variable(model, 0.90 <= v <= 1.10)
@variable(model, -QBAR <= q <= QBAR)

@constraint(model, v == 1.035 + 0.03 * q)     # the grid, one linear sensitivity: accepted
@constraint(model, q == q_droop(v))           # the droop: this is the line that fails
```

## What happens

```
  q_droop(1.010) = -0.5   ← correct, on the steep segment

  @constraint(model, q == q_droop(v))

    Cannot evaluate `<=` between a variable and a number.
```

The function is a correct implementation of the curve — on a *number* it returns exactly
the right answer. The model still does not build.

`v` is a decision variable. The solver has not chosen it yet, so `v <= 0.90` has no truth
value to branch on and the body of the `if` never runs. Branching is not an algebraic
constraint, and rewriting the `if` does not change that.

## What to do instead

The curve has to reach the solver whole, as algebra. Three exact ways to write it, all of
them in this repository:

| encoding | idea | class |
|:--|:--|:--|
| **Big-M** | a binary per segment activates that segment's law and voltage range | MILP |
| **Lambda / SOS2** | the operating point is a blend of two adjacent breakpoints | MILP |
| **Heaviside** | segment masks built from unit steps, no integers at all | NLP |

All three return the same dispatch. See [`examples/single_phase`](../single_phase),
[`examples/three_phase`](../three_phase), and the
[tutorial](https://epsrlab-ub.github.io/SmartInverter-3P-DOPF.jl/dev/tutorial_voltvar/).
