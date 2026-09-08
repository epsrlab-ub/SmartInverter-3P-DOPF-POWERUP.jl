# =====================================================================================
#  Why "if-else" cannot go straight into a solver
#
#  The IEEE 1547 Volt-VAr curve, written as a Julia `if`-`else` and handed to a solver.

# =====================================================================================

using JuMP

const VBP  = [0.88, 0.90, 0.97, 1.00, 1.02, 1.10]   # IEEE 1547 breakpoints, p.u.
const QBAR = 1.0                                    # reactive capability, p.u.
const A1   = -QBAR / (VBP[3] - VBP[2])              # slope of the first sloped segment
const A2   = -QBAR / (VBP[5] - VBP[4])              # slope of the second

"The Volt-VAr droop, written the way anyone writes it first."
function q_droop(v)
    if     v <= VBP[2];  return  QBAR               # saturated, full injection
    elseif v <= VBP[3];  return  A1 * (v - VBP[3])  # sloping down to zero
    elseif v <= VBP[4];  return  0.0                # dead-band
    elseif v <= VBP[5];  return  A2 * (v - VBP[4])  # sloping down again
    else                 return -QBAR               # saturated, full absorption
    end
end

# On a number it is perfect. Nothing is wrong with this implementation of the curve.
println("\n  q_droop(1.010) = ", q_droop(1.010), "   ← correct, on the steep segment\n")

# Now ask a solver to pick the voltage instead.
model = Model()
@variable(model, 0.90 <= v <= 1.10)
@variable(model, -QBAR <= q <= QBAR)

@constraint(model, v == 1.035 + 0.03 * q)     # the grid, one linear sensitivity: accepted
try
    @constraint(model, q == q_droop(v))       # the droop: this is the line that fails
    println("  built without error — unexpected!")
catch err
    # Delete the try/catch and the script simply crashes here. JuMP's message runs to a
    # page; its first line is the whole story.
    println("  @constraint(model, q == q_droop(v))\n")
    println("    ", first(split(sprint(showerror, err), '\n')))
    println("""

  `v` is a decision variable, not a number. The solver has not chosen it yet, so
  `v <= 0.90` has no truth value to branch on and the body of the `if` never runs.
  Branching is not an algebraic constraint, and rewriting the `if` does not change that.

  The curve has to be handed to the solver whole, as algebra. That is what Big-M,
  Lambda/SOS2 and Heaviside each do, in examples/single_phase and examples/three_phase.
""")
end
