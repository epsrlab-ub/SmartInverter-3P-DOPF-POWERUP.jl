# =====================================================================================
#  Experiment 2: What changes when voltage becomes a JuMP decision variable?
#
#  In Experiment 1, v was a known numerical value.
#  Here, v is a JuMP decision variable whose value will be chosen by an optimizer.
#
#  The goal is to discover why a normal Julia if-else can no longer be used directly.
# =====================================================================================

using JuMP


# -------------------------------------------------------------------------------------
# IEEE 1547 Volt-VAr curve parameters
# -------------------------------------------------------------------------------------

const VBP  = [0.88, 0.90, 0.97, 1.00, 1.02, 1.10]
const QBAR = 1.0

const A1 = -QBAR / (VBP[3] - VBP[2])
const A2 = -QBAR / (VBP[5] - VBP[4])


# -------------------------------------------------------------------------------------
# Same Volt-VAr droop function used in Experiment 1
# -------------------------------------------------------------------------------------

function q_droop(v)

    if v <= VBP[2]

        return QBAR

    elseif v <= VBP[3]

        return A1 * (v - VBP[3])

    elseif v <= VBP[4]

        return 0.0

    elseif v <= VBP[5]

        return A2 * (v - VBP[4])

    else

        return -QBAR

    end

end


# -------------------------------------------------------------------------------------
# Helper function for workshop pauses
# -------------------------------------------------------------------------------------

function pause()
    println()
    print("Press ENTER to continue...")
    readline()
    println()
end


# =====================================================================================
# STEP 1: Create a JuMP model
# =====================================================================================

println("\n============================================================")
println("  Experiment 2: Let the optimizer choose the voltage")
println("============================================================\n")

println("Instead of giving v a numerical value ourselves,")
println("we will ask an optimizer to choose v for us.")

pause()


model = Model()

@variable(model, 0.90 <= v <= 1.10)
@variable(model, -QBAR <= q <= QBAR)


# -------------------------------------------------------------------------------------
# Inspect the JuMP variable
# -------------------------------------------------------------------------------------

println("We create the voltage using:")
println()
println("    @variable(model, 0.90 <= v <= 1.10)")
println()

println("So v is no longer a numerical value.")

println("Rather, we want the optimizer to determine the numerical value of v.")

pause()


# =====================================================================================
# STEP 2: Let the participant define a grid equation
# =====================================================================================

println("============================================================")
println("  First, define a simple grid equation")
println("============================================================\n")

println("Before adding the Volt-VAr controller, we need a simple")
println("relationship describing how reactive power affects voltage.")
println()

println("Assume a linear grid model of the form:")
println()
println("    v = v₀ + kq")
println()

println("where:")
println()
println("    v₀ = voltage when q = 0")
println("    k  = voltage sensitivity to reactive power")
println()

println("You will choose the two parameters.\n")


print("Choose v₀ [p.u.]: ")
v0 = parse(Float64, readline())

print("Choose k [p.u./p.u.]: ")
k = parse(Float64, readline())


println()
println("You defined the grid equation:")
println()
println("    v = $(v0) + $(k)q")
println()

print("Press ENTER to add your equation to the JuMP model...")
readline()


@constraint(model, v == v0 + k * q)


println()
println("SUCCESS!")
println()
println("JuMP accepted your algebraic equation:")
println()
println("    v == $(v0) + $(k)q")
println()

println("Even though v and q are unknown decision variables,")
println("JuMP can represent an algebraic relationship between them.")

pause()


# =====================================================================================
# STEP 3: Introduce the Volt-VAr relationship
# =====================================================================================

println("============================================================")
println("  Now add the Volt-VAr controller")
println("============================================================\n")

println("We now have the GRID relationship:")
println()
println("    v = $(v0) + $(k)q")
println()

println("But the inverter must also follow the Volt-VAr controller:")
println()
println("    q = q_droop(v)")
println()

println("The operating point must satisfy BOTH relationships:")
println()
println("    Grid:        v = v₀ + kq")
println("    Controller:  q = q_droop(v)")
println()

pause()


# =====================================================================================
# STEP 4: The key experiment
# =====================================================================================

println("Now for the key experiment.")
println()

println("We will try to use the SAME q_droop function")
println("that worked in Experiment 1.")
println()

println("We want to add:")
println()
println("    @constraint(model, q == q_droop(v))")
println()

println("Remember that q_droop contains ordinary Julia logic:")
println()
println("    if v <= 0.90")
println("        ...")
println("    elseif v <= 0.97")
println("        ...")
println("    elseif v <= 1.00")
println("        ...")
println("    end")
println()


println("    q == q_droop(v)")
println()


# =====================================================================================
# STEP 5: Deliberately trigger the failure
# =====================================================================================

try

    @constraint(model, q == q_droop(v))

    println()
    println("The constraint was built successfully.")

catch err

    println()
    println("============================================================")
    println("  FAILED")
    println("============================================================\n")

    println("JuMP could not build the Volt-VAr constraint.")
    println()

    println("The key part of the error message is:")
    println()
    println("    ", first(split(sprint(showerror, err), '\n')))

    pause()


    # =================================================================================
    # STEP 6: Explain why
    # =================================================================================

    println("============================================================")
    println("  Why did it fail?")
    println("============================================================\n")

    println("The grid equation worked:")
    println()
    println("    v = $(v0) + $(k)q")
    println()
    println("because it is an algebraic relationship that JuMP")
    println("can pass to an optimization solver.")
    println()

    println("But q_droop(v) contains if-else statements.")
    println()
    println("For example:")
    println()
    println("    if v <= 0.90")
    println()

    println("At this point, v does not yet have a numerical value.")
    println()
    println("    typeof(v) = ", typeof(v))
    println()
    println("The optimizer is supposed to determine v only AFTER")
    println("the complete optimization model has been constructed.")
    println()

    println("Therefore Julia cannot decide which branch of")
    println("q_droop(v) should execute.")

    println("The issue is that:")
    println()
    println("    if-else  = programming logic")
    println()
    println("while an optimization solver requires")
    println()
    println("    algebraic constraints")
    println()

    println("So the Volt-VAr curve must be reformulated into")
    println("an optimization-compatible mathematical representation.")
    println()

    println("Next, we will look at formulations such as:")
    println()
    println("    1. Big-M")
    println("    2. Lambda / SOS2")
    println("    3. Heaviside")
    println()

end
  
print("Press ENTER to exit...")
readline()
