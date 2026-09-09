# =====================================================================================
#  Experiment 1: A Julia if-else works perfectly when voltage is a known number
#
#  Participants choose a voltage value, evaluate the IEEE 1547 Volt-VAr droop,
#  and see their operating point highlighted on the full droop curve.
# =====================================================================================

using Plots

# -------------------------------------------------------------------------------------
# IEEE 1547 Volt-VAr curve parameters
# -------------------------------------------------------------------------------------

const VBP  = [0.88, 0.90, 0.97, 1.00, 1.02, 1.10]   # voltage breakpoints, p.u.
const QBAR = 1.0                                     # reactive capability, p.u.

const A1 = -QBAR / (VBP[3] - VBP[2])                # first sloped segment
const A2 = -QBAR / (VBP[5] - VBP[4])                # second sloped segment


# -------------------------------------------------------------------------------------
# Volt-VAr droop function
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
# Participant interaction
# -------------------------------------------------------------------------------------

println("\n============================================================")
println("  Experiment 1: You choose the voltage")
println("============================================================\n")

println("The inverter follows an IEEE 1547 Volt-VAr curve.")
println("Choose any voltage between 0.88 and 1.10 p.u.\n")


# Keep asking until the participant enters a valid voltage
while true

    print("Enter voltage v [p.u.]: ")

    global v

    try

        v = parse(Float64, readline())

        if 0.88 <= v <= 1.10
            break
        else
            println("\nPlease enter a voltage between 0.88 and 1.10 p.u.\n")
        end

    catch

        println("\nPlease enter a valid numerical value, for example 1.01.\n")

    end

end


# -------------------------------------------------------------------------------------
# Evaluate the droop curve
# -------------------------------------------------------------------------------------

q = q_droop(v)


println("\n------------------------------------------------------------")
println("Result")
println("------------------------------------------------------------")

println("  v = ", round(v, digits=4), " p.u.")
println("  q = ", round(q, digits=4), " p.u.")


# Identify the part of the Volt-VAr curve
if v <= VBP[2]

    println("  Region: maximum reactive power injection")

elseif v <= VBP[3]

    println("  Region: reactive power injection slope")

elseif v <= VBP[4]

    println("  Region: deadband")

elseif v <= VBP[5]

    println("  Region: reactive power absorption slope")

else

    println("  Region: maximum reactive power absorption")

end


println("\nJulia can evaluate this if-else because v is already a known number.")
println("typeof(v) = ", typeof(v))


# -------------------------------------------------------------------------------------
# Generate the complete Volt-VAr curve
# -------------------------------------------------------------------------------------

v_curve = range(VBP[1], VBP[end], length=500)

q_curve = q_droop.(v_curve)


# -------------------------------------------------------------------------------------
# Plot the Volt-VAr curve
# -------------------------------------------------------------------------------------

p = plot(
    v_curve,
    q_curve,
    linewidth = 3,
    xlabel = "Voltage, v [p.u.]",
    ylabel = "Reactive Power, q [p.u.]",
    title = "IEEE 1547 Volt-VAr Droop Curve",
    label = "Volt-VAr curve",
    xlims = (0.875, 1.105),
    ylims = (-1.1, 1.1),
    legend = :bottomleft,
    grid = true,
)


# -------------------------------------------------------------------------------------
# Plot the Volt-VAr breakpoints
# -------------------------------------------------------------------------------------

scatter!(
    p,
    VBP,
    q_droop.(VBP),
    markersize = 5,
    label = "Breakpoints",
)



# -------------------------------------------------------------------------------------
# Highlight the selected operating point
# -------------------------------------------------------------------------------------

scatter!(
    p,
    [v],
    [q],
    markersize = 9,
    markerstrokewidth = 2,
    label = "Your point ($(round(v,digits=3)), $(round(q,digits=3)))",
)


# -------------------------------------------------------------------------------------
# Display the plot
# -------------------------------------------------------------------------------------

display(p)


println("\nThe highlighted point is the operating point corresponding to your voltage.")
println("\nPress ENTER when you are finished viewing the plot.")

readline()

# Close the popup window
Plots.GR.inline("closeall")
