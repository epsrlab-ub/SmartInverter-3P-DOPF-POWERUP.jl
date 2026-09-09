# Why "if-else" cannot go straight into a solver

The IEEE 1547 Volt-VAr curve, written as a Julia `if`-`else` and handed to a solver.
That is the whole example: no network, no OPF, no data. It needs only JuMP — no solver,
because it never gets as far as solving anything.

## Run it

```bash
git clone https://github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl


julia --project=. examples\minimal\ifelse_numericworks.jl

julia --project=. examples\minimal\ifelse_variablefails.jl
```
