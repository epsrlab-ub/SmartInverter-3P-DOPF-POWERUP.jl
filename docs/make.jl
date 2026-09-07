using Documenter
using SmartInverterDOPF

# The tutorial is built from precomputed results committed under
# docs/src/assets/results/threephase/, so no optimisation solver is needed here: the
# example blocks only read JSON and draw figures. Regenerate those results with
#     julia --project=examples/three_phase examples/three_phase/generate_results.jl
#     julia --project=examples/three_phase examples/three_phase/scalability.jl
ENV["GKSwstype"] = "100"          # headless GR backend

DocMeta.setdocmeta!(SmartInverterDOPF, :DocTestSetup,
                    :(using SmartInverterDOPF); recursive = true)

makedocs(
    modules  = [SmartInverterDOPF],
    sitename = "SmartInverter-3P-DOPF-POWERUP.jl",
    authors  = "Adedoyin Inaolaji",
    repo     = Remotes.GitHub("epsrlab-ub", "SmartInverter-3P-DOPF-POWERUP.jl"),
    # run @example blocks with docs/src as the working directory, so they can read
    # assets/results/*.json by relative path
    workdir  = joinpath(@__DIR__, "src"),
    format = Documenter.HTML(
        prettyurls  = get(ENV, "CI", nothing) == "true",
        canonical   = "https://epsrlab-ub.github.io/SmartInverter-3P-DOPF-POWERUP.jl",
        mathengine  = Documenter.KaTeX(),
        sidebar_sitename = false,
        assets      = String["assets/custom.css"],
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial_voltvar.md",
        "API reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/epsrlab-ub/SmartInverter-3P-DOPF-POWERUP.jl.git",
    devbranch = "main",
    push_preview = false,
)
