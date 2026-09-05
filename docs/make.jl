using Documenter
using SmartInverterDOPF

# The tutorial is built from precomputed results committed under
# docs/src/assets/results/, so no optimisation solver is needed here: the example
# blocks only read JSON and draw figures. Regenerate those results with
#     julia --project=scripts scripts/generate_results.jl
ENV["GKSwstype"] = "100"          # headless GR backend

DocMeta.setdocmeta!(SmartInverterDOPF, :DocTestSetup,
                    :(using SmartInverterDOPF); recursive = true)

makedocs(
    modules  = [SmartInverterDOPF],
    sitename = "SmartInverterDOPF.jl",
    authors  = "Rahmat Emami Mirak",
    repo     = Remotes.GitHub("epsrlab-ub", "SmartInverter-3P-DOPF.jl"),
    # run @example blocks with docs/src as the working directory, so they can read
    # assets/results/*.json by relative path
    workdir  = joinpath(@__DIR__, "src"),
    format = Documenter.HTML(
        prettyurls  = get(ENV, "CI", nothing) == "true",
        canonical   = "https://epsrlab-ub.github.io/SmartInverter-3P-DOPF.jl",
        mathengine  = Documenter.KaTeX(),
        sidebar_sitename = false,
        assets      = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial_voltvar.md",
        "API reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/epsrlab-ub/SmartInverter-3P-DOPF.jl.git",
    devbranch = "main",
    push_preview = false,
)
