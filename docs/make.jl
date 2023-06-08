using Hurst
using Documenter

ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(Hurst, :DocTestSetup, :(using Hurst); recursive = true)

makedocs(;
    modules = [Hurst],
    authors = "Joseph Bradley <josephbradley16@googlemail.com> and contributors",
    repo = "https://github.com/josephcbradley/Hurst.jl/blob/{commit}{path}#{line}",
    sitename = "Hurst.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://josephcbradley.github.io/Hurst.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md", 
            "Getting Started" => "getting_started.md", 
            "Troubleshooting" => "troubleshooting.md", 
            "Reference" => "ref.md"],
)

deploydocs(; repo = "github.com/josephcbradley/Hurst.jl", devbranch = "main")
