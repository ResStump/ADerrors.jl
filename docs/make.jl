using Documenter, ADerrors

makedocs(modules=[ADerrors], doctest=true,
         pages = [
             "Getting Started" => "tutorial.md", 
             "API" => "api.md",
             "Contents" => "toc.md"
             ], 
         sitename = "ADerrors.jl",
         repo = "https://gitlab.ift.uam-csic.es/alberto/aderrors.jl")

