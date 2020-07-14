using Documenter, ADerrors

makedocs(modules=[ADerrors], doctest=true,
         pages = [
             "ADerrors" => "index.md",
             "Contents" => "toc.md"
             ], 
         sitename = "ADerrors.jl",
         repo = "https://gitlab.ift.uam-csic.es/alberto/aderrors.jl")

