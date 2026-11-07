using Documenter, HTSLib

DocMeta.setdocmeta!(HTSLib, :DocTestSetup, :(using HTSLib); recursive = true)

makedocs(
    modules = [HTSLib, HTSLib.Bindings],
    repo = Documenter.Remotes.GitHub("BioJulia", "HTSLib.jl"),
    sitename = "HTSLib.jl",
    doctest = true,
    pages = [
        "HTSLib" => "index.md",
        "Manual" => [
            "BGZF and tabix" => "tabix.md",
            "SAM, BAM and CRAM" => "alignments.md",
            "Positions and intervals" => "positions.md",
            "Records and iteration" => "records.md",
            "The GenomicFeatures extension" => "genomicfeatures.md",
            "Messages and errors" => "messages.md",
        ],
        "Reference" => "reference.md",
        "Low-level bindings" => "bindings.md",
    ],
    format = Documenter.HTML(edit_link = "develop"),
    authors = "Ciarán O'Mara",
    checkdocs = :public,
)

deploydocs(;
    repo = "github.com/BioJulia/HTSLib.jl.git",
    devbranch = "develop",
    push_preview = true,
    deps = nothing,
    make = nothing,
)
