push!(LOAD_PATH,"../src/")

using Documenter
using ExoFlow

makedocs(
    sitename = "ExoFlow.jl",
    modules = [ExoFlow],
    pages = [
        "Introduction" => "index.md",
        "Usage" => "usage.md",
        "MPI Run" => "mpi.md",
        "Numeric Settings" => "numericsettings.md",
        "Cases" =>[
        "Taylor Green" => "Cases/taylorgreen.md",
        "Lid Driven Cavity Flow" => "Cases/liddriven.md",
        "Cylinder" => "Cases/cylinder.md",
        "Channel" => "Cases/channel.md",
        "Airfoil" => "Cases/airfoil.md"],
        "Reference" => "ref.md"
    ],
)

deploydocs(
    repo = "github.com/carlodev/ExoFlow.jl",
    push_preview = true,
)
