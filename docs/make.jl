using Documenter, ExoFlow

makedocs(
    sitename = "ExoFlow.jl",
    modules = [ExoFlow],
    pages = [
        "Introduction" => "index.md",
        "Usage" => "usage.md",
        "MPI Run" => "mpi.md",
        "Taylor Green" => "Case/taylorgreen.md",
        "Lid Driven Cavity Flow" => "Case/liddriven.md",
        "Airfoil" => "Cylinder/cylinder.md",
        "Channel" => "Case/channel.md",
        "Airfoil" => "Case/airfoil.md",
        "Reference" => "ref.md"
    ],
)

deploydocs(
    repo = "github.com/carlodev/ExoFlow.jl",
    push_preview = true,
)
