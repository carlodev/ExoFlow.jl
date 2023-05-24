using SyntheticEddyMethod

function turbulence_box()
    Δy = 0.05 / 20
    Δz = 0.2 / 10
    σ = 2 * maximum([Δy, Δz])
    σ = 0.02 #eddy dimensions, the same in all the directions
    y = -3:0.1:3
    z = -0.1:Δz:0.1
    # z = -1:0.1:1
    
    Vboxinfo = VirtualBox(collect(y), collect(z), σ)

    # N = Vboxinfo.N #you can override it 
    return Vboxinfo
end

