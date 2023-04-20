using SyntheticEddyMethod

function turbulence_box()
    Δy = 0.05 / 20
    Δz = 0.2 / 10
    σ = 2 * maximum([Δy, Δz])
    σ = 0.005 #eddy dimensions, the same in all the directions
    y = 0:0.01:0.07
    # z = -0.5:Δz:0.5
    z = -0.1:0.1:0.1
    
    Vboxinfo = VirtualBox(collect(y), collect(z), σ)

    # N = Vboxinfo.N #you can override it 
    return Vboxinfo
end

