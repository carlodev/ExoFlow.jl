using SyntheticEddyMethod

function turbulence_box()
    Δy = 0.05 / 20
    Δz = 0.2 / 10
    σ = 2 * maximum([Δy, Δz])
    σ = 0.05 #eddy dimensions, the same in all the directions
    y = -1:Δy:1
    # z = -0.5:Δz:0.5
    z = -2/3*pi:Δz:2/3*pi
    x = 0:σ:2*pi
    Vboxinfo = VirtualBox(collect(x), collect(y), collect(z), σ)

    # N = Vboxinfo.N #you can override it 
    return Vboxinfo
end

