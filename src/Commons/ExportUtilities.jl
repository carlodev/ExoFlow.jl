"""
    conv_VectorValue(v::VectorValue)

Convert VectorValue (Gridap Type) to a Vector
"""
function conv_VectorValue(v::VectorValue)
    [v...]
end

"""
    get_dimension(vv::Vector)

It provides the dimension
"""
function get_dimension(vv::Vector)
    D = 0
    if length(vv) > 0
        D = length(vv[1])
    end
    return D

end

"""
    conv_to_df(vv::Vector)

Convert a Vector{Vector} to a DataFrame. It is used for export velocity field.
"""
function conv_to_df(vv::Vector)
    d = get_dimension(vv)
    n = length(vv)
    x = zeros(n)
    y = zeros(n)
    z = zeros(n)
    for (i, val) in enumerate(vv)
        x[i] = val[1]
        y[i] = val[2]
        if d == 3
            z[i] = val[3]
        end

    end
    df = DataFrame(x=x, y=y, z=z)
    return df
end


"""
    conv_to_df(vv::Vector)

Convert a Vector{Float64} to a DataFrame. It is used for export pressure field.
"""
function conv_to_df(vv::Vector{Float64})
    df = DataFrame(p=vv)
    return df
end

"""
    export_time_step(t::Float64, vv::Vector, fname::String, part::Int64)

Save .csv file, one for each processor
"""
function export_time_step(t::Float64, vv::Vector, fname::String, part::Int64)
    df = conv_to_df(vv)

    dir = joinpath("Results", "$(fname)_$t")
    mkpath(dir)
    filename = joinpath(dir, "$(fname)_$(part)_$t.csv")
    CSV.write(filename, df)
end



function export_nodes(parts, trian)
    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

    #export nodes
    local_unique_idx = map_parts(parts, trian.trians) do part, ttrian
        ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
        visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
        visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
        nodes_tri = unique(visgrid_)

        unique_idx = unique(i -> visgrid_[i], eachindex(visgrid_))
        export_time_step(0.0, nodes_tri, "nodes", part)

        return unique_idx
    end

    return local_unique_idx
end


function export_fields(params::Dict{Symbol,Any}, local_unique_idx::AbstractPData, tt::Float64, uh0, ph0)
    #get unique values in each processor

    @unpack force_params, parts = params

    if force_params !== nothing
        f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

        @unpack Γ,n_Γ = force_params
        n_Γ = - n_Γ #pointing from the body to the outside
        t_Γ = rotation∘n_Γ #extract tangent

        friction = (transpose(∇(uh0))⋅n_Γ) ⋅ t_Γ

        cellfields = Dict("ph" => ph0, "n_Γ" => n_Γ,"friction" => friction)

        fdat = GridapDistributed._prepare_fdata(Γ.trians, cellfields)

        map_parts(parts, Γ.trians, fdat, local_unique_idx) do part, ttrian, cf, unique_idx
            ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
            visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
            pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)
            for field in keys(cellfields)
            field_h = pdata[field][unique_idx]
            export_time_step(tt, field_h, field, part)
            end
        end
    end
end


function rotation(n::VectorValue{2,Float64})
    n1,n2 = [n...]
    VectorValue(-n2,n1)
end

function rotation(n::VectorValue{3,Float64})
    n1,n2,n3 = [n...]
    VectorValue(-n2,n1,n3)
end

function get_nodes(params::Dict{Symbol,Any})
    @unpack force_params,parts = params
    if force_params !== nothing
        @unpack Γ,n_Γ = force_params
        local_unique_idx =  export_nodes(parts, Γ)
    else
        local_unique_idx = nothing
    end

    return local_unique_idx
end

