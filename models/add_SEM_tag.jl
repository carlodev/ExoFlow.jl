
"add the SEM tag to the nodes in front of the airfoil. c: chord of the airfoil, the points are on a circle at 1 chord of distance from the leading edge. 
set a_tol for selecting only one lines of points, it depends on how fine the mesh is. "
function add_SEM_tag!(model; c=1.0, a_tol = 1e-1)

    v1(v) = v[1]
    v2(v) = v[2]
    v3(v) = v[3]

 
    function is_SEM(x)
        x = x[1]
        r2 = (v1.(x)).^2 .+ (v2.(x)).^2
        x[1][1]< -c/3 && isapprox(r2, c.^2, atol=a_tol)
    end

  
    function create_SEM_tag!(model::GridapDistributed.DistributedDiscreteModel)
        map_parts(model.models) do model
            create_SEM_tag!(model)
        end

    end

    function create_SEM_tag!(model)
        labels = get_face_labeling(model)
        model_nodes = DiscreteModel(Polytope{0}, model)
        cell_nodes_coords = get_cell_coordinates(model_nodes)
        cell_node_centre = collect1d(lazy_map(is_SEM, cell_nodes_coords))
        cell_node = findall(cell_node_centre)
        new_entity = num_entities(labels) + 1
        for centre_point in cell_node
            labels.d_to_dface_to_entity[1][centre_point] = new_entity
        end
        add_tag!(labels, "SEM", [new_entity])
    end


    create_SEM_tag!(model)

    return model
end