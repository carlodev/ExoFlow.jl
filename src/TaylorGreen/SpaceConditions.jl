"""
function add_centre_tag!(model, N, diameter)
    eps1 = 2 * diameter/ N / 10 # a value smaller than the side of one element, useful for extracting the centre point of the domain
    v1(v) = v[1]
    v2(v) = v[2]
    

    function is_centre(x)
    norm(v1.(x)) <= eps1 && norm(v2.(x)) <= eps1
    end
    map_parts(model.models) do model
        labels = get_face_labeling(model)
        model_nodes = DiscreteModel(Polytope{0}, model)
        cell_nodes_coords = get_cell_coordinates(model_nodes)
        cell_node_centre = collect1d(lazy_map(is_centre, cell_nodes_coords))
        cell_node = findall(cell_node_centre)
        new_entity = num_entities(labels) + 1
        for centre_point in cell_node
        labels.d_to_dface_to_entity[1][centre_point] = new_entity
        end
        add_tag!(labels, "centre", [new_entity])
    end

   return model
end
"""



#It is a unique case, with pressure that is time dependent and no boundary conditions on the velocity, it is periodic in all dimensions
function CreateTGSpaces(model, params, pa)
    
    model = add_centre_tag!(model, Point(0.0, 0.0)) #(0.0, 0.0) is the centre coordinate
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{params[:D], Float64}, params[:order])
    reffeₚ = ReferenceFE(lagrangian, Float64, params[:order])


    V = TestFESpace(model, reffeᵤ, conformity=:H1)
    U = TransientTrialFESpace(V)

    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags="centre")
    P = TransientTrialFESpace(Q, pa)

    Y = MultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])
    
    return V, Q, U, P, Y, X, model
end

#It is a unique case, with pressure that is time dependent and no boundary conditions on the velocity, it is periodic in all dimensions
function CreateTGSpaces_va(model, params, va)
    
    model = add_centre_tag!(model, Point(0.0, 0.0)) #(0.0, 0.0) is the centre coordinate
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{params[:D], Float64}, params[:order])
    reffeₚ = ReferenceFE(lagrangian, Float64, params[:order])


    V = TestFESpace(model, reffeᵤ, conformity=:H1,dirichlet_tags="centre")
    U = TransientTrialFESpace(V, va)

    Q = TestFESpace(model, reffeₚ, conformity=:H1)
    P = TransientTrialFESpace(Q)

    Y = MultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])
    
    return V, Q, U, P, Y, X, model
end