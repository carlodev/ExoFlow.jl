function initialize_Mat(p_dofs, u_dofs; formulation =:segregated)
    Mat_Tuu = spzeros(Float64, Int32, u_dofs, u_dofs)
    Mat_Tpu =spzeros(Float64, Int32, p_dofs, u_dofs)
    Mat_Auu  = spzeros(Float64, Int32, u_dofs, u_dofs)
    Mat_Aup  = spzeros(Float64, Int32, u_dofs, p_dofs)
    Mat_Apu = spzeros(Float64, Int32, p_dofs, u_dofs)
    Mat_App = spzeros(Float64, Int32, p_dofs, p_dofs)
if formulation == :segregated
    Mat_ML = spzeros(Float64, Int32, u_dofs, u_dofs)
    Mat_inv_ML = spzeros(Float64, Int32, u_dofs, u_dofs)
    return Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, Mat_ML, Mat_inv_ML

elseif formulation ==:coupled
    Mat_Tpp = spzeros(Float64, Int32,p_dofs, p_dofs) 
    Mat_Tup = spzeros(Float64, Int32,u_dofs, p_dofs)
    return Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, Mat_Tpp,  Mat_Tup
else
    error("formulation $formulation not recognized as a valid input, use :segregated or :coupled")
end

end


function update_Mat!(Mat_Tuu::SparseMatrixCSC, 
    Mat_Tpu::SparseMatrixCSC, 
    Mat_Auu::SparseMatrixCSC, 
    Mat_Aup::SparseMatrixCSC,  
    Mat_Apu::SparseMatrixCSC,  
    Mat_App::SparseMatrixCSC, 
    u_adv, as; simplified = false)



cconv(u_adv, ∇u) = u_adv ⋅ (∇u)

Tuu(u, v) = ∫((v + τsu∘(u_adv, h) * (cconv ∘ (u_adv, ∇(v)))) ⊙ u)dΩ
Tpu(u, q) = ∫((τsu∘(u_adv, h)) * (∇(q)) ⊙ u)dΩ
Auu1(u, v) = ∫(ν* ∇(v) ⊙ ∇(u) + (cconv ∘ (u_adv, ∇(u))) ⋅ v  + ((τsu∘(u_adv, h)) *(cconv ∘ (u_adv, ∇(v)))) ⊙(cconv ∘ (u_adv, ∇(u))))dΩ

Auu2(u, v) = ∫(((τb∘(u_adv, h)) * (∇⋅v)) ⊙ (∇ ⋅ u) + 0.5 .*u_adv⋅(v + (τsu∘(u_adv, h))*(cconv ∘ (u_adv, ∇(v))))⋅(∇ ⋅ u))dΩ

Auu(u, v) = Auu1(u, v) + Auu2(u, v)

Aup(p, v) = ∫( - (∇ ⋅ v) * p + ((τsu∘(u_adv, h))*(cconv ∘ (u_adv, ∇(v)))) ⊙ ∇(p))dΩ

Apu(u, q) = ∫( q *(∇ ⋅ u)+  0.5 .* (τsu∘(u_adv, h))⋅(∇(q))⋅u_adv⋅(∇ ⋅ u) + (τsu∘(u_adv, h))* (∇(q)) ⊙ (cconv ∘ (u_adv, ∇(u))))dΩ

App(p, q) =  ∫(((τsu∘(u_adv, h)) * ∇(q)) ⊙ (∇(p)) )dΩ

Tuu_simplified(u,v) = ∫(v' ⊙ u)dΩ
Auu1_simplified(u, v) = ∫((ν + (τb∘(u_adv, h)) ) ⋅v' ⊙ u)dΩ
Auu2_simplified(u, v) = ∫((τb∘(u_adv, h)) ⋅v' ⊙ u)dΩ
Auu_simplified(u, v) = Auu1_simplified(u, v) + Auu2_simplified(u, v)

if simplified
    Mat_Tuu[:,:] = Gridap.FESpaces.assemble_matrix(Tuu_simplified, as.UV, as.U0, as.V)
    Mat_Auu[:,:] = Gridap.FESpaces.assemble_matrix(Auu_simplified, as.UV,as.U0, as.V)


else
 

    Mat_Tuu[:,:] = Gridap.FESpaces.assemble_matrix(Tuu, as.UV, as.U0, as.V)
    Mat_Auu[:,:] = Gridap.FESpaces.assemble_matrix(Auu1, as.UV, as.U0, as.V) +
                   Gridap.FESpaces.assemble_matrix(Auu2, as.UV, as.U0, as.V) 
    


end
    
    Mat_Tpu[:,:] = Gridap.FESpaces.assemble_matrix(Tpu, as.UQ, as.U0, as.Q)
    Mat_Aup[:,:] = Gridap.FESpaces.assemble_matrix(Aup, as.PV, as.P0, as.V)
    Mat_Apu[:,:] = Gridap.FESpaces.assemble_matrix(Apu, as.UQ, as.U0, as.Q)
    Mat_App[:,:] = Gridap.FESpaces.assemble_matrix(App, as.PQ, as.P0, as.Q)
   
    
end


function update_Mat!(Mat_Tuu::SparseMatrixCSC, 
    Mat_Tpu::SparseMatrixCSC, 
    Mat_Auu::SparseMatrixCSC, 
    Mat_Aup::SparseMatrixCSC,  
    Mat_Apu::SparseMatrixCSC,  
    Mat_App::SparseMatrixCSC, 
    u_adv, t::Real; simplified = false)



cconv(u_adv, ∇u) = u_adv ⋅ (∇u)

Tuu(u, v) = ∫((v + τsu∘(u_adv, h) ⋅ (cconv ∘ (u_adv, ∇(v))))' ⋅ u)dΩ
Tpu(u, q) = ∫((τsu∘(u_adv, h)) ⋅ (∇(q))' ⋅ u)dΩ
Auu1(u, v) = ∫(ν* ∇(u) ⊙ ∇(v) + (cconv ∘ (u_adv, ∇(u))) ⋅ v  + ((τsu∘(u_adv, h)) *(cconv ∘ (u_adv, ∇(v))))' ⋅(cconv ∘ (u_adv, ∇(u))))dΩ

Auu2(u, v) = ∫(((τb∘(u_adv, h)) ⋅ (∇⋅v)) ⋅ (∇ ⋅ u) + 0.5 .*u_adv⋅(v + (τsu∘(u_adv, h))⋅(cconv ∘ (u_adv, ∇(v))))⋅(∇ ⋅ u))dΩ

Auu(u, v) = Auu1(u, v) + Auu2(u, v)

Aup(p, v) = ∫( v'⋅∇(p)+ ((τsu∘(u_adv, h))⋅(cconv ∘ (u_adv, ∇(v))))' ⊙ ∇(p))dΩ

Apu(u, q) = ∫( q *(∇ ⋅ u)+  0.5 .* (τsu∘(u_adv, h))⋅(∇(q))⋅u_adv⋅(∇ ⋅ u) + (τsu∘(u_adv, h))⋅ (∇(q))' ⋅ (cconv ∘ (u_adv, ∇(u))))dΩ

App(p, q) =  ∫(((τsu∘(u_adv, h)) ⋅ ∇(q))' ⋅ (∇(p)) )dΩ

Tuu_simplified(u,v) = ∫(v' ⊙ u)dΩ
Auu1_simplified(u, v) = ∫((ν + (τb∘(u_adv, h)) ) ⋅v' ⊙ u)dΩ
Auu2_simplified(u, v) = ∫((τb∘(u_adv, h)) ⋅v' ⊙ u)dΩ
if simplified
    Mat_Tuu[:,:] = Gridap.FESpaces.assemble_matrix(Tuu_simplified, assemUV(t), U(t), V)
    Mat_Auu[:,:] = Gridap.FESpaces.assemble_matrix(Auu1_simplified, assemUV(t), U(t), V) +
                   Gridap.FESpaces.assemble_matrix(Auu2_simplified, assemUV(t), U(t), V)


else
 

    Mat_Tuu[:,:] = Gridap.FESpaces.assemble_matrix(Tuu, assemUV(t), U(t), V)
    Mat_Auu[:,:] = Gridap.FESpaces.assemble_matrix(Auu1, assemUV(t), U(t), V) +
                   Gridap.FESpaces.assemble_matrix(Auu2, assemUV(t), U(t), V) 
    


end
    
    Mat_Tpu[:,:] = Gridap.FESpaces.assemble_matrix(Tpu, assemUQ(t), U(t), Q)
    Mat_Aup[:,:] = Gridap.FESpaces.assemble_matrix(Aup, assemPV(t), P(t), V)
    Mat_Apu[:,:] = Gridap.FESpaces.assemble_matrix(Apu, assemUQ(t), U(t), Q)
    Mat_App[:,:] = Gridap.FESpaces.assemble_matrix(App, assemPQ(t), P(t), Q)
    
end

function inv_lump_vel_mass!(Mat_inv_ML::SparseMatrixCSC, Mat_ML::SparseMatrixCSC)
    inv_ML_vec = 1 ./ sum(Mat_ML, dims=2)[:,1]
        if !isempty(inv_ML_vec[inv_ML_vec.==Inf])
        error("The matrix ML can not be inverted because after lumping zero values are detected")
    end
    idx = Int32.(1:1:size(Mat_ML)[1])
    Mat_inv_ML[:, :] = SparseArrays.sparse(idx, idx, inv_ML_vec)
end


function initialize_vec(ph_init, uh_init)
    vec_pm = get_free_dof_values(ph_init)
    vec_um = get_free_dof_values(uh_init)
    
    p_dofs = length(vec_pm)
    u_dofs = length(vec_um)

    vec_am = spzeros(Float64, Int32, u_dofs)

    vec_sum_pm = spzeros(Float64, Int32, p_dofs)
    return vec_pm, vec_sum_pm, vec_um, vec_am
end



function update_vec_vec_um!(vec_vec_um::Matrix{Float64}, vec_um::Vector{Float64})
    vec_vec_um[:,2:end] = vec_vec_um[:,1:3]
    vec_vec_um[:,1] = vec_um
end



function  update_vec_um!(vec_um::Vector{Float64}, vec_vec_um::Matrix{Float64}, coeff::Vector{Float64})
    vec_um[:] = vec_vec_um * coeff
end