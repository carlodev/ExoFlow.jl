function stretching_y_function(x)
   gamma1 = 2.5
   -tanh.(gamma1 .* (x)) ./ tanh.(gamma1)
end




function mesh_channel(params)
   
   function stretching(x::Point)
      m = zeros(length(x))
      m[1] = x[1]
      m[2] = stretching_y_function(x[2])
      if length(x) > 2
         if params[:periodic]
            m[3] = x[3]
         else
            m[3] = stretching_y_function(x[3])
   
         end
   
      end
      Point(m)
   end
   
   """
   mesh_channel() generate a mesh for a channel; 
   Periodic boundaries in dimensions 1 and 3 if periodic is set
   In dimensions 1 and 3 equally spaced
   In dimension 2 function distributed
   #Arguments
   - D::Integer number of dimensions (2 or 3)
   - N::Integer numer of cells in each dimension, deault value N=32
   - parts :: if distributed
   - printmodel::Boolean if true create vtk file pf the model
   """


   #D = 2 # Number of spatial dimensions
   #N = 32 # Partition (i.e., number of cells per space dimension)

   Lx = params[:Lx] 
   Ly = params[:Ly] 
   Lz = params[:Lz] 

   #generalization of the partition, if needed in the future each dimension can be split on a different number of nodes
   nx = params[:N]
   ny = params[:N]
   nz = params[:N]


   if params[:D] == 2
      domain = (0, Lx, -Ly / 2, Ly / 2)
      partition = (nx, ny)
      periodic_tuple = (params[:periodic], false)
   elseif params[:D] == 3
      domain = (0, Lx, -Ly / 2, Ly / 2, -Lz / 2, Lz / 2)
      partition = (nx, ny, nz)
      periodic_tuple = (params[:periodic], false, params[:periodic])
   end
   model_name = "model$(params[:D])d"

   if params[:parts] === nothing
      model = CartesianDiscreteModel(domain, partition, map=stretching, isperiodic=periodic_tuple)
   else
      model = CartesianDiscreteModel(params[:parts], domain, partition, map=stretching, isperiodic=periodic_tuple)
   end
   
   
   if params[:periodic]
      if params[:D] ==2
         tag_coordinate = Point(Lx/2, 0.0)
      else
         tag_coordinate = Point(Lx/2, 0.0, 0.0)
      end
      model = add_centre_tag!(model, tag_coordinate)
   end
   
   
   return model
end



