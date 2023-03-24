# Numeric Settings
In each section, for each type of simulation, there are numerical hint and explation of how the boundaries and data are managed. 

!!! info "numeric info" 
    In general, the body force should be zero. The only reasonable execption is for the periodic channel.

!!! info "numeric info" 
    In general is better to set a reference speed `u_in = 1` and the desired reynolds, the viscosity is then computed ad `1/Reynods`. 

!!! warning "numeric waning" 
    In `MPI` context the Additive Shwarz Method is shown to better scale than the Geometric Algebraic MultiGrid.

!!! info "numeric info" 
    It is better to avoid using `:julia` as a solver because it uses a direct solver which is slow for large problems.
