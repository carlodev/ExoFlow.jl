# Index

## Conservation
```@docs
conservations_equations
derivative_conservations_equations
```

## Linear Utilities
```@docs
create_ũ_vector
update_ũ
update_ũ_vector!
update_linear!
```

## SUPG
```@docs
h_param
τ
τb
SUPG
SUPG_lin
```

## VMS
```@docs
G_params
τm
τc
VMS
VMS_lin
```

## PETSc
```@docs
petsc_options
```

## Common procedures
```@docs
hf_gen!
add_centre_tag!
printmodel
creation_fe_spaces
creation_op
creation_initial_conditions
create_system_solver
creation_ode_parameters
creation_ode_solver
compute_solution_
compute_solution
iterate_solution
compute_solution_benchmark
solve_case

initialize_parameters
```


## Forces
```@docs
forces_domain
compute_forces
write_forces
```