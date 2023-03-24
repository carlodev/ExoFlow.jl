# OK
function petsc_options(prec::Symbol)
    if prec == :snesgamg
    options = "-snes_type newtonls -snes_linesearch_type basic -snes_linesearch_damping 1.0 -snes_rtol 1.0e-8 -snes_atol 0 -snes_monitor  -snes_max_it 10 \
    -ksp_type fgmres -pc_fieldsplit_type symmetric_multiplicative -ksp_converged_reason -ksp_max_it 50 -ksp_rtol 1e-3 -ksp_atol 1e-5 \
    -fieldsplit_vel_pc_type gamg -fieldsplit_vel_ksp_type gmres -fieldsplit_vel_ksp_converged_reason \
    -fieldsplit_pres_pc_type gamg -fieldsplit_pres_pc_gamg_aggressive_coarsening 0.08 -fieldsplit_pres_pc_gamg_agg_nsmooths 10 -fieldsplit_pres_ksp_type gmres -fieldsplit_pres_ksp_converged_reason"
    
    elseif prec == :sneslu
        options = "-snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 1.0 -snes_rtol 1.0e-14 -snes_atol 0.0 -snes_monitor -pc_type asm -sub_pc_type lu  -ksp_type gmres -ksp_gmres_restart 30  -snes_converged_reason -ksp_converged_reason -ksp_error_if_not_converged true "
    
    elseif prec == :kspgamg
        options =  "-log_view -ksp_type fgmres -pc_fieldsplit_type symmetric_multiplicative -ksp_converged_reason -ksp_max_it 50 -ksp_rtol 1e-3 -ksp_atol 1e-5 \
         -fieldsplit_vel_pc_type gamg -fieldsplit_vel_ksp_type gmres -fieldsplit_vel_ksp_converged_reason \
         -fieldsplit_pres_pc_type gamg -fieldsplit_pres_pc_gamg_aggressive_coarsening 0.08 -fieldsplit_pres_pc_gamg_agg_nsmooths 10 -fieldsplit_pres_ksp_type gmres -fieldsplit_pres_ksp_converged_reason"
        
    
    elseif prec == :ksplu
        options = "-pc_type asm -sub_pc_type lu  -ksp_type gmres -ksp_gmres_restart 30  -snes_converged_reason -ksp_converged_reason -ksp_error_if_not_converged true"
    
    else
        error("petsc_options $prec not recognized as valid")
    end
    
    return options * " -log_view" * " -ksp_monitor_short"
end


function mykspsetup(ksp)
    pc = Ref{GridapPETSc.PETSC.PC}()
    @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[], pc)
    @check_error_code GridapPETSc.PETSC.PCSetType(pc[], GridapPETSc.PETSC.PCFIELDSPLIT)
    # @check_error_code GridapPETSc.PETSC.PCFieldSplitSetType(pc[],GridapPETSc.PETSC.PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE)
    @check_error_code GridapPETSc.PETSC.PCFieldSplitSetType(pc[], GridapPETSc.PETSC.PC_COMPOSITE_SCHUR)
    @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
    @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
end


function mysnessetup(snes)
    ksp = Ref{GridapPETSc.PETSC.KSP}()
    pc = Ref{GridapPETSc.PETSC.PC}()
    @check_error_code GridapPETSc.PETSC.SNESSetFromOptions(snes[])
    @check_error_code GridapPETSc.PETSC.SNESGetKSP(snes[], ksp)
    @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[], pc)
    @check_error_code GridapPETSc.PETSC.PCSetType(pc[], GridapPETSc.PETSC.PCFIELDSPLIT)
    # @check_error_code GridapPETSc.PETSC.PCFieldSplitSetType(pc[],GridapPETSc.PETSC.PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE)
    @check_error_code GridapPETSc.PETSC.PCFieldSplitSetType(pc[], GridapPETSc.PETSC.PC_COMPOSITE_SCHUR)

    @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
    @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
end



