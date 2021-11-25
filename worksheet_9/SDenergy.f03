   program SDenergy
!
   use mqc_gaussian
   use iso_fortran_env
!
!    Variable Declarations ...
!
   implicit none
   character(len=:), allocatable :: fileName
   type(mqc_gaussian_unformatted_matrix_file):: fileInfo
   integer(kind=int64)::iOut=6,iPrint=2,IAlpha,IBeta,JAlpha,JBeta,NAlpha,NBeta
   logical :: UHF
   type(mqc_pscf_wavefunction)::wavefunction
   type(mqc_molecule_data)::moleculeInfo
   type(mqc_twoERIs)::eris,mo_ERIs
   type(mqc_scf_integral ):: mo_core_ham
   type(mqc_scalar)::Vnn,alpha_one_energy,beta_one_energy,alpha_two_energy_ss,&
       beta_two_energy_ss,two_energy_os, final_energy
!
!    Get the user-defined filename from the command line and then call the
!    routine that reads the Gaussian matrix file .
!
   call mqc_get_command_argument(1,fileName)
   call fileInfo %load(filename)
   call fileInfo %getESTObj('wavefunction',wavefunction)
   call fileInfo %getMolData(moleculeInfo)
   call fileInfo %get2ERIs('regular',eris)
!
!    Determine if the wavefunction that we read in was restricted or unrestricted .
!
   if(wavefunction%wf_type.eq.'U') then
       UHF = .true.
       write(*,*) 'Found unrestricted wavefunction'
   elseIf (wavefunction%wf_type.eq.'R') then
       UHF = .false.
       write(*,*) 'Found restricted wavefunction'
   else
       call mqc_error('Unsupported wavefunction type')
   endIf
   if (wavefunction%wf_complex) call mqc_error('Unsupported wavefunction type')
   NAlpha = wavefunction%NAlpha
   NBeta = wavefunction%NBeta
!
!   Compute the nuclear-nuclear repulsion energy.
!
   Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
   call Vnn%print(iOut,'Nuclear Repulsion Energy (au)')
!
!    Transform one and two-electron integrals to MO basis. These are the hij, Jij and Kij terms.
!
   if(iPrint.eq.1) write(iOut,*) 'Transforming MO integrals'
   mo_core_ham = matmul(transpose(wavefunction%MO_Coefficients),&
      matmul(wavefunction%core_Hamiltonian,&
      Wavefunction%MO_Coefficients))
   if(IPrint.ge.2) call mo_core_ham%print(iOut,'MO Basis Core Hamiltonian')
   call twoERI_trans(iOut,iPrint,wavefunction%MO_Coefficients,ERIs,mo_ERIs)
!
!   Compute Slater determinant energy.
!   For restricted , E = 2*sum_{i}hii + 2*sum_{ij}Jij-sum_{ij}Kij
!   For unrestricted , E = sum_{i}hiaia + sum_{i}hibib + sum_{ij}Jiaja + sum {ij}Jibjb +
!                         2*sum_{ij}Jiajb - sum {ij}Kiaja - sum_{ij}Kibjb
!
!   one-electron term
  alpha_one_energy = 0.0
  do IAlpha = 1, NAlpha
      alpha_one_energy = alpha_one_energy + mo_core_ham%at(IAlpha,IAlpha,'alpha')
  endDo
  if (.not.UHF) then
       alpha_one_energy = 2*alpha_one_energy
  else
       beta_one_energy = 0.0
       do IBeta = 1, NBeta
           beta_one_energy = beta_one_energy + mo_core_ham%at(IBeta,IBeta,'beta')
       endDo
   endIf
   call alpha_one_energy%print(iOut,'Alpha one electron energy: haa')
   If(UHF) call beta_one_energy%print(iOut,'Beta one electron energy: hbb')
!
!    two-electron same-spin term
   alpha_two_energy_ss = 0.0
   do IAlpha = 1, NAlpha
       do JAlpha = IAlpha+1, NAlpha
          alpha_two_energy_ss = alpha_two_energy_ss + &
             mo_ERIs%at(IAlpha,IAlpha,JAlpha,JAlpha,'alpha') - &
             mo_ERIs%at(IAlpha,JAlpha,JAlpha,IAlpha,'alpha')
      endDo
   endDo
   if(.not.UHF) then
       alpha_two_energy_ss = 2*alpha_two_energy_ss
   else
       beta_two_energy_ss = 0.0
       do IBeta = 1, NBeta
           do JBeta = IBeta+1, NBeta
              beta_two_energy_ss = beta_two_energy_ss + &
                  mo_ERIs%at(IBeta,IBeta,JBeta,JBeta,'beta') - &
                  mo_ERIs%at(IBeta,JBeta,JBeta,IBeta,'beta')
           endDo
        endDo
     endIf
     call alpha_two_energy_ss%print(iOut,'Alpha same spin two electron energy: Jaa - Kaa')
     if(UHF) call beta_two_energy_ss%print(iOut,'Beta same spin two electron energy: Jbb - Kbb')
!
!      two-electron opposite-spin term
     two_energy_os = 0.0
     do IAlpha = 1, NAlpha
        do IBeta = 1, NBeta
            if(UHF) then
               two_energy_os = two_energy_os + &
               mo_ERIs%at(IAlpha,IAlpha,IBeta,IBeta,'alphaBeta')
            else
               two_energy_os = two_energy_os + &
               mo_ERIs%at(IAlpha,IAlpha,IBeta,IBeta,'alpha')
           endIf
        endDo
     endDo
     call two_energy_os%print(iOut,'Opposite spin two electron energy (Coulomb only): Jab')
!
     if (.not.UHF) then
          final_energy = alpha_one_energy + alpha_two_energy_ss + two_energy_os
     else
          final_energy = alpha_one_energy + beta_one_energy + alpha_two_energy_ss + &
               beta_two_energy_ss + two_energy_os
     endIf
     call final_energy %print(iOut,'Total electronic energy')
!
     final_energy = Vnn + final_energy
     call final_energy %print(iOut,'Total electronic energy + nuclear repulsion energy')
!
     End Program SDEnergy

