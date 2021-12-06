program hartreefock

use mqc_gaussian
use iso_fortran_env

implicit none
character(len=:), allocatable :: command,fileName
type(mqc_gaussian_unformatted_matrix_file)::matFile
integer(kind=int64)::iOut=output_unit,iIn=input_unit,iPrint=0,i,j,nElec,nAlpha,nBeta,nBasis, &
iter =1,max_iter=4096,nBasUse,multi
type(mqc_scalar)::conver,thresh
type(mqc_wavefunction)::wavefunc
type(mqc_molecule_data)::moleculeInfo
type(mqc_twoERIs),dimension(1)::eris
type(mqc_scalar)::Vnn,Energy,half
type(mqc_scf_integral):: Gmat,Fock,Xmat,old_density
!
!****TIMER****'
real :: start, finish,runtime
!**************
call cpu_time(start)
!************
!
Write(IOut,*)
write(IOut,*)
Write(IOut,*)
Write(Iout,*)
Write(IOut,*)
half = 0.5
thresh = 1.0e-8
!
 j = 1
do i=1,command_argument_count()
    if(i .ne.j ) cycle
    call mqc_get_command_argument(i,command)
    if(command.eq.'-f') then
!
!*      -f matrix_file      Input matrix file with initial set of molecular orbitals .
!*
    call mqc_get_command_argument(i+1,fileName)
    j = i+2
   else
       call mqc_error_A('Unrecognised input flag',6,'command',command)
   endIf
   deallocate(command)
endDo
!
!    Recover required data from matrix file .
!
call matFile%load(fileName)
call matFile%getESTObj('wavefunction',wavefunc)
call wavefunc%print(iOut,'all')
call matFile%getMolData(moleculeInfo)
call matFile%get2ERIs('regular',eris(1))
call eris(1)%print(iOut,'AO 2ERIs')
nElec = wavefunc%nElectrons%ival()
nAlpha = wavefunc%nAlpha%ival()
nBeta = wavefunc%nBeta%ival()
nBasis = wavefunc%nbasis%ival()
multi = wavefunc%multiplicity%ival()
!
!     Compute the nuclear-nuclear repulsion energy.
!
Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
call moleculeInfo%print(iOut)
call Vnn%print(iOut,'Nuclear Repulsion Energy (au)')
!
!     Initialize the density matrix
call wavefunc%density_matrix%init(nBasis,nBasis)
call wavefunc%density_matrix%print(iOut,'Initial density matrix')
!
!     Determine othogonalization matrix  
call mqc_scf_transformation_matrix(wavefunc%overlap_matrix,Xmat,nBasUse)
call Xmat%print(iOut,'orthogonalization matrix')
if(nBasUse.ne.nBasis) call mqc_error('Linear dependencies not implemented') 
!
!     enter iterations 
do while (iter.le.max_iter)
    write(iOut,'(A,I3)') 'Iteration #: ', iter
!
!     form G matrix 
Gmat = contraction ([eris],wavefunc%density_matrix)
call Gmat%print(iOut,'G(P)')
!
!     form  Fock matrix 
Fock = wavefunc%core_hamiltonian + Gmat
call Fock%print(iOut,'Fock matrix')
!
!     Compute the energy of the current iteration. Save the energy in the Energy variable.
!
!     Orthogonalize Fock basis
Fock = matmul(dagger(Xmat),matmul(Fock,Xmat))
call Fock%print(iOut,'Orthogonalized Fock matrix')
!
!     Diagonalize Fock matrix 
call Fock%diag(wavefunc%MO_energies,wavefunc%MO_coefficients)
!
!     Back-transform  MO coefficients 
call wavefunc%MO_coefficients%print(iOut,' Orthogonal MO coeffiecnts') 
wavefunc%MO_coefficients = matmul(Xmat,wavefunc%MO_coefficients)
if(wavefunc%wf_type.eq.'G') then
  wavefunc%MO_coefficients = wavefunc%MO_coefficients%orbitals(alphaOrbsIn=[(i,i=1,nbasis,2),(i,i=-2,-1*nBasis,-2)], &
    betaOrbsIn=[(i,i=-2,-1*nBasis,-2),(i,i=1,nBasis,2)])
endIf
call wavefunc%MO_energies%print(iOut,'Orbital energies')
call wavefunc%MO_coefficients%print(iOut,'MO coefficients')
!
!     form density matrix 
old_density = wavefunc%density_matrix 
wavefunc%density_matrix = matmul(wavefunc%MO_coefficients%orbitals('occupied',[nAlpha],[nBeta]), &
dagger(wavefunc%MO_coefficients%orbitals('occupied',[nAlpha],[nBeta])))
call wavefunc%density_matrix%print(iOut,'Density matrix')
!
!     compute Energy at iteration
Energy = contraction(wavefunc%core_hamiltonian,wavefunc%density_matrix)
Energy = Energy + half*contraction(wavefunc%density_matrix,Gmat)
Energy = Energy + Vnn
call Energy%print(iOut,'Energy')
!
!     test convergence and exit or iterate
conver = mqc_integral_norm((wavefunc%density_matrix-old_density),'F')
call conver%print(iOut,'Convergence on density matrix')
if(conver.le.thresh) exit
if(iter.eq.max_iter) call mqc_error('Maximum number of iterations reached')
iter = iter+1

endDo

call energy%print(iOut,'Final energy')
!*********RUNTIME PRINTING SECTION********
call cpu_time(finish)
runtime = finish-start
 print*,'runtime in seconds: ', runtime
!*******************
end program hartreefock

