program workshop_8

use mqc_gaussian
implicit none

character(len=:), allocatable :: fileName
type(mqc_gaussian_unformatted_matrix_file)::matFile
type(mqc_wavefunction)::wavefunction
type(mqc_molecule_data)::molecule
type(mqc_scalar)::Vnn
type(mqc_scf_eigenvalues ):: energy
type(mqc_scf_integral ):: coefficients

call mqc_get_command_argument(1,fileName)
call matFile%load(fileName)
call matFile%getESTObj('wavefunction',wavefunction)
call matFile%getMolData(molecule)
Vnn = molecule%getNucRep()
call Vnn%print(6,'Nuclear Repulsion Energy (au)')
call molecule%print(6)
call wavefunction%print(6,'overlap')
call wavefunction%print(6,'core hamiltonian')
call wavefunction%print(6,'mo coefficients')
call wavefunction%print(6,'orbital energies')

wavefunction%Overlap_Matrix = &
   wavefunction%Overlap_Matrix%orbitals(alphaOrbsIn=[2,6],axis=1)
wavefunction%Overlap_Matrix = &
   wavefunction%Overlap_Matrix%orbitals(alphaOrbsIn=[2,6],axis=2)
call wavefunction%print(6,'overlap')
wavefunction%Core_Hamiltonian = &
   wavefunction%Core_Hamiltonian%orbitals(alphaOrbsIn=[2,6],axis=1)
wavefunction%Core_Hamiltonian = &
   wavefunction%Core_Hamiltonian%orbitals(alphaOrbsIn=[2,6],axis=2)
call wavefunction%print(6,'core hamiltonian')

call wavefunction%Core_Hamiltonian%eigensys(wavefunction%Overlap_Matrix,energy,coefficients)
call energy%print(6,'Huckel orbital energies')
call coefficients %print(6,'Huckel orbital coefficients')

end program workshop_8

