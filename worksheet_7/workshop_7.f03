program workshop_7

use mqc_gaussian
implicit none

character(len=:),allocatable :: fileName
type(mqc_gaussian_unformatted_matrix_file)::matFile
type(mqc_wavefunction)::wavefunction
type(mqc_molecule_data)::molecule
type(mqc_scalar)::Vnn

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

end program workshop_7
