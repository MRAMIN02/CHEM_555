program mulliken 

use MQC_Gaussian
use MQC_Algebra
use MQC_EST

implicit none
character(len=:),allocatable :: fileName
type(mqc_gaussian_unformatted_matrix_file)::fileInfo
type(mqc_molecule_data)::moldata
type(mqc_scf_integral)::density,overlap
type(mqc_matrix)::tmpMat,sMat
type(mqc_vector)::MullikenPop
integer::iout=6,nBasis,nAtoms,i,j
integer,dimension(:),allocatable::basisArray

call mqc_get_command_argument(1,fileName)
call fileInfo%load(fileName)
nBasis = fileInfo%getVal('nBasis')
call fileInfo%getESTObj('overlap',est_integral=overlap)
call fileInfo%getESTObj('density',est_integral=density)
call fileInfo%getMolData(moldata)
nAtoms = moldata%getNumAtoms()
basisArray = fileInfo%getBasisArray('basis2Atom')
sMat = overlap%getBlock('alpha') 
tmpMat = density%getblock('alpha') + density%getblock('beta')
call mullikenPop%init(nAtoms)
do i = 1, nBasis
  do j = 1, nBasis
    call mullikenPop%put(mullikenPop%at(basisArray(i))+0.5*sMat%at(i,j)*tmpMat%at(j,i),basisArray(i))
    call mullikenPop%put(mullikenPop%at(basisArray(j))+0.5*sMat%at(i,j)*tmpMat%at(j,i),basisArray(j))
  endDo
endDo
call mullikenPop%print(iout,'Mulliken Electron Population')
mullikenPop = moldata%nuclear_charges - mullikenPop
call mullikenPop%print(iout,'Mulliken Charges')

end program mulliken
