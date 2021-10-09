program workshop_3

use mqc_gaussian

implicit none
type(mqc_matrix)::A,B,C

A=reshape([2,1,1,0],[2,2])
B=reshape([1,2,1,3],[2,2])
C=matmul(A,B)
call A%print(6,'Matrix A')
call B%print(6,'Matrix B')
call C%print(6,'Result of A.B')

end program workshop_3
