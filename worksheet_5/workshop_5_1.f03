program workshop_5

implicit none

integer::i,j,max_i=5,max_j=5
do i=1,max_i
  do j=max_j,1,-1
  write(*,*) 'i,j equal', i,',',j
  endDo
endDo

end program workshop_5
