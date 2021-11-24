program workshop_5

implicit none

character(len=10)::fruit
fruit='orange'
if (fruit.eq.'pear') then
  write(*,*) 'fruit is peer'
elseIf (fruit.eq.'apple') then
  write(*,*) 'fruit is apple'
elseIf (fruit.eq.'plum') then
  write(*,*) 'fruit is plum'
elseIf (fruit.eq.'orange') then
  write(*,*) 'fruit is orange' 
else
  write(*,*) 'we have a mystery fruit'
endIf

end program workshop_5
