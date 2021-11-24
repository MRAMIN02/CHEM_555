program workshop_5

implicit none

character(len=256)::fruit
integer::nargs

nargs=command_argument_count()
if(nargs.ne.1) then
  write(*,*) 'A single fruit must be specified with the command'
  stop
endIf
call get_command_argument(1,fruit)

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
