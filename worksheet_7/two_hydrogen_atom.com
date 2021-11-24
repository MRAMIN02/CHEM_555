%chk=two_hydrogen_atoms.chk
#p

test

0 1
H
H 1 1.0
H 1 R1  2 90.0
H 3 1.0 1 90.0 2 0.0

R1 = 1.0

