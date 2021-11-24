%chk=two_p_orbitals.chk
#P gen output=matrix

test

0 1
C
C 1 R1

R1 = 5.0

C 0
SP 1 1.00
0.1478600533D+00 0.3995128261D+00 0.6076837186D+00
****

two_p_orbitals.mat

