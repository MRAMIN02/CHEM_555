%chk=two_p_orbitals_1.chk
#P gen output=matrix

test

0 1
C
C 1 R1
C 2 R2 1 180.0
C 3 R3 2 180.0 1 0.0

R1 = 1.0
R2 = 1.1
R3 = 1.2

C 0
SP 1 1.00
0.1478600533D+00 0.3995128261D+00 0.6076837186D+00
****

two_p_orbitals_1.mat

