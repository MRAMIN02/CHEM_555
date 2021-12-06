#P rhf nosymm output=matrix int=noraf scf=conven

c8

0 1
C
C 1 1.0
H 1 1.0 2 110.6
H 1 1.0 2 110.6 3 -120.0
H 1 1.0 2 110.6 3  120.0
H 2 1.0 1 110.6 5  180.0
H 2 1.0 1 110.6 5   60.0
H 2 1.0 1 110.6 5  -60.0

c8.mat

