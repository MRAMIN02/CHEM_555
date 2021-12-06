#p rhf nosymm output=matrix int=noraf scf=conven

h2-sto3g

0 1
H              
H 1 1.2
H 2 1.2 1 180.
H 3 1.2 2 180. 1 0.0
H 4 1.2 3 180. 1 0.0
H 5 1.2 4 180. 1 0.0
H 6 1.2 5 180. 1 0.0
H 7 1.2 6 180. 1 0.0

h2-8.mat

