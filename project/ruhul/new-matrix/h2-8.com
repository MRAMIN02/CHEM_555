#p rhf nosymm output=matrix int=noraf scf=conven

h2-sto3g

0 1
H              
H 1 0.6
H 2 0.6 1 180.
H 3 0.6 2 180. 1 0.0
H 4 0.6 3 180. 1 0.0
H 5 0.6 4 180. 1 0.0
H 6 0.6 5 180. 1 0.0
H 7 0.6 6 180. 1 0.0

h2-8.mat

