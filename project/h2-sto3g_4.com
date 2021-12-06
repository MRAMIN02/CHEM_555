#p rhf nosymm output=matrix int=noraf scf=conven

h2-sto3g

0 1
H              
H 1 1.2
H 2 1.2 1 180. 
H 3 1.2 2 180. 1 0.0

h2-sto3g_4.mat

