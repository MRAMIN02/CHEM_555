#p output=matrix nosymm scf=conven int=noraff

h2-sto3g

0 1
H              
H 1 0.6

h2-sto3g.mat



