#P uhf guess=mix stable =opt nosymm output=matrix int=noraf scf=conven

test

0 1
H
H 1 1.1

h2_uhf.mat


