# 1mm scale
setscale 1
setoption costfunction bbr
setoption optimisationtype brent
setoption tolerance 0.0005000 0.0005000 0.0005000 0.0200000 0.0200000 0.0200000 0.002000 0.002000 0.002000 0.001000 0.001000 0.001000
setoption boundguess 1
clear UU
clear UV
clear U
setrow UU 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UU:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 8
setoption optimisationtype powell
optimise 12 U:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 8
setoption optimisationtype brent
optimise 12 U:2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
sort U

