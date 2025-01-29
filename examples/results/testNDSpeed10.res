
=======================================================================
Compare L2 vs L1 norms, and Cholesky vs Triangular decompositions.
Types: Int = NTL::ZZ, Real = double
TestNormDecomp with m = 1021
Number of replications (different multipliers a): 1000
Timings are in basic clock units (microseconds) 
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10         12         14 
Microseconds:        6521      13954      26165      43909      72486     112838 
Aver. squares:  19172.243  68906.201 148415.328 248996.637 343256.604 464481.866 
Aver. calls BB:         4          7         11         22         39         98 
Max |z_j|:              1          1          1          1          2          3 
Total time for everything: 0.282518 seconds

PRIMAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10         12         14 
Microseconds:       10438      21379      38904      65741     104565     167629 
Aver. squares:  19172.243  68906.201 148415.328 248996.637 343256.604 464481.866 
Aver. calls BB:       172        422        828       1454       2242       3671 
Max |z_j|:            180        336        466        595        679        791 
Total time for everything: 0.415039 seconds

PRIMAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10         12         14 
Microseconds:       10082      47923     835391    2753778    5835211   11085862 
Aver. squares:  53421.824  251337.26 703808.764 1006605.287 1024485.797 1030976.305 
Aver. calls BB:         9         54        689       2220       3334       4592 
Max |z_j|:              2          3          4          4          5          6 
Total time for everything: 20.574816 seconds

PRIMAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10         12         14 
Microseconds:       11885      27063      60347      91072     111269     141446 
Aver. squares:  53421.824  251337.26 703808.764 1006605.287 1024485.797 1030976.305 
Aver. calls BB:       284        853       2166       3143       3208       3234 
Max |z_j|:            333        670       1112       1444       1934       2399 
Total time for everything: 0.449534 seconds

DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10         12         14 
Microseconds:        6752      13620      23964      38804      59869      91051 
Aver. squares:     18.683      6.805      4.963      4.162      3.859      3.743 
Aver. calls BB:         4          7         14         27         62        156 
Max |z_j|:              1          1          1          2          2          2 
Total time for everything: 0.263849 seconds

DUAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10         12         14 
Microseconds:       13508      27588      51799      82548     124621     188700 
Aver. squares:     18.683      6.805      4.963      4.162      3.859      3.743 
Aver. calls BB:       208        477       1055       1802       2825       4566 
Max |z_j|:        1301506  169328969 15517834240 1977540128404 4402260669092 425972927270 
Total time for everything: 0.51811 seconds

DUAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10         12         14 
Microseconds:       10715      42605     313899    3081171   15582314  144769887 
Aver. squares:     48.752     21.984     17.722     15.485      14.31     13.988 
Aver. calls BB:         9         47        383       2959      16448     134192 
Max |z_j|:              2          3          3          4          5          6 
Total time for everything: 163.832716 seconds

DUAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10         12         14 
Microseconds:       14164      29257      54999      90402     134496     207411 
Aver. squares:     48.752     21.984     17.722     15.485      14.31     13.988 
Aver. calls BB:       256        575       1297       2169       3179       5255 
Max |z_j|:        1820680  254304343 24828455645 3949326528329 3311406844471 450998140744 
Total time for everything: 0.560846 seconds

