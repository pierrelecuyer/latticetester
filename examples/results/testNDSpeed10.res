
=======================================================================
Compare L2 vs L1 norms, and Cholesky vs Triangular decompositions.
Types: Int = NTL::ZZ, Real = double
TestNormDecomp with m = 1021
Number of replications (different multipliers a): 1000

*********************************************************
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10         12         14 
Microseconds:        6450      13830      25889      43719      71770     112120 
Aver. squares:  19172.243  68906.201 148415.328 248996.637 343256.604 464481.866 
Aver. calls BB:         4          7         11         22         39         98 
Total time for everything: 0.280053 seconds

*********************************************************
PRIMAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10         12         14 
Microseconds:       10291      20998      38756      65731     104105     167362 
Aver. squares:  19172.243  68906.201 148415.328 248996.637 343256.604 464481.866 
Aver. calls BB:       172        422        828       1454       2242       3671 
Total time for everything: 0.413794 seconds

*********************************************************
PRIMAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10         12         14 
Microseconds:        9932      47286     823962    2716689    5724406   10920108 
Aver. squares:  53421.824  251337.26 703808.764 1006605.287 1024485.797 1030976.305 
Aver. calls BB:         9         54        689       2220       3334       4592 
Total time for everything: 20.248777 seconds

*********************************************************
PRIMAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10         12         14 
Microseconds:       11493      26430      59100      89296     108892     138097 
Aver. squares:  53421.824  251337.26 703808.764 1006605.287 1024485.797 1030976.305 
Aver. calls BB:       284        853       2166       3143       3208       3234 
Total time for everything: 0.439721 seconds

*********************************************************
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10         12         14 
Microseconds:        6734      13448      23816      38654      59358      90120 
Aver. squares:     18.683      6.805      4.963      4.162      3.859      3.743 
Aver. calls BB:         4          7         14         27         62        156 
Total time for everything: 0.26163 seconds

*********************************************************
DUAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10         12         14 
Microseconds:       13232      27406      50007      81111     122046     184529 
Aver. squares:     18.683      6.805      4.963      4.162      3.859      3.743 
Aver. calls BB:       208        477       1055       1802       2825       4566 
Total time for everything: 0.508335 seconds

*********************************************************
DUAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10         12         14 
Microseconds:       10899      42349     312071    3064502   15573982  144380934 
Aver. squares:     48.752     21.984     17.722     15.485      14.31     13.988 
Aver. calls BB:         9         47        383       2959      16448     134192 
Total time for everything: 163.417504 seconds

*********************************************************
DUAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10         12         14 
Microseconds:       14057      29198      54184      89037     131649     202759 
Aver. squares:     48.752     21.984     17.722     15.485      14.31     13.988 
Aver. calls BB:       256        575       1297       2169       3179       5255 
Total time for everything: 0.551558 seconds

