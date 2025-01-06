
=======================================================================
Compare L2 vs L1 norms, and Cholesky vs Triangular decompositions.
*********************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 100
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY

Timings are in basic clock units (microseconds) 
Num dimens:            4          6          8         10         12         14 
Microseconds:        766       1558       2981       5094       8067      12583 
Sums squares:    1948316    6830430   14920230   25362160   34928517   47362376 
Num calls BB:        448        720       1170       2346       4130      10374 
Total time for everything: 0.031781 seconds

*********************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 100
PRIMAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR

Timings are in basic clock units (microseconds) 
Num dimens:            4          6          8         10         12         14 
Microseconds:       1075       2185       4040       6991      10928      17786 
Sums squares:    1948316    6830430   14920230   25362160   34928517   47362376 
Num calls BB:      17427      41938      83160     148227     228770     377332 
Total time for everything: 0.043635 seconds

*********************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 100
PRIMAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY

Timings are in basic clock units (microseconds) 
Num dimens:            4          6          8         10         12         14 
Microseconds:       1033       4979      87268     288208     646189    1081678 
Sums squares:    5442331   25085445   71146801  101150440  102658609  103201855 
Num calls BB:        991       5379      69845     227238     349660     453871 
Total time for everything: 2.110001 seconds

*********************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 100
PRIMAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR

Timings are in basic clock units (microseconds) 
Num dimens:            4          6          8         10         12         14 
Microseconds:       1221       2794       6445       9835      12093      15109 
Sums squares:    5442331   25085445   71146801  101150440  102658609  103201855 
Num calls BB:      28841      85246     218789     315975     321493     323630 
Total time for everything: 0.048135 seconds

*********************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 100
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY

Timings are in basic clock units (microseconds) 
Num dimens:            4          6          8         10         12         14 
Microseconds:        691       1345       2354       3892       5961       9144 
Sums squares:       1874        683        502        423        390        377 
Num calls BB:        442        735       1426       2929       6446      16150 
Total time for everything: 0.025856 seconds

*********************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 100
DUAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR

Timings are in basic clock units (microseconds) 
Num dimens:            4          6          8         10         12         14 
Microseconds:       1385       2809       5332       9338      14851      23780 
Sums squares:       1874        683        502        423        390        377 
Num calls BB:      20924      47343     107344     192439     295186     468503 
Total time for everything: 0.059832 seconds

*********************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 100
DUAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY

Timings are in basic clock units (microseconds) 
Num dimens:            4          6          8         10         12         14 
Microseconds:       1168       4583      32708     326925    1662033   15212026 
Sums squares:       4971       2208       1773       1567       1449       1414 
Num calls BB:        974       4779      37556     305999    1684635   13658944 
Total time for everything: 17.24211 seconds

*********************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 100
DUAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR

Timings are in basic clock units (microseconds) 
Num dimens:            4          6          8         10         12         14 
Microseconds:       1494       3068       6050      10438      16144      26526 
Sums squares:       4971       2208       1773       1567       1449       1414 
Num calls BB:      26252      57946     132118     219199     321037     533459 
Total time for everything: 0.066129 seconds

