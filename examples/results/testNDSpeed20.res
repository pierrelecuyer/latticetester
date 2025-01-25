
=======================================================================
Compare L2 vs L1 norms, and Cholesky vs Triangular decompositions.
Types: Int = NTL::ZZ, Real = NTL::RR
TestNormDecomp with m = 1048573
Number of replications (different multipliers a): 1000
Timings are in basic clock units (microseconds) 
*********************************************************
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10 
Microseconds:       91392     295727     729583    1518998 
Aver. squares:  574748394 6882739584 2.696390723e+10 6.349727849e+10 
Aver. calls BB:         4          6         10         18 
Total time for everything: 2.640739 seconds

*********************************************************
PRIMAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10 
Microseconds:    16727791   66251719          0          0 
Aver. squares:  574748394 6882739584          0          0 
Aver. calls BB:     24216      93236          0          0 
Total time for everything: 82.982281 seconds

*********************************************************
PRIMAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10 
Microseconds:      108385     426770    2867497   45972817 
Aver. squares: 1601226335 2.681784853e+10 1.357848086e+11 3.924602016e+11 
Aver. calls BB:         8         54        680      10777 
Total time for everything: 49.380843 seconds

*********************************************************
PRIMAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10 
Microseconds:    19020298   92880047          0          0 
Aver. squares: 1601226335 2.681784853e+10          0          0 
Aver. calls BB:     40440     189829          0          0 
Total time for everything: 111.902984 seconds

*********************************************************
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10 
Microseconds:       75651     208867     473994     956314 
Aver. squares:    548.739     66.099     25.916     15.301 
Aver. calls BB:         4          6         11         20 
Total time for everything: 1.734311 seconds

*********************************************************
DUAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10 
Microseconds:    24941961   98359657          0          0 
Aver. squares:    548.739     66.099          0          0 
Aver. calls BB:     29755     100178          0          0 
Total time for everything: 123.309902 seconds

*********************************************************
DUAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10 
Microseconds:       91089     342252    2502629   32882394 
Aver. squares:   1523.097    256.456    126.819      85.13 
Aver. calls BB:         8         56        683       8759 
Total time for everything: 35.838621 seconds

*********************************************************
DUAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10 
Microseconds:    26962131  113116951          0          0 
Aver. squares:   1523.097    257.224          0          0 
Aver. calls BB:     42872     146147          0          0 
Total time for everything: 140.087386 seconds

