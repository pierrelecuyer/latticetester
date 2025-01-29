
=======================================================================
Compare L2 vs L1 norms, and Cholesky vs Triangular decompositions.
Types: Int = NTL::ZZ, Real = NTL::RR
TestNormDecomp with m = 1048573
Number of replications (different multipliers a): 1000
Timings are in basic clock units (microseconds) 
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10 
Microseconds:       92378     299087     737623    1535497 
Aver. squares:  574748394 6882739584 2.696390723e+10 6.349727849e+10 
Aver. calls BB:         4          6         10         18 
Max |z_j|:              1          1          1          1 
Total time for everything: 2.669668 seconds

PRIMAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10 
Microseconds:    16880832   66992177          0          0 
Aver. squares:  574748394 6882739584          0          0 
Aver. calls BB:     24216      93236          0          0 
Max |z_j|:          35168     116488          0          0 
Total time for everything: 83.875585 seconds

PRIMAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10 
Microseconds:      108478     427281    2878420   46127812 
Aver. squares: 1601226335 2.681784853e+10 1.357848086e+11 3.924602016e+11 
Aver. calls BB:         8         54        680      10777 
Max |z_j|:              2          3          4          6 
Total time for everything: 49.547191 seconds

PRIMAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10 
Microseconds:    19086380   93192968          0          0 
Aver. squares: 1601226335 2.681784853e+10          0          0 
Aver. calls BB:     40440     189829          0          0 
Max |z_j|:          62198     235466          0          0 
Total time for everything: 112.281916 seconds

DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10 
Microseconds:       76157     210182     477079     961732 
Aver. squares:    548.739     66.099     25.916     15.301 
Aver. calls BB:         4          6         11         20 
Max |z_j|:              1          1          1          1 
Total time for everything: 1.744493 seconds

DUAL lattice,  Norm: L2NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10 
Microseconds:    25595257  100931285          0          0 
Aver. squares:    548.739     66.099          0          0 
Aver. calls BB:     29755     100178          0          0 
Max |z_j|:     7543524485197 9223342208803491621          0          0 
Total time for everything: 126.53489 seconds

DUAL lattice,  Norm: L1NORM,  Decomposition: CHOLESKY
Num dimens:             4          6          8         10 
Microseconds:       91543     344986    2527931   33268288 
Aver. squares:   1523.097    256.456    126.819      85.13 
Aver. calls BB:         8         56        683       8759 
Max |z_j|:              2          3          5          6 
Total time for everything: 36.252888 seconds

DUAL lattice,  Norm: L1NORM,  Decomposition: TRIANGULAR
Num dimens:             4          6          8         10 
Microseconds:    27417609  115429300          0          0 
Aver. squares:   1523.097    257.224          0          0 
Aver. calls BB:     42872     146147          0          0 
Max |z_j|:     13526315159501 9223342400490776417          0          0 
Total time for everything: 142.855114 seconds

