**************************************************************
TestMatrixCreationSpeed, outside loop is on dimension, inside loop on replications. 
Types: Int = long, Real = double
Number of repetitions: 100000
Dimension:                            5        10        20        30        40        50      total 

Timings for the different tasks, in basic clock units (microseconds): 
One matrix, no resize                503       822      1540      2189      2631      2534     10219 
One matrix, resize rows often        386       462       712       930      1148      1360      4998 
One matrix, resize when needed       247       552       963      1129      1508      1568      5967 
One matrix, resize all often        7329     17072     36514     54565     73899     89775    279154 
New matrix for each repet.          6432     16344     36227     52236     69288     84726    265253 

Sums of diagonal entries:
One matrix, no resize             500005   1000035   2000170   3000405   4000740   5001175 
One matrix, resize rows often     500005   1000035   2000170   3000405   4000740   5001175 
One matrix, resize when needed    500005   1000035   2000170   3000405   4000740   5001175 
One matrix, resize all often      500005   1000035   2000170   3000405   4000740   5001175 
New matrix for each repet.        500005   1000035   2000170   3000405   4000740   5001175 

**************************************************************
TestMatrixCreationSpeed, inside loop is on dimension. 
Types: Int = long, Real = double
Number of repetitions: 100000
Timings (total) in basic clock units (microseconds): 
One matrix, no resize               3495
One matrix, resize rows often       4306
One matrix, resize when needed    279681
One matrix, resize all often      286974
New matrix for each repet.        261856

