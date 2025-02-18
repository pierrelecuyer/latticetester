**************************************************************
TestMatrixCreationSpeed, outside loop is on dimension, inside loop on replications. 
Types: Int = NTL::ZZ, Real = double
Number of repetitions: 100000
Dimension:                            5        10        20        30        40        50      total 

Timings for the different tasks, in basic clock units (microseconds): 
One matrix, no resize               2730      2811      4410      5366      5299      6429     27045 
One matrix, resize rows often        950      1790      3129      4357      5158      6353     21737 
One matrix, resize when needed       790      1522      2973      4128      5383      6052     20848 
One matrix, resize all often       13117     28269     81691    131946    192229    264850    712102 
New matrix for each repet.         13699     35076     84632    139407    208210    282660    763684 

Sums of diagonal entries:
One matrix, no resize                  500005         1000035         2000170         3000405         4000740         5001175 
One matrix, resize rows often          500005         1000035         2000170         3000405         4000740         5001175 
One matrix, resize when needed         500005         1000035         2000170         3000405         4000740         5001175 
One matrix, resize all often           500005         1000035         2000170         3000405         4000740         5001175 
New matrix for each repet.             500005         1000035         2000170         3000405         4000740         5001175 

**************************************************************
TestMatrixCreationSpeed, inside loop is on dimension. 
Types: Int = NTL::ZZ, Real = double
Number of repetitions: 100000
Timings (total) in basic clock units (microseconds): 
One matrix, no resize              19619
One matrix, resize rows often      20483
One matrix, resize when needed    767782
One matrix, resize all often      759866
New matrix for each repet.        835051

