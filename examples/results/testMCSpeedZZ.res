**************************************************************
TestMatrixCreationSpeed, outside loop is on dimension, inside loop on replications. 
Types: Int = NTL::ZZ, Real = double
Number of replications: 1000000
Sums of diagonal entries (the sums are not the same for all 5 methods):
Dimension:                                5            10            20            30            40            50  
One matrix, no resize                  2500007500000         7500047500000         17500227500000         32500647500000         52501407500000         77502607500000 
One matrix, resize only rows           2500007500000         7500047500000         17500227500000         32500647500000         52501407500000         77502607500000 
One matrix, resize when needed         2500007500000         5000040000000         10000180000000         15000420000000         20000760000000         25001200000000 
One matrix, resize all often           2500007500000         5000040000000         10000180000000         15000420000000         20000760000000         25001200000000 
New matrix for each repet.             2500007500000         5000040000000         10000180000000         15000420000000         20000760000000         25001200000000 

Timings (total) with the different methods, in basic clock units (microseconds): 
One matrix, no resize             284391 
One matrix, resize only rows      292091 
One matrix, resize when needed    284623 
One matrix, resize all often     7908641 
New matrix for each repet.       8556413 

**************************************************************
TestMatrixCreationSpeed, inside loop is on dimension. 
Types: Int = NTL::ZZ, Real = double
Number of replications: 1000000
Timings (total) in basic clock units (microseconds): 
One matrix, no resize             281340
One matrix, resize only rows      294787
One matrix, resize when needed   7638361
One matrix, resize all often     7607095
New matrix for each repet.       8495962

