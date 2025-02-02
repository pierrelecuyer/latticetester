**************************************************************
TestBasisConstructSpeedTri with m = 1048573
Types: Int = NTL::ZZ, Real = double
LLL with delta = 0.5
Number of replications (different multipliers a): 1000
Total time (includes time for other operations between triangularizations): 39.8568 seconds.

Dimension:           5        10        20        30        40        50        60        70  

Timings for the different tasks, in basic clock units (microseconds): 
LLLPrimal         11020     30630     75784    126755    194028    272800    386172    504242 
mDualBasis        18406     68956    296381    695048   1375940   2401785   3763748   5822523 
LLLDualmDual       2345      7721     34422     63813    105717    188594    339405    572687 
LowTriP            4155      7792     18358     29108     41242     54977     69851     86566 
mDualLow            893      2362     11975     37416     82247    150454    247439    373449 
LLLDualUT          7382     24955     55120     82514    123192    177130    254118    356616 
UppTriP            3863      7617     18854     29270     41396     54973     69536     85109 
UppTriP2           1972      3501      8000     13992     21843     31244     41680     53267 
mDualUp             831      2302     11681     36955     83995    159291    267457    414301 
LLLDualLT          7304     24373     53809     78673    115812    166293    232287    321837 
LowTriDual         7261     23929     55221     86816    127242    175466    230325    293517 
UppTriDual         6503     23514     94215    201257    331340    491562    671524    873623 
UppTriDualOld     14711     62418    275570    651031   1245900   2112957   3295329   4853358 
UppTriDual2        3865     10981     32576     61241     98267    143325    213744    277215 

Sum of squares after each of the LLL in dual, with delta = 0.5: 
LLLDualmDual     151859     15453      7359      6735      6361      6027      5867      5708 
LLLDualUT        152584     15928      8649      7658      7066      6858      6476      6312 
LLLDualLT        152578     16137      9018      7967      7404      7083      6801      6593 

