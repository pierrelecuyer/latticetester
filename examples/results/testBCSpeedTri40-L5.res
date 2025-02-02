**************************************************************
TestBasisConstructSpeedTri with m = 1099511627791
Types: Int = NTL::ZZ, Real = double
LLL with delta = 0.5
Number of replications (different multipliers a): 1000
Total time (includes time for other operations between triangularizations): 66.1228 seconds.

Dimension:           5        10        20        30        40        50        60        70  

Timings for the different tasks, in basic clock units (microseconds): 
LLLPrimal         17056     67535    224808    406621    669468   1008176   1441399   1986495 
mDualBasis        22951    107525    606000   1586903   3114663   5357904   8735251  13112882 
LLLDualmDual       2485      7719     41978     91636    141429    215628    371684    639091 
LowTriP            4354      8074     21428     33671     46946     61583     77690     94271 
mDualLow            959      2447     11505     36044     81561    150934    249679    379418 
LLLDualUT         11294     41166    110638    153111    205156    273312    367328    486884 
UppTriP            3942      7640     20836     34103     46840     61152     77320     94629 
UppTriP2           2004      3376      7400     12982     20260     29104     38767     49751 
mDualUp             811      2326     11096     34141     78909    150710    255054    397538 
LLLDualLT          9646     37877    102498    141079    188595    249021    327039    427511 
LowTriDual         7538     27095     81832    117833    160245    211429    270323    336945 
UppTriDual         6947     24677    102941    227294    385100    573943    792551   1034139 
UppTriDualOld     19157     75321    361897    820844   1496258   2443413   3719883   5352611 
UppTriDual2        3937     10875     33301     62133     98969    148263    206553    268886 

Sum of squares after each of the LLL in dual, with delta = 0.5 :
LLLDualmDual  3.72512e+07    227899     31490     27101     27413     26336     25707     24705 
LLLDualUT     3.738e+07    240555     46798     40303     35597     34554     31690     30720 
LLLDualLT     3.73501e+07    240824     46535     39233     35574     33287     31688     30531 

