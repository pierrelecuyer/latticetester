**************************************************************
TestBasisConstructSpeedTri with m = 1099511627791
Types: Int = NTL::ZZ, Real = double
LLL with delta = 0.99999
Number of replications (different multipliers a): 1000
Total time (includes time for other operations between triangularizations): 131.911 seconds.

Dimension:           5        10        20        30        40        50        60        70  

Timings for the different tasks, in basic clock units (microseconds): 
LLLPrimal         19004     95821    681432   1892975   3302096   4886668   6738622   8908062 
mDualBasis        24169    114329    644393   2220320   5113427   8818025  14296552  21438928 
LLLDualmDual       2799     12317     81433    278134    734704   1539194   2699556   4246538 
LowTriP            4439      8391     22427     43701     68920     90329    113085    138011 
mDualLow            961      2532     11610     33480     75148    145814    248815    387313 
LLLDualUT         11993     51315    268058    678178   1023630   1234750   1479642   1794158 
UppTriP            4014      7964     21054     41383     66794     88915    112627    136650 
UppTriP2           2035      3682      7399     12478     19871     28959     40023     52215 
mDualUp             827      2444     11066     31743     71534    140259    244449    383428 
LLLDualLT         10316     48488    256121    643701    934655   1075884   1237928   1428963 
LowTriDual         7730     28949    126672    287769    406295    470286    540595    619446 
UppTriDual         7140     25814    105717    261390    513903    844453   1224525   1642319 
UppTriDualOld     19282     75931    371689   1018218   2035255   3309216   4952798   7027240 
UppTriDual2        3928     11023     34267     70031    116845    167176    229916    298828 

Sum of squares after each of the LLL in dual, with delta = 0.99999 :
LLLDualmDual  3.71393e+07    225731     25183     14521     12620     12885     13605     14358 
LLLDualUT     3.71393e+07    225646     25035     14714     12992     12555     12299     12079 
LLLDualLT     3.71393e+07    225642     25152     14923     13271     12819     12487     12266 

