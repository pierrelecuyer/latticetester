**************************************************************
TestBasisConstructSpeedTri with m = 1048573
Types: Int = NTL::ZZ, Real = double
LLL with delta = 0.99999
Number of replications (different multipliers a): 1000
Total time (includes time for other operations between triangularizations): 61.0775 seconds.

Dimension:           5        10        20        30        40        50        60        70  

Timings for the different tasks, in basic clock units (microseconds): 
LLLPrimal         12305     43870    217933    427270    659722    913619   1233272   1587766 
mDualBasis        18638     70332    355749    982000   1869483   3268615   5126859   7926519 
LLLDualmDual       2691     11916     73353    248017    626412   1227947   2130378   3479951 
LowTriP            4257      8118     22248     39466     56533     74998     95156    116368 
mDualLow            887      2464     11778     34327     79518    151787    255034    389637 
LLLDualUT          7852     32124    144756    275906    378289    503053    683759    949790 
UppTriP            3900      7815     21151     38830     56343     74482     94811    116582 
UppTriP2           2003      3715      7896     13660     21994     32176     43806     56807 
mDualUp             849      2415     11056     32913     78172    151613    262663    419428 
LLLDualLT          7853     31958    141664    253363    314265    383976    473555    593184 
LowTriDual         7358     26978     97219    159410    203715    256575    318259    387717 
UppTriDual         6585     23962     92416    219642    419265    663189    941664   1242482 
UppTriDualOld     14889     61300    269112    692033   1347812   2306603   3615976   5353747 
UppTriDual2        3902     11215     35333     68702    107825    163304    231034    297249 

Sum of squares after each of the LLL in dual: 
LLLDualmDual     151790     15308      6790      5788      5371      5108      4909      4807 
LLLDualUT        151767     15282      6773      5866      5484      5324      5189      5053 
LLLDualLT        151782     15288      6831      6003      5658      5485      5359      5250 

