****************************************************
Types: Int = NTL::ZZ, Real = double

TestReducersSpeed with m = 1099511627791
Results for `testLoop`
Timings (in microseconds) for different methods for 10 replications 

 dim:               5        10        20        30        40  

LLL5               275       919      2402      3135      4565 
LLL99999           160      1181      5910     14591     20947 
BKZ99999           157      1020      9187     30831     72865 
LLL8+BB            270      1737      6374    540479         0 
LLL99999+BB        195      1290      8268    108144  40961018 
BKZ99999+BB        224      1534     10000     72953   7110022 
L8+BKZ+BB          242      2140      8994     71750   7145572 
L5+99+BKZ+BB       283      1765     10183     70708   6356085 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5            248016       2239        552        421        388 
LLL99999        248016       2197        249        158        138 
BKZ99999        248016       2197        249        150        121 
LLL8+BB         248016       2197        249        147          0 
LLL99999+BB     248016       2197        249        147        116 
BKZ99999+BB     248016       2197        249        147        116 
L8+BKZ+BB       248016       2197        249        147        116 
L5+99+BKZ+BB    248016       2197        249        147        116 

Total time for everything: 62.6894 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = xdouble

TestReducersSpeed with m = 1099511627791
Results for `testLoop`
Timings (in microseconds) for different methods for 10 replications 

 dim:               5        10        20        30        40  

LLL5               450      2071      7558     11158     28172 
LLL99999           299      3507     21845     66826    127558 
BKZ99999           339      3074     38817    155916    400102 
LLL8+BB            444      2969     31024   3811241         0 
LLL99999+BB        385      3521     37568    812071 315710506 
BKZ99999+BB        443      3613     40400    464085  48798144 
L8+BKZ+BB          454      3258     31917    429582  58149551 
L5+99+BKZ+BB       522      4550     38578    478395  53842188 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5            248016       2239        552        421        388 
LLL99999        248016       2197        249        158        138 
BKZ99999        248016       2197        249        150        121 
LLL8+BB         248016       2197        249        147          0 
LLL99999+BB     248016       2197        249        147        116 
BKZ99999+BB     248016       2197        249        147        116 
L8+BKZ+BB       248016       2197        249        147        116 
L5+99+BKZ+BB    248016       2197        249        147        116 

Total time for everything: 483.591 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = quad_float

TestReducersSpeed with m = 1099511627791
Results for `testLoop`
Timings (in microseconds) for different methods for 10 replications 

 dim:               5        10        20        30        40  

LLL5               392      2088      6030     11202     17540 
LLL99999           277      2512     22114     75345     99215 
BKZ99999           368      3369     36737    187809    451304 
LLL8+BB            357      2697     23943   4318803         0 
LLL99999+BB        353      3134     39499    844097 336373475 
BKZ99999+BB        358      3746     43927    518637  61525967 
L8+BKZ+BB          505      4237     35206    434031  57108617 
L5+99+BKZ+BB       554      3471     44445    495754  61317563 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5            248016       2239        552        421        388 
LLL99999        248016       2197        249        158        138 
BKZ99999        248016       2197        249        150        121 
LLL8+BB         248016       2197        249        147          0 
LLL99999+BB     248016       2197        249        147        116 
BKZ99999+BB     248016       2197        249        147        116 
L8+BKZ+BB       248016       2197        249        147        116 
L5+99+BKZ+BB    248016       2197        249        147        116 

Total time for everything: 524.087 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = NTL::RR

TestReducersSpeed with m = 1099511627791
Results for `testLoop`
Timings (in microseconds) for different methods for 10 replications 

 dim:               5        10        20        30        40  

LLL5              2028     11764     55171     77861    142702 
LLL99999          1158     19028    208733    639642    990106 
BKZ99999          1505     21703    281801   1406857   3938710 
LLL8+BB           1530     15812    176528  31690372         0 
LLL99999+BB       1352     19764    221740   6179586 2592938278 
BKZ99999+BB       1647     21493    330931   3778797 453928248 
L8+BKZ+BB         1846     22971    237343   3710493 482443960 
L5+99+BKZ+BB      1897     20758    295789   3637369 366505919 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5            248016       2239        552        421        388 
LLL99999        248016       2197        249        158        138 
BKZ99999        248016       2197        249        153        123 
LLL8+BB         248016       2197        249        147          0 
LLL99999+BB     248016       2197        249        147        116 
BKZ99999+BB     248016       2197        249        147        116 
L8+BKZ+BB       248016       2197        249        147        116 
L5+99+BKZ+BB    248016       2197        249        147        116 

Total time for everything: 3954.01 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


