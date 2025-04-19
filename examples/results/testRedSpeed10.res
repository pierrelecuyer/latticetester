
========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 1e-06.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                  358       924      2191      3570      5404 
LLL99999              280      1124      3964      6510     10045 
BKZ99999-10           299      1194      8486     15703     24361 
L5+9+BKZ-10           396      1624      9986     18190     27820 
LLL5+BB               385      1600      9148     18952     25980 
LLL8+BB               367      1832      9229     18536     27141 
LLL99999+BB           333      1821     10000     20946     30398 
BKZ99999-6+BB         338      1828     12430     26980     41628 
BKZ99999-8+BB         313      1733     13081     27231     43163 
BKZ99999-10+BB        303      1669     13638     29914     45500 
BKZ99999-12+BB        308      1622     14490     31981     46159 
BKZ999-6+BB           309      1596     12181     26913     41001 
BKZ999-8+BB           301      1576     13133     27303     42464 
BKZ999-10+BB          303      1579     13472     30328     45615 
BKZ999-12+BB          303      1561     14500     32417     46115 
L8+BKZ-10+BB          394      2102     14339     31461     48179 
L5+9+BKZ-10+BB        452      2322     15296     32498     48779 
L5+9+BKZ-12+BB        405      2131     15489     34915     49821 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5             43693.62  264162.24   918278.3 1014955.68 1019693.12 
LLL99999         43087.06  249732.08  865143.84 1014955.68 1019693.12 
BKZ99999-10      43087.06  249732.08  864198.86 1014955.68 1019693.12 
L5+9+BKZ-10      43087.06  249732.08  861468.34 1014955.68 1019693.12 
LLL5+BB          43087.06  249732.08  860030.46 1014955.68 1019693.12 
LLL8+BB          43087.06  249732.08  860030.46 1014955.68 1019693.12 
LLL99999+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-6+BB    43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-8+BB    43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-10+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-12+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-6+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-8+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-10+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-12+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
L8+BKZ-10+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
L5+9+BKZ-10+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
L5+9+BKZ-12+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 6         28       2872       6494       6883 
LLL8+BB                 6         23       1499       3303       3387 
LLL99999+BB             6         22       1215       3553       3431 
BKZ99999-6+BB           6         22       1158       3430       3872 
BKZ99999-8+BB           6         22       1118       3420       3823 
BKZ99999-10+BB          6         22       1132       3325       3735 
BKZ99999-12+BB          6         22       1060       3217       3510 
BKZ999-6+BB             6         22       1160       3504       3979 
BKZ999-8+BB             6         22       1174       3370       4035 
BKZ999-10+BB            6         22       1133       3306       3778 
BKZ999-12+BB            6         22       1060       3193       3451 
L8+BKZ-10+BB            6         22       1108       3488       3538 
L5+9+BKZ-10+BB          6         22       1114       3268       3657 
L5+9+BKZ-12+BB          6         22       1076       3316       3463 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          1          2          4          4 
LLL8+BB                 1          1          2          6          5 
LLL99999+BB             1          1          1          6          5 
BKZ99999-6+BB           1          1          1          6          6 
BKZ99999-8+BB           1          1          1          6          7 
BKZ99999-10+BB          1          1          1          6          6 
BKZ99999-12+BB          1          1          1          6          6 
BKZ999-6+BB             1          1          1          6          6 
BKZ999-8+BB             1          1          1          6          6 
BKZ999-10+BB            1          1          1          6          6 
BKZ999-12+BB            1          1          1          7          5 
L8+BKZ-10+BB            1          1          1          6          6 
L5+9+BKZ-10+BB          1          1          1          6          6 
L5+9+BKZ-12+BB          1          1          1          6          6 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          2          5          5          6 
LLL8+BB                 1          1          4          5          4 
LLL99999+BB             1          1          6          9          9 
BKZ99999-6+BB           1          1          5          9          6 
BKZ99999-8+BB           1          1          5          9          7 
BKZ99999-10+BB          1          1          6          5          7 
BKZ99999-12+BB          1          1          6          6         10 
BKZ999-6+BB             1          1          5          9         10 
BKZ999-8+BB             1          1          5          9         10 
BKZ999-10+BB            1          1          5          6         10 
BKZ999-12+BB            1          1          5          6         10 
L8+BKZ-10+BB            1          1          8          6         10 
L5+9+BKZ-10+BB          1          1          5          6         10 
L5+9+BKZ-12+BB          1          1          8          6         10 

Total time for everything: 1.341328 seconds

****************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 1e-06.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                  342       868      1739      2891      4477 
LLL99999              281      1083      3546      5306      7357 
BKZ99999-10           303      1201      8298     20919     35211 
L5+9+BKZ-10           429      1646      9239     23526     40109 
LLL5+BB               375      1416      7700     39667    120540 
LLL8+BB               347      1525      8021     41855    130572 
LLL99999+BB           315      1534      8863     41717    139524 
BKZ99999-6+BB         338      1625     11392     44347    142547 
BKZ99999-8+BB         301      1451     12183     43554    128655 
BKZ99999-10+BB        284      1363     12915     47834    100734 
BKZ99999-12+BB        279      1296     13873     48083    106137 
BKZ999-6+BB           280      1258     11269     45348    132622 
BKZ999-8+BB           272      1242     12064     43596    124834 
BKZ999-10+BB          284      1245     12712     47313     99257 
BKZ999-12+BB          275      1226     13707     48359    111548 
L8+BKZ-10+BB          388      1801     13740     46332    102400 
L5+9+BKZ-10+BB        455      2127     14476     49436    115280 
L5+9+BKZ-12+BB        400      1936     14823     52923     99823 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5                 10.8       4.28       3.36       3.16       3.04 
LLL99999            10.72       4.14       3.16          3        2.9 
BKZ99999-10         10.72       4.14       3.16        2.9        2.8 
L5+9+BKZ-10         10.72       4.14       3.16        2.9       2.84 
LLL5+BB             10.72       4.14       3.16        2.9       2.78 
LLL8+BB             10.72       4.14       3.16        2.9       2.78 
LLL99999+BB         10.72       4.14       3.16        2.9       2.78 
BKZ99999-6+BB       10.72       4.14       3.16        2.9       2.78 
BKZ99999-8+BB       10.72       4.14       3.16        2.9       2.78 
BKZ99999-10+BB      10.72       4.14       3.16        2.9       2.78 
BKZ99999-12+BB      10.72       4.14       3.16        2.9       2.78 
BKZ999-6+BB         10.72       4.14       3.16        2.9       2.78 
BKZ999-8+BB         10.72       4.14       3.16        2.9       2.78 
BKZ999-10+BB        10.72       4.14       3.16        2.9       2.78 
BKZ999-12+BB        10.72       4.14       3.16        2.9       2.78 
L8+BKZ-10+BB        10.72       4.14       3.16        2.9       2.78 
L5+9+BKZ-10+BB      10.72       4.14       3.16        2.9       2.78 
L5+9+BKZ-12+BB      10.72       4.14       3.16        2.9       2.78 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 6         33       1913      21736      81753 
LLL8+BB                 6         28       1509      21913      86871 
LLL99999+BB             6         27       1377      20356      90743 
BKZ99999-6+BB           6         26       1132      14107      76107 
BKZ99999-8+BB           6         26       1135      12706      63275 
BKZ99999-10+BB          6         26       1080      11941      38602 
BKZ99999-12+BB          6         26       1113      11491      38327 
BKZ999-6+BB             6         26       1128      14611      68345 
BKZ999-8+BB             6         26       1138      12662      60156 
BKZ999-10+BB            6         26       1080      12419      39017 
BKZ999-12+BB            6         26       1105      11120      40965 
L8+BKZ-10+BB            6         26       1069      10970      38165 
L5+9+BKZ-10+BB          6         26       1092      11458      44392 
L5+9+BKZ-12+BB          6         26       1076      12260      27377 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          1          7         19         20 
LLL8+BB                 1          1          7         20         21 
LLL99999+BB             1          1          6         18         26 
BKZ99999-6+BB           1          1          6         20         43 
BKZ99999-8+BB           1          1          6         17         31 
BKZ99999-10+BB          1          1          6         21         20 
BKZ99999-12+BB          1          1          6         18         18 
BKZ999-6+BB             1          1          6         21         41 
BKZ999-8+BB             1          1          6         17         32 
BKZ999-10+BB            1          1          6         22         18 
BKZ999-12+BB            1          1          6         20         23 
L8+BKZ-10+BB            1          1          6         15         16 
L5+9+BKZ-10+BB          1          1          6         16         22 
L5+9+BKZ-12+BB          1          1          6         21         13 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          2          5          7          7 
LLL8+BB                 1          1          5          6          7 
LLL99999+BB             1          1          5          6          7 
BKZ99999-6+BB           1          1          4          6          6 
BKZ99999-8+BB           1          1          4          5          5 
BKZ99999-10+BB          1          1          4          7          7 
BKZ99999-12+BB          1          1          5          5          7 
BKZ999-6+BB             1          1          4          6          6 
BKZ999-8+BB             1          1          4          5          5 
BKZ999-10+BB            1          1          4          7          8 
BKZ999-12+BB            1          1          4          5          7 
L8+BKZ-10+BB            1          1          4          5          6 
L5+9+BKZ-10+BB          1          1          5          5          9 
L5+9+BKZ-12+BB          1          1          5          6          6 

Total time for everything: 2.73146 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 0.1e-5.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                  708      2715      9399     20049     35431 
LLL99999              768      4187     21515     45222     77646 
BKZ99999-10           875      5577     55365    127167    213403 
L5+9+BKZ-10          1160      6504     57416    131903    224215 
LLL5+BB               871      3825     41111     95022    120889 
LLL8+BB               906      4703     34242     76899    108646 
LLL99999+BB           908      5254     38629     93641    135430 
BKZ99999-6+BB        1012      6194     57647    141883    233405 
BKZ99999-8+BB         996      6355     64831    145977    246126 
BKZ99999-10+BB        986      6427     71411    174483    273595 
BKZ99999-12+BB        990      6420     80402    193725    286101 
BKZ999-6+BB           979      6020     57578    141917    226202 
BKZ999-8+BB           974      6286     65608    145330    239352 
BKZ999-10+BB          976      6346     70610    175733    273958 
BKZ999-12+BB          973      6296     80585    196744    287845 
L8+BKZ-10+BB         1155      6926     71312    174729    285145 
L5+9+BKZ-10+BB       1300      7428     73618    177940    283534 
L5+9+BKZ-12+BB       1273      7394     80468    201386    298452 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5             43693.62  264162.24   918278.3 1014955.68 1019693.12 
LLL99999         43087.06  249732.08  865143.84 1014955.68 1019693.12 
BKZ99999-10      43087.06  249732.08  864198.86 1014955.68 1019693.12 
L5+9+BKZ-10      43087.06  249732.08  861002.98 1014955.68 1019693.12 
LLL5+BB          43087.06  249732.08  860030.46 1014955.68 1019693.12 
LLL8+BB          43087.06  249732.08  860030.46 1014955.68 1019693.12 
LLL99999+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-6+BB    43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-8+BB    43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-10+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-12+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-6+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-8+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-10+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-12+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
L8+BKZ-10+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
L5+9+BKZ-10+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
L5+9+BKZ-12+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 6         28       2872       6532       6858 
LLL8+BB                 6         23       1500       3340       3433 
LLL99999+BB             6         22       1216       3562       3494 
BKZ99999-6+BB           6         22       1154       3410       3900 
BKZ99999-8+BB           6         22       1116       3465       3907 
BKZ99999-10+BB          6         22       1132       3392       3751 
BKZ99999-12+BB          6         22       1060       3196       3564 
BKZ999-6+BB             6         22       1158       3493       3989 
BKZ999-8+BB             6         22       1172       3410       4026 
BKZ999-10+BB            6         22       1133       3362       3781 
BKZ999-12+BB            6         22       1062       3188       3506 
L8+BKZ-10+BB            6         22       1115       3515       3615 
L5+9+BKZ-10+BB          6         22       1099       3251       3681 
L5+9+BKZ-12+BB          6         22       1094       3307       3486 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          1          2          6          5 
LLL8+BB                 1          1          2          8          7 
LLL99999+BB             1          1          1          7          8 
BKZ99999-6+BB           1          1          1          8          7 
BKZ99999-8+BB           1          1          1          7          8 
BKZ99999-10+BB          1          1          1          7          7 
BKZ99999-12+BB          1          1          1          7          8 
BKZ999-6+BB             1          1          1          8          6 
BKZ999-8+BB             1          1          1          8          7 
BKZ999-10+BB            1          1          1          7          7 
BKZ999-12+BB            1          1          1          7          8 
L8+BKZ-10+BB            1          1          1          7          7 
L5+9+BKZ-10+BB          1          1          1          8          7 
L5+9+BKZ-12+BB          1          1          1          7          8 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          2          5          5          5 
LLL8+BB                 1          1          4          5          4 
LLL99999+BB             1          1          6          9          9 
BKZ99999-6+BB           1          1          5          9          6 
BKZ99999-8+BB           1          1          5          9          7 
BKZ99999-10+BB          1          1          6          5          7 
BKZ99999-12+BB          1          1          6          6         10 
BKZ999-6+BB             1          1          5          9         10 
BKZ999-8+BB             1          1          5          9         10 
BKZ999-10+BB            1          1          5          6         10 
BKZ999-12+BB            1          1          5          6         10 
L8+BKZ-10+BB            1          1          8          6         10 
L5+9+BKZ-10+BB          1          1          5          6         10 
L5+9+BKZ-12+BB          1          1          5          6         10 

Total time for everything: 7.480468 seconds

****************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 0.1e-5.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                  716      2375      6339     13510     25458 
LLL99999              723      3656     16253     27388     41340 
BKZ99999-10           825      5287     54211    161426    289794 
L5+9+BKZ-10          1177      6085     55494    171385    303705 
LLL5+BB               861      3392     29940    263816   1008057 
LLL8+BB               866      3990     28936    272671   1078788 
LLL99999+BB           852      4515     33977    256352   1150470 
BKZ99999-6+BB         952      5742     52224    275950   1119301 
BKZ99999-8+BB         924      5889     61555    291898    916244 
BKZ99999-10+BB        914      5935     68928    312502    770958 
BKZ99999-12+BB        909      5861     77144    324079    880447 
BKZ999-6+BB           910      5456     52051    286577   1030219 
BKZ999-8+BB           907      5677     61605    291081    879205 
BKZ999-10+BB          903      5789     69125    311167    734541 
BKZ999-12+BB          900      5741     77161    320325    876970 
L8+BKZ-10+BB         1100      6305     66446    309439    757027 
L5+9+BKZ-10+BB       1302      7024     71677    334174    873398 
L5+9+BKZ-12+BB       1250      6822     77782    339563    772148 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5                 10.8       4.28       3.36       3.16       3.04 
LLL99999            10.72       4.14       3.16          3        2.9 
BKZ99999-10         10.72       4.14       3.16        2.9       2.82 
L5+9+BKZ-10         10.72       4.14       3.16        2.9       2.82 
LLL5+BB             10.72       4.14       3.16        2.9       2.78 
LLL8+BB             10.72       4.14       3.16        2.9       2.78 
LLL99999+BB         10.72       4.14       3.16        2.9       2.78 
BKZ99999-6+BB       10.72       4.14       3.16        2.9       2.78 
BKZ99999-8+BB       10.72       4.14       3.16        2.9       2.78 
BKZ99999-10+BB      10.72       4.14       3.16        2.9       2.78 
BKZ99999-12+BB      10.72       4.14       3.16        2.9       2.78 
BKZ999-6+BB         10.72       4.14       3.16        2.9       2.78 
BKZ999-8+BB         10.72       4.14       3.16        2.9       2.78 
BKZ999-10+BB        10.72       4.14       3.16        2.9       2.78 
BKZ999-12+BB        10.72       4.14       3.16        2.9       2.78 
L8+BKZ-10+BB        10.72       4.14       3.16        2.9       2.78 
L5+9+BKZ-10+BB      10.72       4.14       3.16        2.9       2.78 
L5+9+BKZ-12+BB      10.72       4.14       3.16        2.9       2.78 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 6         33       1927      21863      82129 
LLL8+BB                 6         28       1513      22121      86648 
LLL99999+BB             6         27       1372      19605      90565 
BKZ99999-6+BB           6         27       1149      14192      76266 
BKZ99999-8+BB           6         27       1097      13446      54793 
BKZ99999-10+BB          6         27       1066      12168      37458 
BKZ99999-12+BB          6         27       1098      10921      40058 
BKZ999-6+BB             6         27       1143      14718      68130 
BKZ999-8+BB             6         27       1109      13406      51242 
BKZ999-10+BB            6         27       1067      12892      36115 
BKZ999-12+BB            6         27       1090      10805      40202 
L8+BKZ-10+BB            6         27       1096      11612      36297 
L5+9+BKZ-10+BB          6         27       1116      13287      44369 
L5+9+BKZ-12+BB          6         27       1130      11237      29516 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          2          9         22         25 
LLL8+BB                 1          1          8         23         26 
LLL99999+BB             1          2          8         19         26 
BKZ99999-6+BB           1          2         10         27         51 
BKZ99999-8+BB           1          2         10         28         39 
BKZ99999-10+BB          1          2         11         29         35 
BKZ99999-12+BB          1          2         10         33         48 
BKZ999-6+BB             1          2          9         29         49 
BKZ999-8+BB             1          2         10         28         44 
BKZ999-10+BB            1          2         11         28         35 
BKZ999-12+BB            1          2         11         33         50 
L8+BKZ-10+BB            1          2         10         29         31 
L5+9+BKZ-10+BB          1          2         10         27         39 
L5+9+BKZ-12+BB          1          2         10         28         32 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          2          5          7          7 
LLL8+BB                 1          1          5          6          6 
LLL99999+BB             1          1          5          6          7 
BKZ99999-6+BB           1          1          4          6          6 
BKZ99999-8+BB           1          1          4          5          5 
BKZ99999-10+BB          1          1          4          7          6 
BKZ99999-12+BB          1          1          5          5          7 
BKZ999-6+BB             1          1          4          6          7 
BKZ999-8+BB             1          1          4          5          5 
BKZ999-10+BB            1          1          4          7          8 
BKZ999-12+BB            1          1          4          5          7 
L8+BKZ-10+BB            1          1          4          5          7 
L5+9+BKZ-10+BB          1          1          4          5          8 
L5+9+BKZ-12+BB          1          1          4          6          5 

Total time for everything: 19.222582 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = NTL::RR
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 0.1e-5.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 3845     17670     54957    103276    164476 
LLL99999             4236     30560    158191    271519    413588 
BKZ99999-10          4391     37134    394100    817929   1234233 
L5+9+BKZ-10          5067     38361    400021    848494   1291744 
LLL5+BB              3965     21377    240865    521513    617350 
LLL8+BB              4383     29106    203679    411922    533379 
LLL99999+BB          4403     34375    247173    521028    675441 
BKZ99999-6+BB        4672     39340    399696    862617   1302335 
BKZ99999-8+BB        4547     40360    444768    878230   1359127 
BKZ99999-10+BB       4430     40732    478382   1059629   1512785 
BKZ99999-12+BB       4331     40680    515220   1148924   1547817 
BKZ999-6+BB          4349     39481    400342    855872   1259371 
BKZ999-8+BB          4273     40209    453501    876523   1319018 
BKZ999-10+BB         4272     40534    473680   1064851   1518441 
BKZ999-12+BB         4226     40561    510423   1166492   1567034 
L8+BKZ-10+BB         5060     40870    461524   1063644   1594702 
L5+9+BKZ-10+BB       5442     41828    481762   1082739   1562665 
L5+9+BKZ-12+BB       5182     41742    497552   1200698   1620241 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5             43693.62  264162.24   918278.3 1014955.68 1019693.12 
LLL99999         43087.06  249732.08  865143.84 1014955.68 1019693.12 
BKZ99999-10      43087.06  249732.08  864198.86 1014955.68 1019693.12 
L5+9+BKZ-10      43087.06  249732.08   861651.5 1014955.68 1019693.12 
LLL5+BB          43087.06  249732.08  860030.46 1014955.68 1019693.12 
LLL8+BB          43087.06  249732.08  860030.46 1014955.68 1019693.12 
LLL99999+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-6+BB    43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-8+BB    43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-10+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ99999-12+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-6+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-8+BB      43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-10+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
BKZ999-12+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
L8+BKZ-10+BB     43087.06  249732.08  860030.46 1014955.68 1019693.12 
L5+9+BKZ-10+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 
L5+9+BKZ-12+BB   43087.06  249732.08  860030.46 1014955.68 1019693.12 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 6         28       2871       6500       6885 
LLL8+BB                 6         23       1499       3310       3412 
LLL99999+BB             6         22       1214       3549       3457 
BKZ99999-6+BB           6         22       1159       3386       3875 
BKZ99999-8+BB           6         22       1113       3421       3852 
BKZ99999-10+BB          6         22       1150       3331       3734 
BKZ99999-12+BB          6         22       1062       3197       3506 
BKZ999-6+BB             6         22       1158       3395       3980 
BKZ999-8+BB             6         22       1159       3360       4014 
BKZ999-10+BB            6         22       1151       3309       3776 
BKZ999-12+BB            6         22       1062       3180       3448 
L8+BKZ-10+BB            6         22       1111       3468       3539 
L5+9+BKZ-10+BB          6         22       1117       3276       3673 
L5+9+BKZ-12+BB          6         22       1087       3303       3459 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          1          2          4          4 
LLL8+BB                 1          1          2          6          6 
LLL99999+BB             1          1          1          6          6 
BKZ99999-6+BB           1          1          1          6          6 
BKZ99999-8+BB           1          1          1          6          6 
BKZ99999-10+BB          1          1          1          5          6 
BKZ99999-12+BB          1          1          1          6          5 
BKZ999-6+BB             1          1          1          6          6 
BKZ999-8+BB             1          1          1          6          6 
BKZ999-10+BB            1          1          1          6          6 
BKZ999-12+BB            1          1          1          6          5 
L8+BKZ-10+BB            1          1          1          6          6 
L5+9+BKZ-10+BB          1          1          1          6          7 
L5+9+BKZ-12+BB          1          1          1          6          6 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          2          5          5          5 
LLL8+BB                 1          1          4          5          4 
LLL99999+BB             1          1          6          9          9 
BKZ99999-6+BB           1          1          5          9          6 
BKZ99999-8+BB           1          1          5          9          7 
BKZ99999-10+BB          1          1          6          5          7 
BKZ99999-12+BB          1          1          6          6         10 
BKZ999-6+BB             1          1          5          9         10 
BKZ999-8+BB             1          1          5          9         10 
BKZ999-10+BB            1          1          5          6         10 
BKZ999-12+BB            1          1          5          6         10 
L8+BKZ-10+BB            1          1          8          6         10 
L5+9+BKZ-10+BB          1          1          5          6         10 
L5+9+BKZ-12+BB          1          1          8          6         10 

Total time for everything: 43.424478 seconds

****************************************************
TestReducersSpeed with m = 1021
Types: Int = NTL::ZZ, Real = NTL::RR
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 0.1e-5.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 3348     12708     39388     89450    168578 
LLL99999             3473     22199    111957    192071    296092 
BKZ99999-10          3559     28649    328600   1075338   2027431 
L5+9+BKZ-10          4270     29537    322483   1119274   2173604 
LLL5+BB              3401     16639    168372   1559641   5828603 
LLL8+BB              3588     21062    164778   1640487   6497942 
LLL99999+BB          3591     25653    208437   1606853   6809049 
BKZ99999-6+BB        3822     31383    327895   1755579   6762201 
BKZ99999-8+BB        3724     31700    370189   1670674   4819753 
BKZ99999-10+BB       3634     32165    404071   1908322   4791739 
BKZ99999-12+BB       3579     32106    432132   1954974   5402812 
BKZ999-6+BB          3547     30975    328564   1797912   6075240 
BKZ999-8+BB          3502     31476    369411   1711019   4885029 
BKZ999-10+BB         3492     32053    399873   1894781   4593722 
BKZ999-12+BB         3451     32083    428646   1945055   5143598 
L8+BKZ-10+BB         4089     30962    388753   1900222   5297273 
L5+9+BKZ-10+BB       4659     32880    400936   1940515   5072855 
L5+9+BKZ-12+BB       4484     32888    429717   2019901   4788574 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5                 10.8       4.28       3.36       3.16       3.04 
LLL99999            10.72       4.14       3.16          3        2.9 
BKZ99999-10         10.72       4.14       3.16        2.9       2.88 
L5+9+BKZ-10         10.72       4.14       3.16        2.9       2.82 
LLL5+BB             10.72       4.14       3.16        2.9       2.78 
LLL8+BB             10.72       4.14       3.16        2.9       2.78 
LLL99999+BB         10.72       4.14       3.16        2.9       2.78 
BKZ99999-6+BB       10.72       4.14       3.16        2.9       2.78 
BKZ99999-8+BB       10.72       4.14       3.16        2.9       2.78 
BKZ99999-10+BB      10.72       4.14       3.16        2.9       2.78 
BKZ99999-12+BB      10.72       4.14       3.16        2.9       2.78 
BKZ999-6+BB         10.72       4.14       3.16        2.9       2.78 
BKZ999-8+BB         10.72       4.14       3.16        2.9       2.78 
BKZ999-10+BB        10.72       4.14       3.16        2.9       2.78 
BKZ999-12+BB        10.72       4.14       3.16        2.9       2.78 
L8+BKZ-10+BB        10.72       4.14       3.16        2.9       2.78 
L5+9+BKZ-10+BB      10.72       4.14       3.16        2.9       2.78 
L5+9+BKZ-12+BB      10.72       4.14       3.16        2.9       2.78 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 6         33       1916      21735      81672 
LLL8+BB                 6         28       1501      22196      88059 
LLL99999+BB             6         26       1375      20209      89946 
BKZ99999-6+BB           6         26       1129      14400      77358 
BKZ99999-8+BB           6         26       1115      11702      45130 
BKZ99999-10+BB          6         26       1056      12052      39125 
BKZ99999-12+BB          6         26       1065      11462      42067 
BKZ999-6+BB             6         26       1129      14886      67701 
BKZ999-8+BB             6         26       1107      12102      46246 
BKZ999-10+BB            6         26       1056      12813      38319 
BKZ999-12+BB            6         26       1070      11388      37994 
L8+BKZ-10+BB            6         26       1063      12158      45155 
L5+9+BKZ-10+BB          6         26       1099      11716      40799 
L5+9+BKZ-12+BB          6         26       1057      12161      29570 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          1          7         18         21 
LLL8+BB                 1          1          5         15         20 
LLL99999+BB             1          1          4         12         20 
BKZ99999-6+BB           1          1          4         13         22 
BKZ99999-8+BB           1          1          4         14         16 
BKZ99999-10+BB          1          1          4          7         18 
BKZ99999-12+BB          1          1          4         14         17 
BKZ999-6+BB             1          1          4         12         20 
BKZ999-8+BB             1          1          4         15         14 
BKZ999-10+BB            1          1          4          6         18 
BKZ999-12+BB            1          1          4         14         23 
L8+BKZ-10+BB            1          1          4         11         14 
L5+9+BKZ-10+BB          1          1          5          9         17 
L5+9+BKZ-12+BB          1          1          5          8         10 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL5+BB                 1          2          5          7          7 
LLL8+BB                 1          1          5          6          6 
LLL99999+BB             1          1          5          6          7 
BKZ99999-6+BB           1          1          5          6          6 
BKZ99999-8+BB           1          1          5          5          5 
BKZ99999-10+BB          1          1          5          7          5 
BKZ99999-12+BB          1          1          5          5          7 
BKZ999-6+BB             1          1          5          6          6 
BKZ999-8+BB             1          1          5          5          5 
BKZ999-10+BB            1          1          5          7          8 
BKZ999-12+BB            1          1          5          5          7 
L8+BKZ-10+BB            1          1          5          5          5 
L5+9+BKZ-10+BB          1          1          5          5          9 
L5+9+BKZ-12+BB          1          1          5          6          7 

Total time for everything: 115.505118 seconds

