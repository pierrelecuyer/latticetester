
========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Num. dimensions:       5        10        20        30        40  

LLL5                  956      3361     11227     20046     33528 
LLL99999              825      4476     33664     91986    161909 
BKZ99999-10           715      4352     39032    135673    339928 
L5+9+BKZ-10           856      3941     24854    100647    282574 
LLL99999+BB           873      5176     39847    229469  45392198 
BKZ99999-6+BB         767      5105     43503    203873  14335919 
BKZ99999-8+BB         751      4940     44373    210764  12329055 
BKZ99999-10+BB        693      4798     45155    212110  10659265 
BKZ99999-12+BB        666      4685     45570    222533   9875383 
BKZ999-6+BB           671      4657     43639    207392  13418660 
BKZ999-8+BB           659      4629     44734    212458  13854416 
BKZ999-10+BB          653      4589     45223    217717  11497435 
BKZ999-12+BB          650      4530     45737    218667  10917121 
L8+BKZ-10+BB          940      5077     35524    180732  13129788 
L5+9+BKZ-10+BB        926      4712     31885    178864  12785387 
L5+9+BKZ-12+BB        790      4487     32008    192186  10475289 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5           1.040622168e+19 4.396980005e+21 2.306416649e+23 9.0331178e+23 1.208432454e+24 
LLL99999       1.040622168e+19 4.205880845e+21 1.146463186e+23 4.33286451e+23 9.83282428e+23 
BKZ99999-10    1.040622168e+19 4.205880845e+21 1.144617381e+23 4.167450442e+23 8.663775181e+23 
L5+9+BKZ-10    1.040622168e+19 4.205880845e+21 1.143017526e+23 4.129312932e+23 8.747641965e+23 
LLL99999+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-6+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-8+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-6+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-8+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-12+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L8+BKZ-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        654      94287   35676811 
BKZ99999-6+BB           5         16        566      51105   11069270 
BKZ99999-8+BB           5         16        553      52183    9426568 
BKZ99999-10+BB          5         16        573      45066    8071678 
BKZ99999-12+BB          5         16        575      44060    7434917 
BKZ999-6+BB             5         16        553      54488   10298185 
BKZ999-8+BB             5         16        559      53281   10650118 
BKZ999-10+BB            5         16        577      46995    8732897 
BKZ999-12+BB            5         16        577      46194    8237243 
L8+BKZ-10+BB            5         16        576      48460   10116347 
L5+9+BKZ-10+BB          5         16        566      46857    9832036 
L5+9+BKZ-12+BB          5         16        559      44873    7941386 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          4 
BKZ99999-6+BB           1          1          1          1          2 
BKZ99999-8+BB           1          1          1          1          2 
BKZ99999-10+BB          1          1          1          1          2 
BKZ99999-12+BB          1          1          1          1          2 
BKZ999-6+BB             1          1          1          1          2 
BKZ999-8+BB             1          1          1          1          2 
BKZ999-10+BB            1          1          1          1          2 
BKZ999-12+BB            1          1          1          1          2 
L8+BKZ-10+BB            1          1          1          1          2 
L5+9+BKZ-10+BB          1          1          1          1          2 
L5+9+BKZ-12+BB          1          1          1          1          2 

Total time for everything: 183.033676 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Num. dimensions:       5        10        20        30        40  

LLL5                  667      2119      5420      7349     11160 
LLL99999              586      2532     13435     32540     46944 
BKZ99999-10           575      2509     17465     67577    168436 
L5+9+BKZ-10           663      2768     17659     65919    162216 
LLL99999+BB           660      3207     19430    208532  73554207 
BKZ99999-6+BB         599      3177     22393    144172  19888233 
BKZ99999-8+BB         551      3031     23161    149543  15309212 
BKZ99999-10+BB        528      2926     23701    154806  16261578 
BKZ99999-12+BB        517      2871     24310    152247  14199518 
BKZ999-6+BB           520      2821     22428    148173  21815602 
BKZ999-8+BB           503      2787     23103    152991  16436182 
BKZ999-10+BB          503      2767     23586    156544  17473434 
BKZ999-12+BB          494      2731     24385    150444  13184568 
L8+BKZ-10+BB          651      3318     22698    147323  13987731 
L5+9+BKZ-10+BB        742      3496     23571    157068  13822829 
L5+9+BKZ-12+BB        674      3355     24001    160298  12090483 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5              36783.7        255      46.82      40.64      36.38 
LLL99999         36395.12     231.94      25.38      14.84       13.5 
BKZ99999-10      36395.12      231.4      25.16      14.28      11.94 
L5+9+BKZ-10      36395.12      231.4      25.14      14.26      11.86 
LLL99999+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ99999-6+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-8+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-10+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ99999-12+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ999-6+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-8+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-10+BB     36395.12      231.4      25.14       14.1       11.5 
BKZ999-12+BB     36395.12      231.4      25.14       14.1       11.5 
L8+BKZ-10+BB     36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-10+BB   36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-12+BB   36395.12      231.4      25.14       14.1       11.5 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        793     129459   58161623 
BKZ99999-6+BB           5         16        726      61994   15542899 
BKZ99999-8+BB           5         16        684      61143   11909034 
BKZ99999-10+BB          5         16        701      57957   12646861 
BKZ99999-12+BB          5         16        693      53608   11006050 
BKZ999-6+BB             5         16        716      66232   17066547 
BKZ999-8+BB             5         16        690      63917   12828353 
BKZ999-10+BB            5         16        697      58726   13659907 
BKZ999-12+BB            5         16        680      52185   10183027 
L8+BKZ-10+BB            5         16        687      56562   10830150 
L5+9+BKZ-10+BB          5         16        680      61158   10721944 
L5+9+BKZ-12+BB          5         16        668      56555    9309786 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          7 
BKZ99999-6+BB           1          1          1          2          6 
BKZ99999-8+BB           1          1          1          2          6 
BKZ99999-10+BB          1          1          1          1          5 
BKZ99999-12+BB          1          1          1          2          5 
BKZ999-6+BB             1          1          1          2          6 
BKZ999-8+BB             1          1          1          2          6 
BKZ999-10+BB            1          1          1          1          5 
BKZ999-12+BB            1          1          1          2          5 
L8+BKZ-10+BB            1          1          1          2          5 
L5+9+BKZ-10+BB          1          1          1          2          5 
L5+9+BKZ-12+BB          1          1          1          2          5 

Total time for everything: 250.969225 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = xdouble
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Num. dimensions:       5        10        20        30        40  

LLL5                 4029     17141     62489    108748    168072 
LLL99999             4016     26857    218400    623810   1042295 
BKZ99999-10          3962     27652    257715    950614   2380315 
L5+9+BKZ-10          3848     21437    158712    678249   1888529 
LLL99999+BB          3990     28155    231451   1701865 455261486 
BKZ99999-6+BB        3833     28812    259173   1426961 141766655 
BKZ99999-8+BB        3668     28847    264728   1478653 121930257 
BKZ99999-10+BB       3567     28867    270204   1476399 104974919 
BKZ99999-12+BB       3495     28690    273164   1557405  97583024 
BKZ999-6+BB          3454     28627    258783   1463121 133718491 
BKZ999-8+BB          3417     28722    267626   1495345 137057386 
BKZ999-10+BB         3401     28653    270835   1526742 113611111 
BKZ999-12+BB         3394     28748    274204   1530173 107447282 
L8+BKZ-10+BB         3927     25915    204279   1279894 116801509 
L5+9+BKZ-10+BB       3912     22809    170665   1256153 143592441 
L5+9+BKZ-12+BB       3754     22494    173304   1274519 105884863 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5           1.040622168e+19 4.396980005e+21 2.306416649e+23 9.0331178e+23 1.208432454e+24 
LLL99999       1.040622168e+19 4.205880845e+21 1.146463186e+23 4.33286451e+23 9.83282428e+23 
BKZ99999-10    1.040622168e+19 4.205880845e+21 1.144617381e+23 4.167450442e+23 8.663775181e+23 
L5+9+BKZ-10    1.040622168e+19 4.205880845e+21 1.143917209e+23 4.166116336e+23 8.805377564e+23 
LLL99999+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-6+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-8+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-6+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-8+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-12+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L8+BKZ-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        654      94287   35676811 
BKZ99999-6+BB           5         16        566      51105   11069270 
BKZ99999-8+BB           5         16        553      52183    9426568 
BKZ99999-10+BB          5         16        573      45066    8071678 
BKZ99999-12+BB          5         16        575      44060    7434917 
BKZ999-6+BB             5         16        553      54488   10298185 
BKZ999-8+BB             5         16        559      53281   10650118 
BKZ999-10+BB            5         16        577      46995    8732897 
BKZ999-12+BB            5         16        577      46194    8237243 
L8+BKZ-10+BB            5         16        573      45842    8987982 
L5+9+BKZ-10+BB          5         16        562      49378   11235264 
L5+9+BKZ-12+BB          5         16        564      45733    8170195 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          4 
BKZ99999-6+BB           1          1          1          1          2 
BKZ99999-8+BB           1          1          1          1          2 
BKZ99999-10+BB          1          1          1          1          2 
BKZ99999-12+BB          1          1          1          1          2 
BKZ999-6+BB             1          1          1          1          2 
BKZ999-8+BB             1          1          1          1          2 
BKZ999-10+BB            1          1          1          1          2 
BKZ999-12+BB            1          1          1          1          2 
L8+BKZ-10+BB            1          1          1          1          2 
L5+9+BKZ-10+BB          1          1          1          1          2 
L5+9+BKZ-12+BB          1          1          1          1          2 

Total time for everything: 1809.058048 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = xdouble
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Num. dimensions:       5        10        20        30        40  

LLL5                 2449      8604     25584     35590     47034 
LLL99999             2239     11529     74700    200807    296395 
BKZ99999-10          2281     12159    101027    431353   1147470 
L5+9+BKZ-10          2387     12015     91761    393732   1053223 
LLL99999+BB          2309     12530     88468   1655998 731590288 
BKZ99999-6+BB        2257     13050    104939   1029850 197843170 
BKZ99999-8+BB        2180     13024    109581   1072028 150750059 
BKZ99999-10+BB       2135     12884    113421   1092629 160217041 
BKZ99999-12+BB       2093     12788    117614   1064620 139152438 
BKZ999-6+BB          2060     12662    105094   1064917 216935477 
BKZ999-8+BB          2045     12655    109731   1100217 163569291 
BKZ999-10+BB         2037     12682    113368   1106535 172210157 
BKZ999-12+BB         2028     12641    117833   1052963 129096865 
L8+BKZ-10+BB         2304     13004    103312   1062214 145521350 
L5+9+BKZ-10+BB       2468     13221    103880   1076645 161227369 
L5+9+BKZ-12+BB       2394     12825    106390   1100690 118705333 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5              36783.7        255      46.82      40.64      36.38 
LLL99999         36395.12     231.94      25.38      14.84       13.5 
BKZ99999-10      36395.12      231.4      25.16      14.28      11.94 
L5+9+BKZ-10      36395.12      231.4      25.14      14.34      11.96 
LLL99999+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ99999-6+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-8+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-10+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ99999-12+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ999-6+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-8+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-10+BB     36395.12      231.4      25.14       14.1       11.5 
BKZ999-12+BB     36395.12      231.4      25.14       14.1       11.5 
L8+BKZ-10+BB     36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-10+BB   36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-12+BB   36395.12      231.4      25.14       14.1       11.5 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        793     129459   58161623 
BKZ99999-6+BB           5         16        726      61994   15542899 
BKZ99999-8+BB           5         16        684      61897   11843113 
BKZ99999-10+BB          5         16        690      57957   12646318 
BKZ99999-12+BB          5         16        693      53608   10942933 
BKZ999-6+BB             5         16        716      66232   17066547 
BKZ999-8+BB             5         16        690      64671   12828353 
BKZ999-10+BB            5         16        686      58726   13659907 
BKZ999-12+BB            5         16        680      53002   10134442 
L8+BKZ-10+BB            5         16        667      59441   11457671 
L5+9+BKZ-10+BB          5         16        678      60102   12710804 
L5+9+BKZ-12+BB          5         16        673      58153    9279093 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          7 
BKZ99999-6+BB           1          1          1          2          6 
BKZ99999-8+BB           1          1          1          2          6 
BKZ99999-10+BB          1          1          1          1          5 
BKZ99999-12+BB          1          1          1          2          5 
BKZ999-6+BB             1          1          1          2          6 
BKZ999-8+BB             1          1          1          2          6 
BKZ999-10+BB            1          1          1          1          5 
BKZ999-12+BB            1          1          1          2          5 
L8+BKZ-10+BB            1          1          1          2          5 
L5+9+BKZ-10+BB          1          1          1          2          6 
L5+9+BKZ-12+BB          1          1          1          2          5 

Total time for everything: 2505.848669 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Num. dimensions:       5        10        20        30        40  

LLL5                 3004     17198     85309    186851    338137 
LLL99999             3321     28918    284509    886385   1682812 
BKZ99999-10          3318     30354    343315   1382139   3814309 
L5+9+BKZ-10          3319     22156    213510   1001975   3066519 
LLL99999+BB          3514     30084    296987   1812136 383589324 
BKZ99999-6+BB        3457     31231    338353   1715721 122404752 
BKZ99999-8+BB        3392     31309    346355   1788986 104553132 
BKZ99999-10+BB       3368     31317    354752   1837760  90824462 
BKZ99999-12+BB       3368     31256    360290   1968582  84498076 
BKZ999-6+BB          3353     30967    337354   1744010 114189955 
BKZ999-8+BB          3346     31178    349700   1807016 117887620 
BKZ999-10+BB         3339     31198    355701   1898421  98010875 
BKZ999-12+BB         3340     31218    361484   1909616  93094557 
L8+BKZ-10+BB         3614     27354    268254   1564256 100441595 
L5+9+BKZ-10+BB       3517     23336    225240   1503933 122927505 
L5+9+BKZ-12+BB       3410     23276    229771   1580076  91054281 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5           1.040622168e+19 4.396980005e+21 2.306416649e+23 9.0331178e+23 1.208432454e+24 
LLL99999       1.040622168e+19 4.205880845e+21 1.146463186e+23 4.33286451e+23 9.83282428e+23 
BKZ99999-10    1.040622168e+19 4.205880845e+21 1.144617381e+23 4.167450442e+23 8.663775181e+23 
L5+9+BKZ-10    1.040622168e+19 4.205880845e+21 1.143917209e+23 4.166116336e+23 8.805377564e+23 
LLL99999+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-6+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-8+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-6+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-8+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-12+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L8+BKZ-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        654      94287   35676811 
BKZ99999-6+BB           5         16        566      51105   11069270 
BKZ99999-8+BB           5         16        553      52183    9426568 
BKZ99999-10+BB          5         16        573      45066    8071678 
BKZ99999-12+BB          5         16        575      44060    7434917 
BKZ999-6+BB             5         16        553      54488   10293946 
BKZ999-8+BB             5         16        559      53194   10650118 
BKZ999-10+BB            5         16        577      46995    8732897 
BKZ999-12+BB            5         16        577      46194    8237243 
L8+BKZ-10+BB            5         16        573      45842    8996138 
L5+9+BKZ-10+BB          5         16        562      49758   11185848 
L5+9+BKZ-12+BB          5         16        564      45664    8079564 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          4 
BKZ99999-6+BB           1          1          1          1          2 
BKZ99999-8+BB           1          1          1          1          2 
BKZ99999-10+BB          1          1          1          1          2 
BKZ99999-12+BB          1          1          1          1          2 
BKZ999-6+BB             1          1          1          1          2 
BKZ999-8+BB             1          1          1          1          2 
BKZ999-10+BB            1          1          1          1          2 
BKZ999-12+BB            1          1          1          1          2 
L8+BKZ-10+BB            1          1          1          1          2 
L5+9+BKZ-10+BB          1          1          1          1          2 
L5+9+BKZ-12+BB          1          1          1          1          2 

Total time for everything: 1562.247491 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Num. dimensions:       5        10        20        30        40  

LLL5                 1687      6311     19558     28919     42665 
LLL99999             1666      8976     61815    164966    241379 
BKZ99999-10          1754      9953     93371    413661   1148828 
L5+9+BKZ-10          1971      9758     85429    380912   1061926 
LLL99999+BB          1865      9981     74400   1428497 624846614 
BKZ99999-6+BB        1865     10808     93287    910808 168856763 
BKZ99999-8+BB        1813     10821     99225    957818 134464712 
BKZ99999-10+BB       1795     10815    105040    994100 137167538 
BKZ99999-12+BB       1770     10753    111938    993329 118490368 
BKZ999-6+BB          1765     10546     92805    936895 184806976 
BKZ999-8+BB          1761     10608     99261    984691 141696247 
BKZ999-10+BB         1757     10649    104828   1010661 144482064 
BKZ999-12+BB         1748     10622    112044    978386 114557814 
L8+BKZ-10+BB         1989     10787     97272    969403 125353136 
L5+9+BKZ-10+BB       2120     10769     97133    981053 139725726 
L5+9+BKZ-12+BB       2070     10678    101554   1022903 103594240 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5              36783.7        255      46.82      40.64      36.38 
LLL99999         36395.12     231.94      25.38      14.84       13.5 
BKZ99999-10      36395.12      231.4      25.16      14.28      11.96 
L5+9+BKZ-10      36395.12      231.4      25.14      14.34      11.92 
LLL99999+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ99999-6+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-8+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-10+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ99999-12+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ999-6+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-8+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-10+BB     36395.12      231.4      25.14       14.1       11.5 
BKZ999-12+BB     36395.12      231.4      25.14       14.1       11.5 
L8+BKZ-10+BB     36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-10+BB   36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-12+BB   36395.12      231.4      25.14       14.1       11.5 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        793     129460   58178894 
BKZ99999-6+BB           5         16        727      61892   15470688 
BKZ99999-8+BB           5         16        683      61332   12334879 
BKZ99999-10+BB          5         16        681      58005   12572744 
BKZ99999-12+BB          5         16        669      54628   10821560 
BKZ999-6+BB             5         16        716      65612   16916204 
BKZ999-8+BB             5         16        689      64748   12989651 
BKZ999-10+BB            5         16        686      58810   13281881 
BKZ999-12+BB            5         16        669      52984   10469421 
L8+BKZ-10+BB            5         16        671      60083   11459510 
L5+9+BKZ-10+BB          5         16        689      60103   12842552 
L5+9+BKZ-12+BB          5         16        674      58309    9463081 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          3          8 
BKZ99999-6+BB           1          1          1          2          6 
BKZ99999-8+BB           1          1          1          2          6 
BKZ99999-10+BB          1          1          1          2          5 
BKZ99999-12+BB          1          1          1          2          5 
BKZ999-6+BB             1          1          1          2          6 
BKZ999-8+BB             1          1          1          2          6 
BKZ999-10+BB            1          1          1          2          6 
BKZ999-12+BB            1          1          1          2          5 
L8+BKZ-10+BB            1          1          1          2          6 
L5+9+BKZ-10+BB          1          1          1          2          5 
L5+9+BKZ-12+BB          1          1          1          2          5 

Total time for everything: 2155.459792 seconds

