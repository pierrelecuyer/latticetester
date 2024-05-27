****************************************************
Types: Int = long, Real = double

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               4285     10247     33253    147147    324167 
LLL8                712      1433      6965    117275    342408 
LLL99               514      1018      4156     96121    482061 
LLL99999            504       795      2042     21023    120216 
LLL99999-new       3678     11285     49495    544196   1712473 
UppTri              578       624      1011      2256      4124 
mDualUT             171       286       453      1282      3148 
LLL5-dual          1915      5488     19258    102770    259866 
LLL8-dual           507      1055      4290     58936    135470 
LLL99-dual          481       873      3037     53032    221360 
LLL99999-dual       423       774      1860     13096     58821 
LLL99999-new       1852      5740     24318    230484    691536 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL8              5.66997e+11    7.13174e+12    6.25167e+13    4.48781e+14    1.04244e+15 
LLL99              5.6698e+11    7.12985e+12    6.24462e+13    4.32408e+14     1.0006e+15 
LLL99999           5.6698e+11    7.12985e+12    6.24462e+13    4.31436e+14     9.9619e+14 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              271343          36505           8736           4981           4320 
LLL8-dual              270865          36134           8302           3943           3503 
LLL99-dual             270865          36134           8287           3817           3344 
LLL99999-dual          270865          36134           8287           3817           3326 
LLL99999-new           270865          36134           8296           3841           3389 

Total time: 5.97753 seconds


Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               4181     10117     33059    147267    324240 
LLL8                678      1419      6959    117427    342395 
LLL99               507      1032      4150     96151    482069 
LLL99999            472       796      2019     20952    120282 
LLL99999-new       3703     11279     49509    544197   1712739 
UppTri              474       612       992      2088      4101 
mDualUT             199       251       435      1295      3122 
LLL5-dual          1925      5480     19270    103130    260031 
LLL8-dual           514      1077      4336     59076    135507 
LLL99-dual          459       856      3004     53123    221390 
LLL99999-dual       436       770      1850     13097     58942 
LLL99999-new       1868      5753     24338    231123    691783 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL8              5.66997e+11    7.13174e+12    6.25167e+13    4.48781e+14    1.04244e+15 
LLL99              5.6698e+11    7.12985e+12    6.24462e+13    4.32408e+14     1.0006e+15 
LLL99999           5.6698e+11    7.12985e+12    6.24462e+13    4.31436e+14     9.9619e+14 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              271343          36505           8736           4981           4320 
LLL8-dual              270865          36134           8302           3943           3503 
LLL99-dual             270865          36134           8287           3817           3344 
LLL99999-dual          270865          36134           8287           3817           3326 
LLL99999-new           270865          36134           8296           3841           3389 

Total time: 5.96562 seconds


****************************************************
Types: Int = NTL::ZZ, Real = double

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               6567     13452     37379    118076    230871 
LLL8               1463      2787     10302     92076    207200 
LLL99              1115      2075      6818     74019    241349 
LLL99999            891      1734      4539     26512     86445 
LLL99999-new       5466     14343     53684    362495    826793 
UppTri             2617      3604      6564     18441     33438 
mDualUT             607      1003      2236     10327     29032 
LLL5-dual          5200     11632     30593     79727    140022 
LLL8-dual          1219      2878     10985     89654    165059 
LLL99-dual          965      2187      7435     76740    228285 
LLL99999-dual       866      1822      4926     26800     75676 
LLL99999-new       4591     11788     40358    236479    468735 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL8              5.66997e+11    7.13174e+12    6.25167e+13    4.48781e+14    1.04244e+15 
LLL99              5.6698e+11    7.12985e+12    6.24462e+13    4.32408e+14     1.0006e+15 
LLL99999           5.6698e+11    7.12985e+12    6.24462e+13    4.31436e+14     9.9619e+14 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              539778          68409          16137           9018           7967 
LLL8-dual              538847          67503          15306           7032           6250 
LLL99-dual             538847          67503          15270           6773           5925 
LLL99999-dual          538847          67503          15270           6770           5894 
LLL99999-new           538847          67498          15288           6831           6001 

Total time: 4.36742 seconds


Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               6398     13239     37118    117666    230661 
LLL8               1474      2777     10569     92652    208521 
LLL99              1085      2109      6863     75486    245090 
LLL99999            905      1764      4565     27141     87787 
LLL99999-new       5462     14486     54432    368774    841575 
UppTri             2663      3643      6737     18787     33705 
mDualUT             513       864      2123     10865     31288 
LLL5-dual          5021     11510     30612     80267    140668 
LLL8-dual          1236      2929     10830     90111    166323 
LLL99-dual          962      2202      7456     76878    229059 
LLL99999-dual       854      1836      4941     26455     73788 
LLL99999-new       4565     11655     40854    239252    471224 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL8              5.66997e+11    7.13174e+12    6.25167e+13    4.48781e+14    1.04244e+15 
LLL99              5.6698e+11    7.12985e+12    6.24462e+13    4.32408e+14     1.0006e+15 
LLL99999           5.6698e+11    7.12985e+12    6.24462e+13    4.31436e+14     9.9619e+14 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              539778          68409          16137           9018           7967 
LLL8-dual              538847          67503          15306           7032           6250 
LLL99-dual             538847          67503          15270           6773           5925 
LLL99999-dual          538847          67503          15270           6770           5894 
LLL99999-new           538847          67498          15288           6831           6001 

Total time: 4.36957 seconds


****************************************************
Types: Int = NTL::ZZ, Real = xdouble

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              20999     45706    128746    368547    633430 
LLL8               3015      7225     32535    328875    651354 
LLL99              2315      5034     19490    259746    817218 
LLL99999           1882      3728     10490     55505    166519 
LLL99999-new      18666     51783    194449   1283385   2761518 
UppTri             2921      3869      6967     19440     35317 
mDualUT             627      1063      2322     10814     30352 
LLL5-dual         16423     38298    105501    254617    402458 
LLL8-dual          2621      7058     33086    303755    498999 
LLL99-dual         2082      5052     20153    262348    779744 
LLL99999-dual      1766      3827     10914     57608    161302 
LLL99999-new      15177     40954    144448    835748   1581393 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76925e+14    1.08841e+15 
LLL8              5.67181e+11    7.13297e+12    6.26165e+13    4.50306e+14    1.04391e+15 
LLL99              5.6698e+11    7.13082e+12    6.24909e+13    4.34266e+14    1.00553e+15 
LLL99999           5.6698e+11    7.13082e+12    6.24868e+13    4.34138e+14    1.00444e+15 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              539778          68409          16137           9018           7967 
LLL8-dual              538847          67498          15340           7073           6348 
LLL99-dual             538847          67498          15282           6813           5988 
LLL99999-dual          538847          67498          15282           6806           5972 
LLL99999-new           538847          67498          15288           6831           6001 

Total time: 13.6851 seconds


Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              20382     45374    127954    365118    629412 
LLL8               2960      7221     32519    327295    646234 
LLL99              2215      5023     19454    259430    812495 
LLL99999           1804      3734     10460     55291    164177 
LLL99999-new      18501     51704    194177   1280884   2751748 
UppTri             2892      3820      6884     19378     34994 
mDualUT             574       937      2233     11289     32581 
LLL5-dual         16177     37963    104771    253487    400783 
LLL8-dual          2622      7041     32822    302169    496308 
LLL99-dual         2063      5021     20052    260790    774555 
LLL99999-dual      1762      3796     10886     56474    157652 
LLL99999-new      14992     40652    144093    832324   1571393 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76925e+14    1.08841e+15 
LLL8              5.67181e+11    7.13297e+12    6.26165e+13    4.50306e+14    1.04391e+15 
LLL99              5.6698e+11    7.13082e+12    6.24909e+13    4.34266e+14    1.00553e+15 
LLL99999           5.6698e+11    7.13082e+12    6.24868e+13    4.34138e+14    1.00444e+15 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              539778          68409          16137           9018           7967 
LLL8-dual              538847          67498          15340           7073           6348 
LLL99-dual             538847          67498          15282           6813           5988 
LLL99999-dual          538847          67498          15282           6806           5972 
LLL99999-new           538847          67498          15288           6831           6001 

Total time: 13.5737 seconds


****************************************************
Types: Int = NTL::ZZ, Real = quad_float

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              11733     26104     74450    256491    531202 
LLL8               2678      5580     21722    203039    470791 
LLL99              2229      4360     14565    160611    534674 
LLL99999           1986      3606      9471     46344    135605 
LLL99999-new      11249     30218    116154    809799   1948895 
UppTri             2551      3655      6582     18637     33439 
mDualUT             550       986      2235     10320     28997 
LLL5-dual          9914     21773     57290    150836    273822 
LLL8-dual          2441      5575     22268    187939    328791 
LLL99-dual         2097      4376     15228    166556    494144 
LLL99999-dual      1918      3677      9883     48622    133422 
LLL99999-new       9630     23864     81862    491804    970468 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL8              5.67181e+11    7.13297e+12    6.26165e+13    4.50262e+14    1.04388e+15 
LLL99              5.6698e+11    7.13082e+12    6.24909e+13    4.34169e+14    1.00555e+15 
LLL99999           5.6698e+11    7.13082e+12    6.24868e+13    4.34066e+14     1.0044e+15 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              539778          68409          16137           9021           7967 
LLL8-dual              538847          67498          15344           7071           6354 
LLL99-dual             538847          67498          15282           6813           5988 
LLL99999-dual          538847          67498          15282           6806           5972 
LLL99999-new           538847          67498          15288           6831           6002 

Total time: 9.16777 seconds


Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              11494     25795     73980    254007    528428 
LLL8               2689      5574     21734    201973    470924 
LLL99              2232      4352     14511    160782    536565 
LLL99999           1990      3633      9450     46689    135092 
LLL99999-new      11223     30317    116003    810624   1953231 
UppTri             2469      3574      6523     18644     33281 
mDualUT             449       851      2084     10620     30963 
LLL5-dual          9713     21494     56541    150236    273209 
LLL8-dual          2468      5594     21887    187436    328490 
LLL99-dual         2111      4404     15027    165684    491916 
LLL99999-dual      1920      3656      9753     47789    130075 
LLL99999-new       9506     23775     81270    491501    965838 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL8              5.67181e+11    7.13297e+12    6.26165e+13    4.50262e+14    1.04388e+15 
LLL99              5.6698e+11    7.13082e+12    6.24909e+13    4.34169e+14    1.00555e+15 
LLL99999           5.6698e+11    7.13082e+12    6.24868e+13    4.34066e+14     1.0044e+15 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              539778          68409          16137           9021           7967 
LLL8-dual              538847          67498          15344           7071           6354 
LLL99-dual             538847          67498          15282           6813           5988 
LLL99999-dual          538847          67498          15282           6806           5972 
LLL99999-new           538847          67498          15288           6831           6002 

Total time: 9.11249 seconds


