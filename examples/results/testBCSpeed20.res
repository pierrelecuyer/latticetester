****************************************************
Types: Int = long, Real = double

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               4396     10413     33813    147678    325905 
LLL8                824      1503      7072    118212    344804 
LLL99               623      1037      4160     97013    484988 
LLL99999            536       821      2057     21373    121697 
LLL99999-new       3754     11332     50083    548472   1722832 
UppTri              501       622       927      2020      4012 
mDualUT             209       272       452      1307      3118 
LLL5-dual          2012      5522     19423    103345    261220 
LLL8-dual           544      1068      4343     59408    136982 
LLL99-dual          494       867      3032     53497    223310 
LLL99999-dual       480       765      1872     13439     60058 
LLL99999-dnew      1895      5747     24551    232049    694604 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
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
LLL99999-dnew          270865          36134           8296           3841           3389 

Total time: 6.00407 seconds


****************************************************
Types: Int = NTL::ZZ, Real = double

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               6454     13150     36752    115737    225875 
LLL8               1428      2736     10261     90592    203944 
LLL99              1052      1991      6748     73260    239376 
LLL99999            877      1664      4479     26314     85776 
LLL99999-new       5342     14313     53776    359121    824510 
UppTri             2500      3473      6498     18118     33303 
mDualUT             456       841      2103     10889     32153 
LLL5-dual          4995     11378     30047     78292    140922 
LLL8-dual          1195      2823     10671     88285    164148 
LLL99-dual          930      2153      7297     75193    225996 
LLL99999-dual       825      1744      4882     25716     72891 
LLL99999-dnew      4420     11559     40198    233219    468637 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
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
LLL99999-dnew          538847          67498          15288           6831           6001 

Total time: 4.29407 seconds


****************************************************
Types: Int = NTL::ZZ, Real = xdouble

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              19221     43241    122836    351060    603249 
LLL8               2819      6826     31191    314013    620449 
LLL99              2132      4734     18717    248915    779389 
LLL99999           1787      3563     10134     53397    157545 
LLL99999-new      17719     49691    186262   1229097   2646469 
UppTri             2710      3618      6608     18233     32758 
mDualUT             525       889      2169     10784     31222 
LLL5-dual         15428     36249    100438    241876    382827 
LLL8-dual          2538      6730     31456    288498    473814 
LLL99-dual         1991      4798     19279    249599    738533 
LLL99999-dual      1727      3650     10479     54538    151032 
LLL99999-dnew     14346     39033    138003    794357   1502005 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
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
LLL99999-dnew          538847          67498          15288           6831           6001 

Total time: 13.0076 seconds


****************************************************
Types: Int = NTL::ZZ, Real = quad_float

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              11201     25173     72592    249558    516528 
LLL8               2602      5424     21208    197435    455925 
LLL99              2184      4218     14193    156569    520265 
LLL99999           1945      3509      9228     45396    131678 
LLL99999-new      11012     29635    113856    796581   1927552 
UppTri             2406      3445      6410     18151     32917 
mDualUT             440       825      2041     10498     30930 
LLL5-dual          9480     20983     55409    146729    270261 
LLL8-dual          2400      5437     21404    183021    322100 
LLL99-dual         2060      4299     14706    161891    483026 
LLL99999-dual      1885      3517      9527     46777    128139 
LLL99999-dnew      9246     23167     79770    479400    948041 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
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
LLL99999-dnew          538847          67498          15288           6831           6002 

Total time: 8.93549 seconds


****************************************************
Types: Int = NTL::ZZ, Real = NTL::RR

TestBasisConstructionSpeed with m = 1048573
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              75076    209343    749954   2633612   4802095 
LLL8               8581     27448    175087   2400292   5526790 
LLL99              6241     19636    108808   1955108   6971480 
LLL99999           4995     13827     55527    544954   2182031 
LLL99999-new      77474    269447   1323004  10659878  24091581 
UppTri             3084      3983      7129     19070     34299 
mDualUT             600       947      2262     10840     31210 
LLL5-dual         56097    149437    475360   1311291   2421321 
LLL8-dual          8098     27016    158932   1858844   3448582 
LLL99-dual         5997     19369     99273   1604261   5340611 
LLL99999-dual      4895     13443     49379    367736   1430478 
LLL99999-dnew     57690    181036    771131   5294149  10445348 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
 dim:                4              6             10             20             30  

LLL5              5.68354e+11    7.21003e+12    6.60639e+13    5.76925e+14    1.08841e+15 
LLL8              5.66997e+11    7.13174e+12    6.25167e+13    4.48786e+14    1.04223e+15 
LLL99              5.6698e+11    7.12985e+12    6.24462e+13    4.32412e+14    1.00055e+15 
LLL99999           5.6698e+11    7.12985e+12    6.24462e+13    4.31437e+14    9.96135e+14 
LLL99999-new       5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 
LLL5-dual              539778          68409          16137           9018           7967 
LLL8-dual              538847          67503          15306           7028           6248 
LLL99-dual             538847          67503          15270           6774           5925 
LLL99999-dual          538847          67503          15270           6771           5892 
LLL99999-dnew          538847          67498          15288           6831           6002 

Total time: 100.682 seconds


