****************************************************
Types: Int = long, Real = double

TestBasisConstructionSpeed with m = 1048573
Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5            5616     10612     31245    105933    214270 
LLL9            1022      1748      7734    118093    372838 
LLL99999         717      1122      3102     43973    218399 
LLL99999-R      4182     11058     44942    387797   1118333 
UppTri          1103      1332      2329      6418     11530 
mDualUT          341       462       887      3336      7821 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL9           5.66997e+11    7.13014e+12     6.2494e+13    4.39388e+14    1.02323e+15 
LLL99999        5.6698e+11    7.13014e+12    6.24587e+13    4.32131e+14    9.98214e+14 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 2.76812 seconds


TestBasisConstructionSpeed with m = 1048573
Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5            5117     11176     33024    106674    214031 
LLL9             954      1758      7909    119427    371288 
LLL99999         730      1043      3016     44125    218521 
LLL99999-R      4309     11603     46458    391374   1139340 
UppTri          1029      1302      2399      6273     11501 
mDualUT          301       427       963      3285      7917 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL9           5.66997e+11    7.13014e+12     6.2494e+13    4.39388e+14    1.02323e+15 
LLL99999        5.6698e+11    7.13014e+12    6.24587e+13    4.32131e+14    9.98214e+14 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 2.77692 seconds


****************************************************
Types: Int = NTL::ZZ, Real = double

TestBasisConstructionSpeed with m = 1048573
Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5           11063     21464     56124    151761    272166 
LLL9            2236      4719     16913    156738    358011 
LLL99999        1702      3276      9191     70202    218738 
LLL99999-R      8556     21888     73276    434404    940328 
UppTri          4570      6551     12612     38968     71766 
mDualUT          927      1733      4330     22226     64139 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL9           5.66997e+11    7.13014e+12     6.2494e+13    4.39388e+14    1.02323e+15 
LLL99999        5.6698e+11    7.13014e+12    6.24587e+13    4.32131e+14    9.98214e+14 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 3.15397 seconds


TestBasisConstructionSpeed with m = 1048573
Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5           10674     21594     53680    151924    276161 
LLL9            2207      4646     17004    162722    365844 
LLL99999        1628      3274      9482     74274    224368 
LLL99999-R      8839     22143     75033    449532    968622 
UppTri          4576      6380     12776     38891     73065 
mDualUT          846      1494      4054     21784     63944 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL9           5.66997e+11    7.13014e+12     6.2494e+13    4.39388e+14    1.02323e+15 
LLL99999        5.6698e+11    7.13014e+12    6.24587e+13    4.32131e+14    9.98214e+14 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 3.16708 seconds


****************************************************
Types: Int = NTL::ZZ, Real = xdouble

TestBasisConstructionSpeed with m = 1048573
Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5           36025     74261    214023    619767   1107923 
LLL9            5190     12714     63342    745240   1688071 
LLL99999        3711      7711     27298    263908    929924 
LLL99999-R     32837     86210    321384   2134605   4591602 
UppTri          5310      7544     14237     43967     79398 
mDualUT          972      1957      4735     24801     68700 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76925e+14    1.08841e+15 
LLL9           5.67195e+11    7.13256e+12    6.24994e+13    4.38174e+14    1.02391e+15 
LLL99999        5.6698e+11    7.13082e+12    6.24673e+13    4.34019e+14    1.00614e+15 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 13.3284 seconds


TestBasisConstructionSpeed with m = 1048573
Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5           36901     78107    213405    620460   1104272 
LLL9            5347     13235     63575    741737   1654657 
LLL99999        3948      7892     28685    259816    896989 
LLL99999-R     32710     89254    323808   2083909   4537669 
UppTri          5408      7516     14521     42547     79277 
mDualUT         1298      1773      4457     22772     70839 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76925e+14    1.08841e+15 
LLL9           5.67195e+11    7.13256e+12    6.24994e+13    4.38174e+14    1.02391e+15 
LLL99999        5.6698e+11    7.13082e+12    6.24673e+13    4.34019e+14    1.00614e+15 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 13.0887 seconds


****************************************************
Types: Int = NTL::ZZ, Real = quad_float

TestBasisConstructionSpeed with m = 1048573
Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5           31099     69673    214150    790361   1659058 
LLL9            6156     14698     71023    825905   2127615 
LLL99999        5045      9864     33016    292949   1081899 
LLL99999-R     28772     82879    345309   2505654   6163289 
UppTri          5269      9539     14450     47544     84868 
mDualUT         1150      1922      4672     25617     72250 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL9           5.67195e+11    7.13256e+12    6.24994e+13    4.38183e+14    1.02394e+15 
LLL99999        5.6698e+11    7.13082e+12    6.24673e+13     4.3402e+14    1.00619e+15 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 16.7505 seconds


TestBasisConstructionSpeed with m = 1048573
Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5           33995     76723    235048    830289   1765873 
LLL9            7034     15162     72695    899167   2290478 
LLL99999        6225     11520     35527    327282   1123288 
LLL99999-R     31568     90883    371349   2732942   6601615 
UppTri          5810      7927     15365     51275     92959 
mDualUT         1098      2021      4978     25719     76863 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76704e+14    1.08841e+15 
LLL9           5.67195e+11    7.13256e+12    6.24994e+13    4.38183e+14    1.02394e+15 
LLL99999        5.6698e+11    7.13082e+12    6.24673e+13     4.3402e+14    1.00619e+15 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 17.8972 seconds


****************************************************
Types: Int = NTL::ZZ, Real = NTL::RR

TestBasisConstructionSpeed with m = 1048573
Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5          191344    538111   1892369   6778652  12926668 
LLL9           21086     78086    507326   8720618  22447352 
LLL99999       16437     47715    215323   3391877  12622656 
LLL99999-R    199891    694274   3425295  28693635  64866256 
UppTri          7255     12120     21984     58115    110890 
mDualUT         1428      3085      7508     29627     91384 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76925e+14    1.08841e+15 
LLL9           5.66997e+11    7.13014e+12     6.2494e+13    4.39369e+14    1.02314e+15 
LLL99999        5.6698e+11    7.13014e+12    6.24587e+13     4.3213e+14    9.98265e+14 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 168.822 seconds


TestBasisConstructionSpeed with m = 1048573
Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5          172993    511774   1806295   6443917  11946016 
LLL9           21771     72398    503514   8331321  20939380 
LLL99999       14566     44377    206566   3023164  11969138 
LLL99999-R    186632    663376   3240820  25474108  58777422 
UppTri          6905     10816     20001     55778     97774 
mDualUT         1326      2288      5597     29259     84928 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5           5.68354e+11    7.21003e+12    6.60639e+13    5.76925e+14    1.08841e+15 
LLL9           5.66997e+11    7.13014e+12     6.2494e+13    4.39369e+14    1.02314e+15 
LLL99999        5.6698e+11    7.13014e+12    6.24587e+13     4.3213e+14    9.98265e+14 
LLL99999-R      5.6698e+11    7.13068e+12    6.24642e+13    4.33817e+14    1.00437e+15 

Total time: 154.745 seconds


