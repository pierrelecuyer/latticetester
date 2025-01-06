****************************************************
TestBasisConstructSpeedLLL with m = 1048573
Types: Int = long, Real = double
Number of replications (different multipliers a): 1000
Total time: 26.6629 seconds.

Timings for different methods, in basic clock units (microseconds) 
 dim:                5        10        20        30        40        50        60        70  

LLL5               5414     19696     59537    106355    173264    255401    367868    499011 
LLL8                867      4268     47463    109142    182840    280537    409545    555017 
LLL99               602      2546     38875    153903    277207    431279    617108    870645 
LLL99999            518      1178      8909     39689     43129     69648    110552    153408 
LLL99999-pnew      5002     27731    206175    531115    917430   1392377   1991637   2722665 
UppTri              573       982      2209      4310      6654      9402     12652     16491 
mDualUT             210       386      1364      3334      6489     11172     17557     25991 
LLL5-dual          4645     19654     69850    143580    260318    433668    685026   1011708 
LLL8-dual           702      4015     44515     87866    143712    225178    340226    493654 
LLL99-dual          565      2444     38625    142537    228304    346239    500512    706069 
LLL99999-dual       508      1141      7451     32392     71655    125612    192617    291967 
LLL99999-dnew      4484     24679    169153    413214    679676   1048147   1555293   2207682 

Sums of square lengths of shortest basis vectors (must be the same for all flexible types): 
 dim:                5            10            20            30            40            50            60            70 

LLL5             2.55405e+12   6.60639e+13   5.76704e+14   1.08841e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL8             2.53724e+12   6.25167e+13   4.48781e+14   1.04244e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99            2.53662e+12   6.24462e+13   4.32408e+14    1.0006e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999         2.53662e+12   6.24462e+13   4.31436e+14    9.9619e+14   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999-pnew     2.5372e+12   6.24642e+13   4.33817e+14   1.00437e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL5-dual             152578         16137          9018          7967          7404          7083          6801          6593 
LLL8-dual             151770         15306          7021          6249          5888          5675          5524          5400 
LLL99-dual            151760         15270          6776          5922          5558          5340          5201          5094 
LLL99999-dual         151760         15270          6773          5892          5529          5311          5174          5065 
LLL99999-dnew         151782         15288          6831          6001          5655          5484          5357          5249 

****************************************************
TestBasisConstructSpeedLLL with m = 1048573
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 1000
Total time: 18.5455 seconds.

Timings for different methods, in basic clock units (microseconds) 
 dim:                5        10        20        30        40        50        60        70  

LLL5               9220     29552     74496    126434    191720    268629    367980    487655 
LLL8               2047      8566     59466    114287    174591    246840    338104    440091 
LLL99              1554      5829     48374    134812    213292    302129    408204    540765 
LLL99999           1299      4046     19700     53852     64528     93987    135335    180886 
LLL99999-pnew      7927     37977    203095    416276    642128    897106   1205708   1559468 
UppTri             3763      7085     19249     34602     48036     61776     77209     96440 
mDualUT             777      2193     10897     31662     70352    132911    224344    349253 
LLL5-dual          7327     24360     52372     79301    115000    161420    222307    303459 
LLL8-dual          1910      9121     58101     96339    138933    192306    264202    357165 
LLL99-dual         1461      6435     50343    128834    181199    244219    326501    431697 
LLL99999-dual      1267      4460     19889     46256     78103    118524    171125    240381 
LLL99999-dnew      6581     30542    142205    252629    312484    381088    468148    579042 

Sums of square lengths of shortest basis vectors (must be the same for all flexible types): 
 dim:                5            10            20            30            40            50            60            70 

LLL5             2.55405e+12   6.60639e+13   5.76704e+14   1.08841e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL8             2.53724e+12   6.25167e+13   4.48781e+14   1.04244e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99            2.53662e+12   6.24462e+13   4.32408e+14    1.0006e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999         2.53662e+12   6.24462e+13   4.31436e+14    9.9619e+14   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999-pnew     2.5372e+12   6.24642e+13   4.33817e+14   1.00437e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL5-dual             152578         16137          9018          7967          7404          7083          6801          6593 
LLL8-dual             151770         15306          7032          6251          5884          5675          5522          5398 
LLL99-dual            151760         15270          6777          5923          5563          5344          5197          5088 
LLL99999-dual         151760         15270          6774          5889          5534          5319          5169          5057 
LLL99999-dnew         151782         15288          6831          6003          5658          5485          5359          5250 

****************************************************
TestBasisConstructSpeedLLL with m = 1048573
Types: Int = NTL::ZZ, Real = xdouble
Number of replications (different multipliers a): 1000
Total time: 85.1714 seconds.

Timings for different methods, in basic clock units (microseconds) 
 dim:                5        10        20        30        40        50        60        70  

LLL5              32188    122451    344985    586839    874499   1187942   1572402   2016671 
LLL8               4552     31191    308347    605627    904469   1243672   1632139   2039366 
LLL99              3379     18747    244520    761848   1205552   1658379   2125730   2697181 
LLL99999           2623     10191     52245    153153    251007    361952    513181    643524 
LLL99999-pnew     31401    184979   1204434   2581144   3858051   5153143   6597478   8172043 
UppTri             3867      7167     19241     34474     48156     62701     78325     95105 
mDualUT             805      2262     10697     30956     69802    133587    222567    340151 
LLL5-dual         25299    100870    243461    389741    576792    824942   1149441   1580528 
LLL8-dual          4172     31455    282937    459577    654218    904533   1225483   1632345 
LLL99-dual         3156     19309    244198    717224    945192   1219812   1567409   2002167 
LLL99999-dual      2535     10572     53248    145198    247999    383906    563633    802974 
LLL99999-dnew     25243    138145    779684   1465785   1792603   2148920   2584122   3121656 

Sums of square lengths of shortest basis vectors (must be the same for all flexible types): 
 dim:                5            10            20            30            40            50            60            70 

LLL5             2.55405e+12   6.60639e+13   5.76925e+14   1.08841e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL8             2.53809e+12   6.26165e+13   4.50306e+14   1.04391e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99             2.5372e+12   6.24909e+13   4.34266e+14   1.00553e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999          2.5372e+12   6.24868e+13   4.34138e+14   1.00444e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999-pnew     2.5372e+12   6.24642e+13   4.33817e+14   1.00437e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL5-dual             152578         16137          9018          7967          7404          7083          6801          6593 
LLL8-dual             151783         15340          7073          6348          6015          5837          5678          5569 
LLL99-dual            151790         15282          6813          5988          5676          5506          5380          5256 
LLL99999-dual         151783         15282          6806          5972          5664          5494          5361          5233 
LLL99999-dnew         151782         15288          6831          6003          5658          5485          5359          5250 

****************************************************
TestBasisConstructSpeedLLL with m = 1048573
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 1000
Total time: 110.667 seconds.

Timings for different methods, in basic clock units (microseconds) 
 dim:                5        10        20        30        40        50        60        70  

LLL5              24021     99753    353932    733665   1272238   1947090   2825425   3934405 
LLL8               4905     28268    272662    638990   1121399   1751636   2537746   3469002 
LLL99              4027     18726    215261    723172   1344939   2106629   3002978   4172116 
LLL99999           3487     11784     58729    172717    331356    556326    881242   1256544 
LLL99999-pnew     25019    163175   1152807   2800464   4808667   7181116   9996667  13347248 
UppTri             3640      6988     18836     34055     48051     62355     78139     95455 
mDualUT             778      2119     10445     30901     70314    129456    215341    333272 
LLL5-dual         18526     72804    190779    347491    582633    927576   1409880   2074678 
LLL8-dual          4596     28584    253577    419702    646764    980261   1445389   2078333 
LLL99-dual         3868     19375    224406    654416    897354   1236698   1708458   2346297 
LLL99999-dual      3400     12133     60247    168214    316121    533159    836174   1254197 
LLL99999-dnew     19733    108865    661443   1257541   1586231   2003046   2555593   3289728 

Sums of square lengths of shortest basis vectors (must be the same for all flexible types): 
 dim:                5            10            20            30            40            50            60            70 

LLL5             2.55405e+12   6.60639e+13   5.76704e+14   1.08841e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL8             2.53809e+12   6.26165e+13   4.50262e+14   1.04388e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99             2.5372e+12   6.24909e+13   4.34169e+14   1.00555e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999          2.5372e+12   6.24868e+13   4.34066e+14    1.0044e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999-pnew     2.5372e+12   6.24642e+13   4.33817e+14   1.00437e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL5-dual             152578         16137          9021          7967          7404          7083          6801          6593 
LLL8-dual             151783         15344          7071          6354          6019          5841          5679          5570 
LLL99-dual            151790         15282          6813          5988          5675          5507          5379          5257 
LLL99999-dual         151783         15282          6806          5972          5663          5494          5361          5235 
LLL99999-dnew         151782         15288          6831          6002          5657          5484          5358          5250 

****************************************************
TestBasisConstructSpeedLLL with m = 1048573
Types: Int = NTL::ZZ, Real = NTL::RR
Number of replications (different multipliers a): 1000
Total time: 805.781 seconds.

Timings for different methods, in basic clock units (microseconds) 
 dim:                5        10        20        30        40        50        60        70  

LLL5             139811    759833   2654773   4837652   7569294  10743124  14720528  19404128 
LLL8              17415    177119   2438013   5579544   8622462  12237204  16416752  20824134 
LLL99             13000    109802   1984440   7060797  11493073  16189262  21251850  27260266 
LLL99999           8673     55815    552004   2214690   2591103   3872637   5593570   7339529 
LLL99999-pnew    156786   1340414  10808684  24381604  37307309  50891678  66103442  82409240 
UppTri             4536      8004     20640     36857     51237     67050     84939    105112 
mDualUT             896      2389     11525     33082     73692    139064    232886    359231 
LLL5-dual         98331    482273   1340375   2481711   4284268   6899907  10485575  15221220 
LLL8-dual         15325    161057   1890888   3490875   5715262   8964776  13575096  19572592 
LLL99-dual        11342    100413   1629559   5406322   8138413  11863457  17017163  23916277 
LLL99999-dual      8497     49686    369187   1452131   3189665   5692628   9138820  14128263 
LLL99999-dnew    109924    782968   5384373  10600710  13491160  17275745  22268279  28689385 

Sums of square lengths of shortest basis vectors (must be the same for all flexible types): 
 dim:                5            10            20            30            40            50            60            70 

LLL5             2.55405e+12   6.60639e+13   5.76925e+14   1.08841e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL8             2.53724e+12   6.25167e+13   4.48786e+14   1.04223e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99            2.53662e+12   6.24462e+13   4.32412e+14   1.00055e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999         2.53662e+12   6.24462e+13   4.31437e+14   9.96135e+14   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL99999-pnew     2.5372e+12   6.24642e+13   4.33817e+14   1.00437e+15   1.09951e+15   1.09951e+15   1.09951e+15   1.09951e+15 
LLL5-dual             152578         16137          9018          7967          7404          7082          6801          6593 
LLL8-dual             151770         15306          7025          6254          5888          5682          5527          5404 
LLL99-dual            151760         15270          6777          5928          5562          5349          5201          5090 
LLL99999-dual         151760         15270          6773          5896          5535          5323          5173          5059 
LLL99999-dnew         151782         15288          6831          6004          5659          5484          5359          5250 
