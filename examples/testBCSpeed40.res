****************************************************
Types: Int = NTL::ZZ, Real = double

TestBasisConstructionSpeed with m = 1099511627791
Number of replications (different multipliers a): 1000

Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               9242     25049     83750    333426    655256 
LLL8               1537      3069     11462    158289    522590 
LLL99              1134      2213      7469     87437    526064 
LLL99999            943      1819      4962     29718    131241 
LLL99999-new       7793     26854    121571   1049870   3069098 
UppTri             2796      3735      6817     18961     36757 
mDualUT             623       990      2275     10395     28771 
LLL5-dual          6393     16605     49661    161330    248582 
LLL8-dual          1232      2713     10717    135600    370773 
LLL99-dual          954      2011      6900     78008    437591 
LLL99999-dual       833      1688      4685     27606    108689 
LLL99999-new       5591     17122     63642    433944   1204842 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5               6.1394e+20    7.34139e+22    4.39089e+24     2.1623e+26    9.03316e+26 
LLL8              6.12959e+20     7.2748e+22    4.18166e+24    1.18649e+26    5.14269e+26 
LLL99             6.12959e+20    7.27262e+22    4.17755e+24    1.13805e+26    4.28598e+26 
LLL99999          6.12959e+20    7.27262e+22    4.17755e+24    1.13522e+26    4.21957e+26 
LLL99999-new      6.12965e+20    7.27284e+22    4.18009e+24    1.13838e+26    4.28475e+26 
LLL5-dual         5.57425e+08    6.44999e+06         240824          46535          39233 
LLL8-dual         5.55584e+08    6.38118e+06         225870          26055          17516 
LLL99-dual        5.55584e+08    6.38118e+06         225626          24969          14845 
LLL99999-dual     5.55584e+08    6.38118e+06         225626          24938          14649 
LLL99999-new      5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 10.4896 seconds


Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               9300     25082     84704    340482    667309 
LLL8               1556      3072     11732    162852    534903 
LLL99              1100      2235      7605     90512    540318 
LLL99999            939      1871      5018     30960    134453 
LLL99999-new       7857     27377    124002   1073802   3139609 
UppTri             2784      3799      6913     19350     37604 
mDualUT             533       922      2219     11151     31635 
LLL5-dual          6245     16574     50006    164702    258590 
LLL8-dual          1221      2733     10765    137720    377636 
LLL99-dual          939      2061      7036     79677    444837 
LLL99999-dual       869      1703      4753     28309    110231 
LLL99999-new       5551     17042     65486    443721   1225500 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5               6.1394e+20    7.34139e+22    4.39089e+24     2.1623e+26    9.03316e+26 
LLL8              6.12959e+20     7.2748e+22    4.18166e+24    1.18649e+26    5.14269e+26 
LLL99             6.12959e+20    7.27262e+22    4.17755e+24    1.13805e+26    4.28598e+26 
LLL99999          6.12959e+20    7.27262e+22    4.17755e+24    1.13522e+26    4.21957e+26 
LLL99999-new      6.12965e+20    7.27284e+22    4.18009e+24    1.13838e+26    4.28475e+26 
LLL5-dual         5.57425e+08    6.44999e+06         240824          46535          39233 
LLL8-dual         5.55584e+08    6.38118e+06         225870          26055          17516 
LLL99-dual        5.55584e+08    6.38118e+06         225626          24969          14845 
LLL99999-dual     5.55584e+08    6.38118e+06         225626          24938          14649 
LLL99999-new      5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 10.679 seconds


****************************************************
Types: Int = NTL::ZZ, Real = xdouble

TestBasisConstructionSpeed with m = 1099511627791
Number of replications (different multipliers a): 1000

Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              33835    102149    354772   1305824   2304106 
LLL8               3365      8597     47661    709543   2281163 
LLL99              2551      5998     28259    375708   2375961 
LLL99999           2067      4440     15186     81905    273384 
LLL99999-new      30379    115080    544679   4489999  13009422 
UppTri             3113      4149      7426     19903     38450 
mDualUT             667      1110      2481     10908     30446 
LLL5-dual         20893     56889    173647    537152    750440 
LLL8-dual          2677      6555     32888    475230   1287649 
LLL99-dual         1974      4781     19187    258843   1591640 
LLL99999-dual      1705      3627     10693     57555    190023 
LLL99999-new      18738     59852    226988   1517039   4207809 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5               6.1394e+20    7.34139e+22    4.39089e+24     2.1623e+26    9.03316e+26 
LLL8              6.13062e+20    7.27787e+22    4.18791e+24    1.18477e+26    5.21214e+26 
LLL99             6.12965e+20    7.27284e+22    4.17884e+24    1.14105e+26    4.32757e+26 
LLL99999          6.12965e+20    7.27284e+22    4.17884e+24    1.13984e+26    4.31391e+26 
LLL99999-new      6.12965e+20    7.27284e+22    4.18009e+24    1.13838e+26    4.28475e+26 
LLL5-dual         5.57425e+08    6.44999e+06         240824          46535          39233 
LLL8-dual         5.55666e+08    6.38448e+06         226197          26169          17695 
LLL99-dual        5.55584e+08    6.38118e+06         225731          25096          14864 
LLL99999-dual     5.55584e+08    6.38118e+06         225686          25064          14821 
LLL99999-new      5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 40.2523 seconds


Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              32811    101627    354434   1302461   2283144 
LLL8               3298      8543     47637    708516   2285356 
LLL99              2474      5942     28165    375637   2386331 
LLL99999           2030      4446     15313     81382    274210 
LLL99999-new      30200    115066    545167   4500519  13046056 
UppTri             3135      4083      7415     20059     39158 
mDualUT             589       997      2392     11335     32695 
LLL5-dual         20571     56821    173369    536383    754506 
LLL8-dual          2631      6588     32845    475300   1283801 
LLL99-dual         2014      4761     19317    258580   1588910 
LLL99999-dual      1727      3631     10773     57136    190129 
LLL99999-new      18521     59737    226963   1520957   4214457 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5               6.1394e+20    7.34139e+22    4.39089e+24     2.1623e+26    9.03316e+26 
LLL8              6.13062e+20    7.27787e+22    4.18791e+24    1.18477e+26    5.21214e+26 
LLL99             6.12965e+20    7.27284e+22    4.17884e+24    1.14105e+26    4.32757e+26 
LLL99999          6.12965e+20    7.27284e+22    4.17884e+24    1.13984e+26    4.31391e+26 
LLL99999-new      6.12965e+20    7.27284e+22    4.18009e+24    1.13838e+26    4.28475e+26 
LLL5-dual         5.57425e+08    6.44999e+06         240824          46535          39233 
LLL8-dual         5.55666e+08    6.38448e+06         226197          26169          17695 
LLL99-dual        5.55584e+08    6.38118e+06         225731          25096          14864 
LLL99999-dual     5.55584e+08    6.38118e+06         225686          25064          14821 
LLL99999-new      5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 40.2558 seconds


****************************************************
Types: Int = NTL::ZZ, Real = quad_float

TestBasisConstructionSpeed with m = 1099511627791
Number of replications (different multipliers a): 1000

Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              18647     55328    212969   1001117   2176849 
LLL8               3015      6613     29564    496787   1744297 
LLL99              2393      5030     19096    271610   1698811 
LLL99999           2114      4111     12153     71768    240637 
LLL99999-new      16969     63690    340967   3271993  10035382 
UppTri             2664      3676      6762     18849     36907 
mDualUT             588       963      2261     10392     28997 
LLL5-dual         12062     31896     95080    305354    471185 
LLL8-dual          2530      5325     22108    287490    776633 
LLL99-dual         2094      4259     14580    165770    938547 
LLL99999-dual      1927      3611      9727     49377    151708 
LLL99999-new      11596     34607    129720    905303   2504572 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5               6.1394e+20    7.34139e+22    4.39089e+24    2.16224e+26    9.03316e+26 
LLL8              6.13062e+20    7.27787e+22    4.18791e+24    1.18481e+26    5.21214e+26 
LLL99             6.12965e+20    7.27284e+22    4.17884e+24    1.14105e+26    4.32719e+26 
LLL99999          6.12965e+20    7.27284e+22    4.17884e+24    1.13984e+26     4.3136e+26 
LLL99999-new      6.12965e+20    7.27284e+22    4.18009e+24    1.13838e+26    4.28475e+26 
LLL5-dual         5.57425e+08    6.44999e+06         240824          46520          39206 
LLL8-dual         5.55666e+08    6.38448e+06         226197          26160          17696 
LLL99-dual        5.55584e+08    6.38118e+06         225731          25092          14865 
LLL99999-dual     5.55584e+08    6.38118e+06         225686          25060          14819 
LLL99999-new      5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 28.9578 seconds


Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              18012     55245    213729   1005639   2190400 
LLL8               2962      6623     29810    498788   1757744 
LLL99              2403      5023     19334    273091   1712202 
LLL99999           2122      4121     12224     72395    239629 
LLL99999-new      17100     64087    344338   3300195  10139938 
UppTri             2623      3609      6769     18987     37176 
mDualUT             481       827      2130     10843     31402 
LLL5-dual         11953     31778     95530    306797    481274 
LLL8-dual          2510      5387     22307    288834    780821 
LLL99-dual         2106      4285     14845    166834    944381 
LLL99999-dual      1937      3627      9779     49646    152990 
LLL99999-new      11604     34626    130943    912109   2522452 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5               6.1394e+20    7.34139e+22    4.39089e+24    2.16224e+26    9.03316e+26 
LLL8              6.13062e+20    7.27787e+22    4.18791e+24    1.18481e+26    5.21214e+26 
LLL99             6.12965e+20    7.27284e+22    4.17884e+24    1.14105e+26    4.32719e+26 
LLL99999          6.12965e+20    7.27284e+22    4.17884e+24    1.13984e+26     4.3136e+26 
LLL99999-new      6.12965e+20    7.27284e+22    4.18009e+24    1.13838e+26    4.28475e+26 
LLL5-dual         5.57425e+08    6.44999e+06         240824          46520          39206 
LLL8-dual         5.55666e+08    6.38448e+06         226197          26160          17696 
LLL99-dual        5.55584e+08    6.38118e+06         225731          25092          14865 
LLL99999-dual     5.55584e+08    6.38118e+06         225686          25060          14819 
LLL99999-new      5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 29.1605 seconds


****************************************************
Types: Int = NTL::ZZ, Real = NTL::RR

TestBasisConstructionSpeed with m = 1099511627791
Number of replications (different multipliers a): 1000

Results for `testLoopResize` (many objects are created or resized)
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5             110418    403963   1691115   7758197  14274431 
LLL8               9181     31627    192920   3977701  14530412 
LLL99              6487     21265    117498   2177951  15134094 
LLL99999           5255     14886     61463    591345   3401129 
LLL99999-new     111602    492588   2828804  28148862  86366665 
UppTri             3428      4408      7849     21029     40380 
mDualUT             637      1172      2576     11068     30733 
LLL5-dual         78325    251328    899792   3227223   4856891 
LLL8-dual          8857     28751    181151   3394959   9759539 
LLL99-dual         6487     20044    109398   1886727  12326181 
LLL99999-dual      5209     14222     56839    477751   2430073 
LLL99999-new      80574    294324   1372974  11415558  32533503 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5               6.1394e+20    7.34139e+22    4.39089e+24     2.1623e+26    9.03316e+26 
LLL8              6.12959e+20     7.2748e+22    4.18166e+24    1.18649e+26    5.14269e+26 
LLL99             6.12959e+20    7.27262e+22    4.17755e+24    1.13802e+26    4.28547e+26 
LLL99999          6.12959e+20    7.27262e+22    4.17755e+24    1.13518e+26    4.21878e+26 
LLL99999-new      6.12965e+20    7.27284e+22    4.18009e+24    1.13838e+26    4.28475e+26 
LLL5-dual         5.57425e+08    6.44999e+06         240824          46535          39233 
LLL8-dual         5.55584e+08    6.38118e+06         225870          26055          17530 
LLL99-dual        5.55584e+08    6.38118e+06         225626          24969          14844 
LLL99999-dual     5.55584e+08    6.38118e+06         225626          24938          14648 
LLL99999-new      5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 268.427 seconds


Results for `testLoop No Resize`
Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5             110212    401277   1683971   7748011  14247710 
LLL8               9119     31552    193107   3970570  14493895 
LLL99              6406     21238    117610   2172465  15086233 
LLL99999           5270     14924     61648    589518   3386648 
LLL99999-new     111522    491010   2827169  28085728  86065910 
UppTri             3392      4364      7812     20852     40469 
mDualUT             665      1078      2518     11518     32981 
LLL5-dual         77889    249400    895945   3223564   4855877 
LLL8-dual          8812     28679    181238   3388159   9732400 
LLL99-dual         6416     20003    109727   1884215  12295617 
LLL99999-dual      5180     14221     56800    475666   2421495 
LLL99999-new      79652    292893   1369947  11392025  32415901 

Sums of square lengths of shortest basis vector (must be the same across all implementations):
 dim:                4              6             10             20             30  

LLL5               6.1394e+20    7.34139e+22    4.39089e+24     2.1623e+26    9.03316e+26 
LLL8              6.12959e+20     7.2748e+22    4.18166e+24    1.18649e+26    5.14269e+26 
LLL99             6.12959e+20    7.27262e+22    4.17755e+24    1.13802e+26    4.28547e+26 
LLL99999          6.12959e+20    7.27262e+22    4.17755e+24    1.13518e+26    4.21878e+26 
LLL99999-new      6.12965e+20    7.27284e+22    4.18009e+24    1.13838e+26    4.28475e+26 
LLL5-dual         5.57425e+08    6.44999e+06         240824          46535          39233 
LLL8-dual         5.55584e+08    6.38118e+06         225870          26055          17530 
LLL99-dual        5.55584e+08    6.38118e+06         225626          24969          14844 
LLL99999-dual     5.55584e+08    6.38118e+06         225626          24938          14648 
LLL99999-new      5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 267.621 seconds


