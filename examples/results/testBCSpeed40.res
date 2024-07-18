****************************************************
Types: Int = NTL::ZZ, Real = double

TestBasisConstructionSpeed with m = 1099511627791
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5               9091     24566     82749    328300    647623 
LLL8               1511      2998     11397    156730    517296 
LLL99              1116      2177      7360     87042    524321 
LLL99999            925      1790      4831     29675    131211 
LLL99999-new       7748     26824    121051   1039301   3067810 
UppTri             2750      3652      6641     18378     36596 
mDualUT             546       908      2180     10859     32350 
LLL5-dual          6138     16210     48542    158453    245504 
LLL8-dual          1221      2648     10351    133130    367624 
LLL99-dual          930      1992      6734     77102    436197 
LLL99999-dual       832      1648      4594     27493    107773 
LLL99999-dnew      5409     16729     63151    429619   1194126 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
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
LLL99999-dnew     5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 10.3851 seconds


****************************************************
Types: Int = NTL::ZZ, Real = xdouble

TestBasisConstructionSpeed with m = 1099511627791
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              30410     94460    329990   1216369   2131915 
LLL8               3082      8019     44609    666244   2138287 
LLL99              2280      5573     26416    354011   2238758 
LLL99999           1880      4178     14290     77524    256900 
LLL99999-new      28019    107329    509878   4223794  12234819 
UppTri             2782      3768      6679     18241     35637 
mDualUT             543       969      2213     10619     30864 
LLL5-dual         19193     52727    161235    500551    702831 
LLL8-dual          2499      6172     30696    443699   1206456 
LLL99-dual         1921      4474     18146    242668   1492085 
LLL99999-dual      1627      3477     10174     54439    179844 
LLL99999-dnew     17324     55768    211805   1418582   3945662 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
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
LLL99999-dnew     5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 37.7162 seconds


****************************************************
Types: Int = NTL::ZZ, Real = quad_float

TestBasisConstructionSpeed with m = 1099511627791
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5              17039     53273    206973    970207   2105888 
LLL8               2850      6309     28680    477732   1717423 
LLL99              2283      4868     18494    261753   1674710 
LLL99999           2050      4043     11710     69712    234713 
LLL99999-new      16512     62140    331701   3183983   9905085 
UppTri             2468      3498      6471     18240     36065 
mDualUT             438       831      2072     10496     30935 
LLL5-dual         11477     30716     92100    296887    462244 
LLL8-dual          2434      5229     21239    277801    752263 
LLL99-dual         2051      4164     14112    160233    915277 
LLL99999-dual      1892      3509      9334     48047    147595 
LLL99999-dnew     11133     33523    126093    875814   2438464 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
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
LLL99999-dnew     5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 28.2943 seconds


****************************************************
Types: Int = NTL::ZZ, Real = NTL::RR

TestBasisConstructionSpeed with m = 1099511627791
Number of replications (different multipliers a): 1000

Timings for different methods, in basic clock units (microseconds) 

 dim:           4         6        10        20        30  

LLL5             103977    379353   1606215   7347408  13459361 
LLL8               8631     29797    183393   3753702  13700725 
LLL99              6094     20158    111680   2053418  14219736 
LLL99999           4981     14190     58471    555611   3196322 
LLL99999-new     105081    464263   2686001  26542168  81056553 
UppTri             3160      4121      7254     19301     37502 
mDualUT             574      1023      2312     10890     31047 
LLL5-dual         73228    236681    850026   3050290   4567482 
LLL8-dual          8290     27339    171811   3215065   9214799 
LLL99-dual         6085     19103    103698   1786435  11654407 
LLL99999-dual      4943     13697     53830    451547   2293912 
LLL99999-dnew     75153    278579   1303697  10793264  30658544 

Sums of square lengths of shortest basis vector (must be the same for all flexible types):
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
LLL99999-dnew     5.55584e+08    6.38118e+06         225642          25152          14923 

Total time: 252.774 seconds


