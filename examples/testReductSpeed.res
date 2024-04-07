****************************************************
Types: Int = long, Real = double

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 20 replications 

 dim:               5        10        20        30        40  

LLL5               328       889      2494      5956     11265 
LLL99999           202      1116      6842     20124     30113 
BKZ99999           249      1165     10122     59642    102257 
Direct BB            0         0         0         0         0 
Pairwise+BB          0         0         0         0         0 
LLL5+BB            198       799      7351    331226  11539412 
LLL8+BB            149       958      5942    171672  10772047 
LLL99999+BB        144      1096      8311     93186   4940680 
BKZ99999+BB        177       954     11293     93667   3000573 
LLL8+BKZ+BB        208       989     10394    102810   2447330 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5              2835        349        181        161        152 
LLL99999          2835        300        134        117        112 
BKZ99999          2835        300        132        115        109 
Direct BB            0          0          0          0          0 
Pairwise+BB          0          0          0          0          0 
LLL5+BB           2835        300        132        113        107 
LLL8+BB           2835        300        132        113        107 
LLL99999+BB       2835        300        132        113        107 
BKZ99999+BB       2835        300        132        113        107 
LLL8+BKZ+BB       2835        300        132        113        107 

Total time for everything: 33.8004 seconds

Interestingly, we see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = double

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 20 replications 

 dim:               5        10        20        30        40  

LLL5               424       830      1700      2916      4789 
LLL99999           253      1246      5852     10809     13797 
BKZ99999           265      1354      9634     32820     73071 
Direct BB            0         0         0         0         0 
Pairwise+BB          0         0         0         0         0 
LLL5+BB            381      1359     10601    396208  14325057 
LLL8+BB            349      1579      9204    196078  13155528 
LLL99999+BB        328      1586      9744    104032   6216338 
BKZ99999+BB        352      1752     13080     86246   3768803 
LLL8+BKZ+BB        429      2005     12824     90816   3093337 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5              2835        349        181        161        152 
LLL99999          2835        300        134        117        112 
BKZ99999          2835        300        132        115        109 
Direct BB            0          0          0          0          0 
Pairwise+BB          0          0          0          0          0 
LLL5+BB           2835        300        132        113        107 
LLL8+BB           2835        300        132        113        107 
LLL99999+BB       2835        300        132        113        107 
BKZ99999+BB       2835        300        132        113        107 
LLL8+BKZ+BB       2835        300        132        113        107 

Total time for everything: 41.6874 seconds

Interestingly, we see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = quad_float

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 20 replications 

 dim:               5        10        20        30        40  

LLL5               921      3151      8979     13240     22667 
LLL99999           842      5038     33118     52391     66043 
BKZ99999           965      5993     60346    203422    394843 
Direct BB            0         0         0         0         0 
Pairwise+BB          0         0         0         0         0 
LLL5+BB           1096      4475     50750   3115973 115271371 
LLL8+BB           1073      4885     30126   1559655 106355228 
LLL99999+BB       1183      6400     37413    745912  50261477 
BKZ99999+BB       1178      8036     60705    609267  30167919 
LLL8+BKZ+BB       1385      8258     66584    655558  28576601 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5              2835        349        181        161        152 
LLL99999          2835        300        134        117        112 
BKZ99999          2835        300        132        115        109 
Direct BB            0          0          0          0          0 
Pairwise+BB          0          0          0          0          0 
LLL5+BB           2835        300        132        113        107 
LLL8+BB           2835        300        132        113        107 
LLL99999+BB       2835        300        132        113        107 
BKZ99999+BB       2835        300        132        113        107 
LLL8+BKZ+BB       2835        300        132        113        107 

Total time for everything: 338.514 seconds

Interestingly, we see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = xdouble

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 20 replications 

 dim:               5        10        20        30        40  

LLL5               983      3122      7900     17949     19279 
LLL99999           766      4451     27440     59535     61564 
BKZ99999           841      4882     42726    169356    308713 
Direct BB            0         0         0         0         0 
Pairwise+BB          0         0         0         0         0 
LLL5+BB            936      3994     57878   2780370 103813282 
LLL8+BB            981      4652     30355   1438497  97345179 
LLL99999+BB        975      5316     41611    709211  45217713 
BKZ99999+BB        947      6327     54590    524777  27480393 
LLL8+BKZ+BB       1121      6683     68081    543364  25553419 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5              2835        349        181        161        152 
LLL99999          2835        300        134        117        112 
BKZ99999          2835        300        132        115        109 
Direct BB            0          0          0          0          0 
Pairwise+BB          0          0          0          0          0 
LLL5+BB           2835        300        132        113        107 
LLL8+BB           2835        300        132        113        107 
LLL99999+BB       2835        300        132        113        107 
BKZ99999+BB       2835        300        132        113        107 
LLL8+BKZ+BB       2835        300        132        113        107 

Total time for everything: 306.457 seconds

Interestingly, we see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = NTL::RR

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 20 replications 

 dim:               5        10        20        30        40  

LLL5              3618     16154     46691     91777    159168 
LLL99999          3862     28189    219009    456358    593116 
BKZ99999          4258     33485    442301   1695383   3549224 
Direct BB            0         0         0         0         0 
Pairwise+BB          0         0         0         0         0 
LLL5+BB           4408     19149    366608  22159803 850350472 
LLL8+BB           4406     25388    206910  11755216 828830926 
LLL99999+BB       4125     29192    287703   5608078 387201544 
BKZ99999+BB       5285     34959    496761   4676544 233667805 
LLL8+BKZ+BB       5295     38014    441864   4846998 158366862 

Sums of square lengths of shortest basis vector:
 dim:              5         10         20         30         40 

LLL5              2835        349        181        161        152 
LLL99999          2835        300        134        117        112 
BKZ99999          2835        300        132        115        111 
Direct BB            0          0          0          0          0 
Pairwise+BB          0          0          0          0          0 
LLL5+BB           2835        300        132        113        107 
LLL8+BB           2835        300        132        113        107 
LLL99999+BB       2835        300        132        113        107 
BKZ99999+BB       2835        300        132        113        107 
LLL8+BKZ+BB       2835        300        132        113        107 

Total time for everything: 2516.82 seconds

Interestingly, we see that LLL or BKZ alone do not always find a shortest vector.


