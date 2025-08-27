Types: Int = NTL::ZZ, Real = double 

=========================================================
Program TestFOMSearch with m = 1099511627791
Norm type: L2
Total number of projections with t:  446
Total number of projections with t0: 305

==================================
FOM experiments in primal lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 8.30528

Best 3 multipliers `a` found, and their FOMs:
  315661629133, 0.241259
  971774789182, 0.233459
  413866867400, 0.215464

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 7.4937

Best 3 multipliers `a` found, and their FOMs:
  918927574849, 0.264833
  1089166524987, 0.259433
  242005322085, 0.255688

--------------------------------------
2 bis. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 7.45615

Best 3 multipliers `a` found, and their FOMs:
  918927574849, 0.264833
  1089166524987, 0.259433
  242005322085, 0.255688

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 5.44752

Best 3 multipliers `a` found, and their FOMs:
  918927574849, 0.264833
  1089166524987, 0.259433
  242005322085, 0.255688

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 9.31555

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.038833

Best 3 multipliers `a` found, and their FOMs:
  918927574849, 0.264833
  1089166524987, 0.259433
  242005322085, 0.255688

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 6.5132

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.048484

Best 3 multipliers `a` found, and their FOMs:
  918927574849, 0.264833
  1089166524987, 0.259433
  242005322085, 0.255688

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.356

Best 3 multipliers `a` found, and their FOMs:
  918927574849, 0.264833
  1089166524987, 0.259433
  242005322085, 0.255688

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 3.17375

Best 3 multipliers `a` found, and their FOMs:
  918927574849, 0.264833
  1089166524987, 0.259433
  242005322085, 0.255688

--------------------------------------
8. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 3.17727

Best 3 multipliers `a` found, and their FOMs:
  918927574849, 0.264833
  1089166524987, 0.259433
  242005322085, 0.255688

==================================
FOM experiments in dual lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 6.22166

Best 3 multipliers `a` found, and their FOMs:
  772303778413, 0.223998
  1005828199026, 0.211894
  220160789285, 0.208582

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 6.16408

Best 3 multipliers `a` found, and their FOMs:
  452005645564, 0.269388
  107648311535, 0.267524
  353522122454, 0.25459

--------------------------------------
2 bis. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 6.22542

Best 3 multipliers `a` found, and their FOMs:
  452005645564, 0.269388
  107648311535, 0.267524
  353522122454, 0.25459

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.36964

Best 3 multipliers `a` found, and their FOMs:
  452005645564, 0.269388
  107648311535, 0.267524
  353522122454, 0.25459

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 7.13311

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.030392

Best 3 multipliers `a` found, and their FOMs:
  452005645564, 0.269388
  107648311535, 0.267524
  353522122454, 0.25459

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 5.51915

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.040224

Best 3 multipliers `a` found, and their FOMs:
  452005645564, 0.269388
  107648311535, 0.267524
  353522122454, 0.25459

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.07891

Best 3 multipliers `a` found, and their FOMs:
  452005645564, 0.269388
  107648311535, 0.267524
  353522122454, 0.25459

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.95787

Best 3 multipliers `a` found, and their FOMs:
  452005645564, 0.269388
  107648311535, 0.267524
  353522122454, 0.25459

--------------------------------------
8. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.95782

Best 3 multipliers `a` found, and their FOMs:
  452005645564, 0.269388
  107648311535, 0.267524
  353522122454, 0.25459

***     DONE     ***

=========================================================
Program TestFOMSearch with m = 1099511627791
Norm type: L1
Total number of projections with t:  54
Total number of projections with t0: 50

==================================
FOM experiments in primal lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 561.281

Best 3 multipliers `a` found, and their FOMs:
  399450937967, 0.354632
  182815190781, 0.323092
  113816658047, 0.306782

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 28.576

Best 3 multipliers `a` found, and their FOMs:
  893774817428, 0.37854
  473774920134, 0.371557
  98649027740, 0.366513

--------------------------------------
2 bis. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 28.5285

Best 3 multipliers `a` found, and their FOMs:
  893774817428, 0.37854
  473774920134, 0.371557
  98649027740, 0.366513

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 3.0753

Best 3 multipliers `a` found, and their FOMs:
  893774817428, 0.37854
  473774920134, 0.371557
  98649027740, 0.366513

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 3.29596

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 1.80666

Best 3 multipliers `a` found, and their FOMs:
  893774817428, 0.37854
  473774920134, 0.371557
  98649027740, 0.366513

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 2.45869

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 2.3311

Best 3 multipliers `a` found, and their FOMs:
  893774817428, 0.37854
  473774920134, 0.371557
  98649027740, 0.366513

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.17898

Best 3 multipliers `a` found, and their FOMs:
  473774920134, 0.388023
  459735993627, 0.380207
  893774817428, 0.37854

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 0.788018

Best 3 multipliers `a` found, and their FOMs:
  473774920134, 0.392551
  459735993627, 0.380207
  893774817428, 0.37854

--------------------------------------
8. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 0.790523

Best 3 multipliers `a` found, and their FOMs:
  473774920134, 0.392551
  459735993627, 0.380207
  893774817428, 0.37854

==================================
FOM experiments in dual lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 545.729

Best 3 multipliers `a` found, and their FOMs:
  451314873491, 0.348062
  182815190781, 0.346219
  178647025976, 0.330939

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 25.118

Best 3 multipliers `a` found, and their FOMs:
  878116461874, 0.377746
  577626986586, 0.375141
  657702149524, 0.371948

--------------------------------------
2 bis. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 25.1051

Best 3 multipliers `a` found, and their FOMs:
  878116461874, 0.377746
  577626986586, 0.375141
  657702149524, 0.371948

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.25826

Best 3 multipliers `a` found, and their FOMs:
  878116461874, 0.377746
  577626986586, 0.375141
  657702149524, 0.371948

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 2.41778

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 0.579561

Best 3 multipliers `a` found, and their FOMs:
  878116461874, 0.377746
  577626986586, 0.375141
  657702149524, 0.371948

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 1.93518

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 0.579675

Best 3 multipliers `a` found, and their FOMs:
  878116461874, 0.377746
  577626986586, 0.375141
  657702149524, 0.371948

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.15366

Best 3 multipliers `a` found, and their FOMs:
  878116461874, 0.377746
  577626986586, 0.375141
  657702149524, 0.374809

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 0.763446

Best 3 multipliers `a` found, and their FOMs:
  878116461874, 0.377746
  577626986586, 0.375141
  657702149524, 0.374809

--------------------------------------
8. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 0.764039

Best 3 multipliers `a` found, and their FOMs:
  878116461874, 0.377746
  577626986586, 0.375141
  657702149524, 0.374809

***     DONE     ***
