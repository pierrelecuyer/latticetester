Types: Int = NTL::ZZ, Real = double 

=========================================================
Program TestFOMSearch with m = 1048573
Norm type: L2
Total number of projections with t:  454
Total number of projections with t0: 305

==================================
FOM experiments in primal lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 20.6956

Best 3 multipliers `a` found, and their FOMs:
  301710, 0.237731
  879684, 0.229267
  157833, 0.22367

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 5.14877

Best 3 multipliers `a` found, and their FOMs:
  861454, 0.270944
  837499, 0.267057
  315393, 0.262839

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 5.14128

Best 3 multipliers `a` found, and their FOMs:
  861454, 0.270944
  837499, 0.267057
  315393, 0.262839

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 3.00835

Best 3 multipliers `a` found, and their FOMs:
  861454, 0.270944
  837499, 0.267057
  315393, 0.262839

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 5.42121

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.048656

Best 3 multipliers `a` found, and their FOMs:
  861454, 0.270944
  837499, 0.267057
  315393, 0.262839

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 3.70492

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.094876

Best 3 multipliers `a` found, and their FOMs:
  861454, 0.270944
  837499, 0.267057
  315393, 0.262839

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.82358

Best 3 multipliers `a` found, and their FOMs:
  861454, 0.270944
  837499, 0.267057
  315393, 0.262839

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.81528

Best 3 multipliers `a` found, and their FOMs:
  861454, 0.270944
  837499, 0.267057
  315393, 0.262839

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.81549

Best 3 multipliers `a` found, and their FOMs:
  861454, 0.270944
  837499, 0.267057
  315393, 0.262839

==================================
FOM experiments in dual lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 19.464

Best 3 multipliers `a` found, and their FOMs:
  301710, 0.250676
  473811, 0.237418
  158478, 0.229681

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.90161

Best 3 multipliers `a` found, and their FOMs:
  464079, 0.265395
  489009, 0.260139
  773657, 0.258507

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.90407

Best 3 multipliers `a` found, and their FOMs:
  464079, 0.265395
  489009, 0.260139
  773657, 0.258507

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.68571

Best 3 multipliers `a` found, and their FOMs:
  464079, 0.265395
  489009, 0.260139
  773657, 0.258507

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 4.35935

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.044076

Best 3 multipliers `a` found, and their FOMs:
  464079, 0.265395
  489009, 0.260139
  773657, 0.258507

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 3.35838

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.065192

Best 3 multipliers `a` found, and their FOMs:
  464079, 0.265395
  489009, 0.260139
  773657, 0.258507

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.69777

Best 3 multipliers `a` found, and their FOMs:
  464079, 0.265395
  489009, 0.260139
  773657, 0.258507

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.7229

Best 3 multipliers `a` found, and their FOMs:
  464079, 0.265395
  489009, 0.260139
  773657, 0.258507

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.72462

Best 3 multipliers `a` found, and their FOMs:
  464079, 0.265395
  489009, 0.260139
  773657, 0.258507

***     DONE     ***
