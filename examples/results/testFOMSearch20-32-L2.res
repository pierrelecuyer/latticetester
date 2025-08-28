Types: Int = NTL::ZZ, Real = double 

=========================================================
Program TestFOMSearch with m = 1048573
Norm type: L2
Vector t = [32 64 16 12 10]
Total number of projections with t:  486
Vector t0 = [6 64 16 12]
Total number of projections with t0: 335

==================================
FOM experiments in primal lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 20.8637

Best 3 multipliers `a` found, and their FOMs:
  189965, 0.216723
  635286, 0.207607
  643649, 0.205

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 5.0259

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 5.04756

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.96946

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 5.10456

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.057718

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 3.4319

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.096437

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.7159

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.7832

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.78457

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

==================================
FOM experiments in dual lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 19.6412

Best 3 multipliers `a` found, and their FOMs:
  759294, 0.226986
  62787, 0.219858
  567536, 0.208627

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.59049

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.63467

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.55781

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 4.12492

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.047982

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 3.1665

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.063546

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.45402

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.61494

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.6134

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

***     DONE     ***
