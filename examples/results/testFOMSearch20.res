Types: Int = NTL::ZZ, Real = double 

=========================================================
Program TestFOMSearch with m = 1048573
Norm type: L2
Vector t = [24 64 16 12 10]
Total number of projections with t:  478
Vector t0 = [6 64 16 12]
Total number of projections with t0: 335

==================================
FOM experiments in primal lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 6.15321

Best 3 multipliers `a` found, and their FOMs:
  189965, 0.216723
  635286, 0.207607
  643649, 0.205

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.53445

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.57292

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.87617

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 4.54155

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.024926

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 3.38191

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.034133

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.68084

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.74291

Best 3 multipliers `a` found, and their FOMs:
  850223, 0.256577
  871334, 0.253085
  762482, 0.244952

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.74166

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
Running Time in seconds: 5.14868

Best 3 multipliers `a` found, and their FOMs:
  759294, 0.226986
  62787, 0.219858
  567536, 0.208627

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.038

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.1135

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.5377

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 3.97956

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.020179

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 3.13993

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.024648

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.47524

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.60426

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.60396

Best 3 multipliers `a` found, and their FOMs:
  489009, 0.258975
  773657, 0.258507
  410357, 0.246321

***     DONE     ***

=========================================================
Program TestFOMSearch with m = 1048573
Norm type: L1
Vector t = [12 64 8 6 5]
Total number of projections with t:  102
Vector t0 = [6 64 8 6]
Total number of projections with t0: 96

==================================
FOM experiments in primal lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 770.233

Best 3 multipliers `a` found, and their FOMs:
  189381, 0.267971
  567536, 0.267515
  964624, 0.249977

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 371.58

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.297565
  838315, 0.282751
  263003, 0.280823

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 372.408

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.297565
  838315, 0.282751
  263003, 0.280823

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.18892

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.297565
  838315, 0.282751
  263003, 0.280823

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 2.3981

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 1.40042

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.297565
  838315, 0.282751
  263003, 0.280823

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 1.98771

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 1.40172

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.297565
  838315, 0.282751
  263003, 0.280823

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.63468

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.297565
  838315, 0.282751
  263003, 0.280823

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.03186

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.297565
  838315, 0.282751
  263003, 0.280823

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.0328

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.297565
  838315, 0.282751
  263003, 0.280823

==================================
FOM experiments in dual lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 304.032

Best 3 multipliers `a` found, and their FOMs:
  189381, 0.267971
  964624, 0.249977
  759294, 0.225558

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 16.3265

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.309098
  758184, 0.301387
  1017459, 0.300101

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 16.3355

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.309098
  758184, 0.301387
  1017459, 0.300101

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.90071

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.309098
  758184, 0.301387
  1017459, 0.300101

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 2.10579

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 0.187404

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.309098
  758184, 0.301387
  1017459, 0.300101

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 1.79863

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 0.187175

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.309098
  758184, 0.301387
  1017459, 0.300101

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.57012

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.309098
  758184, 0.301387
  1017459, 0.300101

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 0.975515

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.309098
  758184, 0.301387
  1017459, 0.300101

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 0.974909

Best 3 multipliers `a` found, and their FOMs:
  597282, 0.309098
  758184, 0.301387
  1017459, 0.300101

***     DONE     ***
