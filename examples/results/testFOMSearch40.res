Types: Int = NTL::ZZ, Real = double 

=========================================================
Program TestFOMSearch with m = 1099511627791
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
Running Time in seconds: 8.39733

Best 3 multipliers `a` found, and their FOMs:
  971774789182, 0.233459
  196316975508, 0.21016
  896091630760, 0.190414

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 6.6111

Best 3 multipliers `a` found, and their FOMs:
  205016098870, 0.244236
  610039476260, 0.24007
  971774789182, 0.233459

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 6.70235

Best 3 multipliers `a` found, and their FOMs:
  205016098870, 0.244236
  610039476260, 0.24007
  971774789182, 0.233459

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.86856

Best 3 multipliers `a` found, and their FOMs:
  205016098870, 0.244236
  610039476260, 0.24007
  971774789182, 0.233459

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 7.992

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.042963

Best 3 multipliers `a` found, and their FOMs:
  205016098870, 0.244236
  610039476260, 0.24007
  971774789182, 0.233459

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 5.73719

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.067785

Best 3 multipliers `a` found, and their FOMs:
  205016098870, 0.244236
  610039476260, 0.24007
  971774789182, 0.233459

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 3.54856

Best 3 multipliers `a` found, and their FOMs:
  205016098870, 0.244236
  610039476260, 0.24007
  971774789182, 0.233459

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.61829

Best 3 multipliers `a` found, and their FOMs:
  205016098870, 0.244236
  610039476260, 0.24007
  971774789182, 0.233459

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.61734

Best 3 multipliers `a` found, and their FOMs:
  205016098870, 0.244236
  610039476260, 0.24007
  971774789182, 0.233459

==================================
FOM experiments in dual lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 6.27511

Best 3 multipliers `a` found, and their FOMs:
  847599096686, 0.196615
  211064359948, 0.196244
  470575592282, 0.195213

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 5.60592

Best 3 multipliers `a` found, and their FOMs:
  728999740898, 0.244919
  327380726433, 0.237775
  365585266176, 0.230093

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 5.67751

Best 3 multipliers `a` found, and their FOMs:
  728999740898, 0.244919
  327380726433, 0.237775
  365585266176, 0.230093

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 4.02501

Best 3 multipliers `a` found, and their FOMs:
  728999740898, 0.244919
  327380726433, 0.237775
  365585266176, 0.230093

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 6.27509

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.040279

Best 3 multipliers `a` found, and their FOMs:
  728999740898, 0.244919
  327380726433, 0.237775
  365585266176, 0.230093

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 50
Running Time in seconds: 4.84765

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 50
Number retained: 3
Running Time in seconds: 0.047297

Best 3 multipliers `a` found, and their FOMs:
  728999740898, 0.244919
  327380726433, 0.237775
  365585266176, 0.230093

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 3.5028

Best 3 multipliers `a` found, and their FOMs:
  728999740898, 0.244919
  327380726433, 0.237775
  365585266176, 0.230093

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.56962

Best 3 multipliers `a` found, and their FOMs:
  728999740898, 0.244919
  327380726433, 0.237775
  365585266176, 0.230093

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.57058

Best 3 multipliers `a` found, and their FOMs:
  728999740898, 0.244919
  327380726433, 0.237775
  365585266176, 0.230093

***     DONE     ***

=========================================================
Program TestFOMSearch with m = 1099511627791
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
Running Time in seconds: 568.274

Best 3 multipliers `a` found, and their FOMs:
  1078308853597, 0.242476
  1013555717205, 0.224149
  80085906043, 0.217217

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 33.4676

Best 3 multipliers `a` found, and their FOMs:
  149613909046, 0.312387
  739601306150, 0.30702
  404526173392, 0.306972

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 33.4804

Best 3 multipliers `a` found, and their FOMs:
  149613909046, 0.312387
  739601306150, 0.30702
  404526173392, 0.306972

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 3.56035

Best 3 multipliers `a` found, and their FOMs:
  149613909046, 0.312387
  739601306150, 0.30702
  404526173392, 0.306972

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 3.86501

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 1.27607

Best 3 multipliers `a` found, and their FOMs:
  149613909046, 0.312387
  739601306150, 0.30702
  404526173392, 0.306972

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 3.05233

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 1.2732

Best 3 multipliers `a` found, and their FOMs:
  149613909046, 0.312387
  739601306150, 0.30702
  404526173392, 0.306972

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.023

Best 3 multipliers `a` found, and their FOMs:
  149613909046, 0.312387
  739601306150, 0.30702
  404526173392, 0.306972

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.38377

Best 3 multipliers `a` found, and their FOMs:
  149613909046, 0.312387
  739601306150, 0.30702
  404526173392, 0.306972

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.39896

Best 3 multipliers `a` found, and their FOMs:
  149613909046, 0.312387
  739601306150, 0.30702
  404526173392, 0.306972

==================================
FOM experiments in dual lattices 

--------------------------------------
1. BKZ + BB  Discard: false,  from List: false
Number of multipliers a examined: 1000
Number retained: 3
Running Time in seconds: 542.612

Best 3 multipliers `a` found, and their FOMs:
  857901435329, 0.234897
  135541956334, 0.230048
  1044039895847, 0.223538

--------------------------------------
2. BKZ + BB.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 28.4009

Best 3 multipliers `a` found, and their FOMs:
  574505091240, 0.303645
  827087580236, 0.285447
  458614882408, 0.284415

--------------------------------------
2'. BKZ + BB, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 28.7714

Best 3 multipliers `a` found, and their FOMs:
  574505091240, 0.303645
  827087580236, 0.285447
  458614882408, 0.284415

--------------------------------------
3. LLL only.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.8574

Best 3 multipliers `a` found, and their FOMs:
  574505091240, 0.303645
  464380231974, 0.292343
  827087580236, 0.285447

--------------------------------------
4a. LLL only, stage 1, with vector t.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 3.14929

--------------------------------------
4b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 1.42027

Best 3 multipliers `a` found, and their FOMs:
  574505091240, 0.303645
  827087580236, 0.285447
  458614882408, 0.284415

--------------------------------------
5a. LLL only, stage 1 with vector t0.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 20
Running Time in seconds: 2.70801

--------------------------------------
5b. LLL + BB, stage 2.  Discard: true,  from List: true
Number of multipliers a examined: 20
Number retained: 3
Running Time in seconds: 1.44248

Best 3 multipliers `a` found, and their FOMs:
  574505091240, 0.303645
  827087580236, 0.285447
  458614882408, 0.284415

--------------------------------------
6. BKZ + BB, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 2.07103

Best 3 multipliers `a` found, and their FOMs:
  574505091240, 0.303645
  827087580236, 0.285447
  458614882408, 0.284415

--------------------------------------
7. LLL only, only non succ.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.39113

Best 3 multipliers `a` found, and their FOMs:
  574505091240, 0.303645
  464380231974, 0.292343
  827087580236, 0.285447

--------------------------------------
7'. LLL only, only non succ, with proj = lat.  Discard: true,  from List: false
Number of multipliers a examined: 100000
Number retained: 3
Running Time in seconds: 1.39721

Best 3 multipliers `a` found, and their FOMs:
  574505091240, 0.303645
  464380231974, 0.292343
  827087580236, 0.285447

***     DONE     ***
