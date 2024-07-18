Types: NTL::ZZ, double 

=============================================================
Results for a small MRG example with m=13, k=3, a=(7,0,4). 
===================================================
We build the lattice and look at projections over successive coordinates.

Figure of merit primal succ, with BB.
coordinates      merit        sqlen   minmerit 
{1,2,3,4}       0.44285        1  0.44285
{1,2,3,4,5}     0.770305        7  0.44285
{1,2,3,4,5,6}   0.712803        11  0.44285
{1,2,3,4,5,6,7} 0.664498        15  0.44285
{1,...,8}       0.711607        25  0.44285
FOM value: 0.44285

===================================================
Then we look at other primal projections, over pairs and triples.

Figure of merit primal non-succ, with BB.
coordinates      merit        sqlen    minmerit  
{1,2}           1        1  1
{1,3}           1        1  1
{1,4}           1        1  1
{1,5}           1        1  1
{1,2,3}         1        1  1
{1,2,4}         1        1  1
{1,2,5}         1        1  1
{1,3,4}         2.23607        5  1
{1,3,5}         1        1  1
{1,4,5}         1        1  1
FOM value: 1

===================================================
We now examine the m-dual lattice, first for successive coordinates.

Figure of merit for dual.
coordinates      merit        sqlen   minmerit  
{1,2,3,4}       0.66143        29  0.66143
{1,2,3,4,5}     0.853946        24  0.66143
{1,2,3,4,5,6}   0.607881        8  0.607881
{1,2,3,4,5,6,7} 0.700048        8  0.607881
{1,...,8}       0.764366        8  0.607881
FOM value: 0.607881

===================================================
Then we look at the m-duals of the other projections.

Figure of merit dual non-succ, with BB.
coordinates      merit        sqlen   minmerit  
{1,2}           1        169  1
{1,3}           1        169  1
{1,4}           1        169  1
{1,5}           1        169  1
{1,2,3}         1        169  1
{1,2,4}         1        169  1
{1,2,5}         1        169  1
{1,3,4}         0.414243        29  0.414243
{1,3,5}         1        169  0.414243
{1,4,5}         1        169  0.414243
FOM value: 0.414243

===================================================
A closer look at the projection over coordinates {1,3,4}: 
Full basis B before taking projection {1,3,4}: 
[[1 0 0 4 2 1 10 13]
[0 1 0 0 4 2 1 10]
[0 0 1 7 10 9 13 1]
[0 0 0 13 0 0 0 0]
[0 0 0 0 13 0 0 0]
[0 0 0 0 0 13 0 0]
[0 0 0 0 0 0 13 0]
[0 0 0 0 0 0 0 13]
]
Basis for projection {1,3,4}: 
[[1 0 4]
[0 1 7]
[0 0 13]
]
Basis for m-dual projection {1,3,4}: 
[[13 0 0]
[0 13 0]
[-4 -7 1]
]

=============================================================
Results for a MRG example with m = 9223372036854773561, a = (1145902849652723, 0, -1184153554609676). 
===================================================
We build the lattice and look at projections over successive coordinates.

Figure of merit primal succ, with BB.
coordinates      merit        sqlen   minmerit 
{1,2,3,4}       1.52588e-05        1  1.52588e-05
{1,2,3,4,5}     0.000223483        1.12465e+08  1.52588e-05
{1,2,3,4,5,6}   0.00817336        1.02613e+15  1.52588e-05
{1,2,3,4,5,6,7} 0.358334        1.09841e+21  1.52588e-05
{1,...,8}       0.678151        4.67515e+23  1.52588e-05
FOM value: 1.52588e-05

===================================================
Then we look at other primal projections, over pairs and triples.

Figure of merit primal non-succ, with BB.
coordinates      merit        sqlen    minmerit  
{1,2}           1        1  1
{1,3}           1        1  1
{1,4}           1        1  1
{1,5}           1        1  1
{1,2,3}         1        1  1
{1,2,4}         1        1  1
{1,2,5}         1        1  1
{1,3,4}         10119.9        1.02412e+08  1
{1,3,5}         1        1  1
{1,4,5}         1        1  1
FOM value: 1

===================================================
We now examine the m-dual lattice, first for successive coordinates.

Figure of merit for dual.
coordinates      merit        sqlen   minmerit  
{1,2,3,4}       4.91481e-07        9.56894e+15  4.91481e-07
{1,2,3,4,5}     0.000332039        9.56894e+15  4.91481e-07
{1,2,3,4,5,6}   0.0249593        9.56894e+15  4.91481e-07
{1,2,3,4,5,6,7} 0.541513        9.56894e+15  4.91481e-07
{1,...,8}       0.51637        8.92522e+13  4.91481e-07
FOM value: 4.91481e-07

===================================================
Then we look at the m-duals of the other projections.

Figure of merit dual non-succ, with BB.
coordinates      merit        sqlen   minmerit  
{1,2}           1        8.50706e+37  1
{1,3}           1        8.50706e+37  1
{1,4}           1        8.50706e+37  1
{1,5}           1        8.50706e+37  1
{1,2,3}         1        8.50706e+37  1
{1,2,4}         1        8.50706e+37  1
{1,2,5}         1        8.50706e+37  1
{1,3,4}         1.06058e-11        9.56894e+15  1.06058e-11
{1,3,5}         1        8.50706e+37  1.06058e-11
{1,4,5}         1        8.50706e+37  1.06058e-11
FOM value: 1.06058e-11

===================================================
A closer look at the projection over coordinates {1,3,4}: 
Full basis B before taking projection {1,3,4}: 
[[1 0 0 9222187883300163885 6629950407641799081 7721957468801369126 18037989607192828976 3020738662083495906]
[0 1 0 0 9222187883300163885 6629950407641799081 7721957468801369126 18037989607192828976]
[0 0 1 1145902849652723 6821832294614679088 11315916339925500873 10351102864336227360 11198630209724120468]
[0 0 0 9223372036854773561 0 0 0 0]
[0 0 0 0 9223372036854773561 0 0 0]
[0 0 0 0 0 9223372036854773561 0 0]
[0 0 0 0 0 0 9223372036854773561 0]
[0 0 0 0 0 0 0 9223372036854773561]
]
Basis for projection {1,3,4}: 
[[1 0 9222187883300163885]
[0 1 1145902849652723]
[0 0 9223372036854773561]
]
Basis for m-dual projection {1,3,4}: 
[[9223372036854773561 0 0]
[0 9223372036854773561 0]
[-9222187883300163885 -1145902849652723 1]
]
