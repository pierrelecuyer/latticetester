Types: Int = int64_t, Real = double
TestBasisConstructionSmall 

Initial Korobov lattice basis (triangular) = 
[[1 33 79 82 80]
[0 101 0 0 0]
[0 0 101 0 0]
[0 0 0 101 0]
[0 0 0 0 101]
]
Square length of first basis vector: 20455

Basis after LLL with delta=0.5: 
[[9 -6 4 31 13]
[5 -37 -9 6 -4]
[-1 -33 22 19 21]
[21 -14 43 5 -37]
[-28 -15 10 27 -18]
]
Square length of first basis vector: 1263

Basis after LLL with delta=0.99999: 
[[9 -6 4 31 13]
[-15 10 27 -18 12]
[-10 -27 18 -12 8]
[-28 -15 10 27 -18]
[21 -14 43 5 -37]
]
Square length of first basis vector: 1263

After `upperTriangularBasis`: 
[[1 -68 79 82 80]
[0 101 0 0 0]
[0 0 101 0 0]
[0 0 0 101 0]
[0 0 0 0 101]
]

m-dual of upper-triangular basis: 
[[101 0 0 0 0]
[68 1 0 0 0]
[-79 0 1 0 0]
[-82 0 0 1 0]
[-80 0 0 0 1]
]

m-dual basis after LLL with delta=0.99999: 
[[-1 0 -1 0 1]
[-2 0 0 1 -1]
[0 3 0 1 -1]
[1 -2 -2 0 -1]
[-1 0 -1 -2 -2]
]
Square length of first dual basis vector: 3

Lattice projection over coordinates {1,3,5}.
In the following basisProj matrices, we need 5 rows and 3 columns
 to make the projection, then 3 rows and 3 columns for the basis.
 When part of matrix is not used, it must be ignored.

basisProj after projectMatrix (the generating vectors): 
[[1 79 80]
[0 0 0]
[0 101 0]
[0 0 0]
[0 0 101]
]
Basis for this projection, obtained with LLL: 
[[5 -9 -4]
[9 4 13]
[35 38 -28]
[0 0 0]
[0 0 0]
]
Upper-triangular basis for this proj. (first 3 rows):
[[1 79 80]
[0 101 0]
[0 0 101]
[0 0 0]
[0 0 0]
]
Triangular basis for this projection, with `buildProjection`: 
[[1 79 80]
[0 101 0]
[0 0 101]
]
Triangular basis for m-dual of this projection: 
[[101 0 0]
[-79 1 0]
[-80 0 1]
]
m-dual basis of proj after LLL with delta=0.99999: 
[[-1 -1 1]
[-4 0 -5]
[6 -7 -2]
]
Square length of first m-dual basis vector: 3

Triangular basis for m-dual of this projection, with `buildProjectionDual`: 
[[101 0 0]
[79 1 0]
[80 0 1]
]
We now look at the direct projection of the dual over the coordinates in proj.
Generating vectors for the projection of the dual: 
[[-1 -1 1]
[-2 0 -1]
[0 0 -1]
[1 -2 -1]
[-1 -1 -2]
]
Reduced basis for this projection (first 3 rows), after LLL with delta=0.99999: 
[[0 0 -1]
[0 -1 0]
[-1 0 0]
[0 0 0]
[0 0 0]
]
Square length of first m-dual basis vector: 1

We see that the dual of the projection differs from the projection of the dual! 

