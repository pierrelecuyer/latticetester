TestBasisConstructionSmall 
Types: Int = NTL::ZZ, Real = double

Initial Korobov lattice basis (triangular) = 
[[1 12 144 707 316]
[0 1021 0 0 0]
[0 0 1021 0 0]
[0 0 0 1021 0]
[0 0 0 0 1021]
]
Square length of first basis vector: 620586

Basis after LLL with delta=0.5: 
[[-78 85 -1 -12 -144]
[13 156 -170 2 24]
[-7 -84 13 156 -170]
[-107 -263 -93 -95 -119]
[-332 100 179 106 251]
]
Square length of first basis vector: 34190

Basis after LLL with delta=0.99999: 
[[-78 85 -1 -12 -144]
[13 156 -170 2 24]
[-7 -84 13 156 -170]
[-94 -107 -263 -93 -95]
[-248 87 23 276 249]
]
Square length of first basis vector: 34190

After `lowerTriangularBasis`: 
[[1021 0 0 0 0]
[0 1021 0 0 0]
[0 0 1021 0 0]
[0 0 0 1021 0]
[-42 -504 78 -85 1]
]

After `upperTriangularBasis`: 
[[1 12 144 -314 316]
[0 1021 0 0 0]
[0 0 1021 0 0]
[0 0 0 1021 0]
[0 0 0 0 1021]
]

m-dual of upper-triangular basis: 
[[1021 0 0 0 0]
[-12 1 0 0 0]
[-144 0 1 0 0]
[314 0 0 1 0]
[-316 0 0 0 1]
]

m-dual basis after LLL with delta=0.99999: 
[[-2 0 0 1 1]
[2 2 2 1 0]
[2 -1 -1 4 -2]
[-1 -1 5 -1 0]
[-2 2 1 2 -5]
]
Square length of first dual basis vector: 6

Lattice projection over coordinates {1,3,5}.
In the following basisProj matrices, we need 5 rows and 3 columns
 to make the projection, then 3 rows and 3 columns for the basis.
 When part of matrix is not used, it must be ignored.

basisProj after projectMatrix (the generating vectors): 
[[1 144 316]
[0 0 0]
[0 1021 0]
[0 0 0]
[0 0 1021]
]
Basis for this projection, obtained with LLL: 
[[71 14 -26]
[-42 78 1]
[-78 -1 -144]
[0 0 0]
[0 0 0]
]
Upper-triangular basis for this proj. (first 3 rows):
[[1 144 316]
[0 1021 0]
[0 0 1021]
[0 0 0]
[0 0 0]
]
Triangular basis for this projection, with `buildProjection`: 
[[1 144 316]
[0 1021 0]
[0 0 1021]
]
Triangular basis for m-dual of this projection: 
[[1021 0 0]
[-144 1 0]
[-316 0 1]
]
m-dual basis of proj after LLL with delta=0.99999: 
[[-2 -1 -6]
[-2 12 1]
[11 6 -6]
]
Square length of first m-dual basis vector: 41

Triangular basis for m-dual of this projection, with `buildProjectionDual`: 
[[1021 0 0]
[144 1 0]
[316 0 1]
]
We now look at the direct projection of the dual over the coordinates in proj.
Generating vectors for the projection of the dual: 
[[-2 0 1]
[2 2 0]
[2 -1 -2]
[-1 5 0]
[-2 1 -5]
]
Reduced basis for this projection (first 3 rows), after LLL with delta=0.99999: 
[[0 1 0]
[0 0 1]
[-1 0 0]
[0 0 0]
[0 0 0]
]
Square length of first m-dual basis vector: 1

We see that the dual of the projection differs from the projection of the dual! 

