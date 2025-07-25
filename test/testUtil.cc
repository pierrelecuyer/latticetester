#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <iostream>
#include <cstdint>
#include <cmath>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "../include/latticetester/Util.h"

NTL::RR r, s;
double d, e;

namespace LatticeTester {

TEST_CASE("testing the swap function") {
    r = 39150490341.1029501;
    s = -13481480.190341;
    NTL::RR r_copy, s_copy;
    r_copy = r;
    s_copy = s;
    swap9(r, s);
    CHECK (r == s_copy);
    CHECK (s == r_copy);
}

TEST_CASE("testing the sign function and the absolute value") {
    r = 5.0764951903940321;
    CHECK(sign(r) == 1);
    CHECK(abs(r) == r);
    r = -9891243901.9410134;
    CHECK(sign(r) == -1);
    CHECK(abs(r) == -r);
}

TEST_CASE("testing root and logarithmic functions") {
    d = 5.0764951903940321;
    // Should return the root for positive input and -1 for negative input
    CHECK(mysqrt(d) == sqrt(d));
    CHECK(mysqrt(-d) == -1);
    d = 27.9019512;
    // Should return the logarithm in base 2
    CHECK(Lg(d) == log(d)/log(2));
}

TEST_CASE("testing the round function") {
    CHECK(Round(1.3) == 1);
    CHECK(Round(-1.5) == -1);
    CHECK(Round(491409134.4914804814) == 491409134);
    CHECK(Round(-491409134.4914804814) == -491409134);   
}

TEST_CASE("testing the factorial function") {
    CHECK(Factorial(0) == 1);
    CHECK(Factorial(1) == 1);
    CHECK(Factorial(2) == 2);
    CHECK(Factorial(3) == 6);
    CHECK(Factorial(10) == 3628800); 
}

TEST_CASE("testing arithmetic functions") {
    NTL::ZZ a, b, q;
    a = 97;
    b = 3;
    Quotient(a, b, q);
    CHECK(q == a/b);
    b = 7;
    Modulo(a, b, q);
    CHECK(q == a % b);
    a = -818;
    Modulo(a, b, q);
    CHECK(q == a % b);
    ModuloPos(a, b, q);
    CHECK(q == a % b);
    NTL::Vec<NTL::ZZ> v, w;
    v.SetLength(3);
    w.SetLength(3);
    v[0] = 27;
    v[1] = 139;
    v[2] = -123;
    w[0] = v[0] % b;
    w[1] = v[1] % b;
    w[2] = v[2] % b;
    ModuloVec(v, b);
    CHECK(v == w);
    a = 9;
    b = 7;
    ModuloTowardZero(a, b, q);
    CHECK(q  == 2);
    a = 13;
    b = 7;
    ModuloTowardZero(a, b, q);
    CHECK(q  == -1);
    r = 39150490341.1029501;
    s = -13481480.190341;
    NTL::RR qu, rem;
    Divide(qu, rem, r, s);
    CHECK(qu == trunc(r/s));
    CHECK(rem == r - trunc(r/s)*s);
    DivideRound(qu, r, s);
    CHECK(qu == -2904);
}

TEST_CASE("testing Euclidean algorithm functions") {
    CHECK(gcd(9 ,7) == 1);
    CHECK(gcd(4158, 15561) == 63);
    NTL::ZZ a, b, c, d, e, f, g;
    a = 4158; b = 15561;
    Euclide(a, b, c, d, e, f, g);
    CHECK(c*a + d*b == 63);
    CHECK(e*a + f*b == 0);    
    a = 4158; b = -15561;
    Euclide(a, b, c, d, e, f, g);
    CHECK(c*a + d*b == -63);
    CHECK(e*a + f*b == 0);    
}

TEST_CASE("testing matrix and vector functions") {
    NTL::Mat<NTL::ZZ> M, MT, N;
    NTL::Vec<NTL::RR> u, v, w;
    NTL::Vec<NTL::ZZ> e, f, g;
    NTL::RR *a;
    NTL::RR *b = new NTL::RR[5]; 
    NTL::RR c;
    NTL::ZZ m;
    m = 8;
    double d;
    *b = NTL::RR(-3.1412);
    c = -3.1412;
    M.SetDims(3, 3);
    MT.SetDims(3, 3);
    N.SetDims(3, 3);    
    M[0][0] = 1; M[0][1] = -2; M[0][2] = 3;
    M[1][0] = -4; M[1][1] = 5; M[1][2] = -6;
    M[2][0] = 7; M[2][1] = -8; M[2][2] = 9;   
    N[0][0] = 1; N[1][0] = -2; N[2][0] = 3;
    N[0][1] = -4; N[1][1] = 5; N[2][1] = -6;
    N[0][2] = 7; N[1][2] = -8; N[2][2] = 9;
    TransposeMatrix<NTL::ZZ>(M, MT); 
    CHECK(MT == N);
    CreateVect<NTL::Vec<NTL::RR>>(v, 4); // Code in Util.h seems to be errorneous
    w.SetLength(5); w[0] = 0; w[1] = 0; w[2] = 0; w[3] = 0; w[4] = 0;
    CHECK(v==w);
    v[0] = 1; v[1] = 2; v[2] = 3; v[3] = 4; v[4] = 5;
    w[0] = 0; w[1] = 0; w[2] = 0; w[3] = 4; w[4] = 5;
    SetZero<NTL::Vec<NTL::RR>>(v, 3);
    CHECK(v==w);
    SetValue<NTL::RR*>(&a, 3, b);
    CHECK(*a == c);
    v[0] = 1; v[1] = 2; v[2] = 3; v[3] = 4; v[4] = 5;
    u.SetLength(5); u[0] = -5; u[1] = -4; u[2] = -3; u[3] = -2; u[4] = -1;
    ProdScal(v, u, 5, d);
    CHECK(d == -35);
    Invert(v, w, 5);
    CHECK(w==u);
    CalcNorm(u, 5, d, L2NORM);
    CHECK(d == 55);
    CalcNorm(u, 5, d, L1NORM);
    CHECK(d == 15);
    CopyPartVec(v, u, 3);
    w[0] = v[0]; w[1]= v[1]; w[2] = v[2]; w[3]= u[3]; w[4] = u[4];
    CHECK(u == w);    
    CopyPartVec(v, u, 4);
    w[0] = v[0]; w[1]= v[1]; w[2] = v[2]; w[3]= v[3]; w[4] = u[4];
    CHECK(u == w); 
    CopyPartMat(M, MT, 2, 2);
    N[0][0] = M[0][0]; N[0][1] = M[0][1]; N[1][0] = M[1][0]; N[1][1] = M[1][1];
    N[0][2] = MT[0][2]; N[1][2] = MT[1][2]; N[2][0] = MT[2][0]; N[2][1] = MT[2][1]; N[2][2] = MT[2][2];
    CHECK(MT == N); 
    copy(M, MT, 2, 2); // This is a duplicate of CopyPartMat
    CHECK(MT == N);
    N = MT;
    std::cout << "N:" << N << "\n";
    std::cout << "MT:" << MT << "\n";
    copy(MT, M); // This is a duplicate of CopyPartMat
    CHECK(M == N);
    u[0] = v[0]+w[0]*5; u[1] = v[1]+w[1]*5; u[2] = v[2]+w[2]*5; u[3] = v[3]; u[4] = v[4];
    ModifVect(v, w, 5, 3);
    CHECK(v == u); 
    e.SetLength(5); f.SetLength(5); g.SetLength(5);
    e[0] = 3; e[1] = 7; e[2] = -9; e[3] = 12; e[4] = 19;
    f[0] = 1; f[1] = -5; f[2] = -3; f[3] = 21; f[4] = 4;
    g[0] = (e[0]+f[0]*7) % m;  g[1] = (e[1]+f[1]*7) % m; g[2] = (e[2]+f[2]*7) % m; g[3] = (e[3]+f[3]*7) % m; g[4] = e[4];
    ModifVectModulo(e, f, 7, m, 4);
    CHECK(e == g);
    g[0] = -e[0]; g[1] = -e[1]; g[2] = e[2]; g[3]=e[3]; g[4]=e[4];
    ChangeSign(e, 2);
    CHECK(e == g);
    std::vector<std::int64_t> V;
    V.resize(7); V[0] = 2; V[1] = 7; V[2] = 10; V[3] = 5; V[4] = -25; V[5] = 17; V[6] = 8;
    CHECK(GCD2vect(V, 3, 5)==5);
    NTL::RR **A = new NTL::RR*[3];
    CreateMatr(A, 3);
    CHECK(A[2][2]==0);
    CreateMatr(M, 4);
    CHECK(M[3][3]==0);
    CopyMatr(M, MT, 2);
    CreateMatr(N, 4); N[0][0] = MT[0][0]; N[0][1] = MT[0][1]; N[1][0]= MT[1][0]; N[1][1] = MT[1][1];
    CHECK(N == M);
    M[0][0] = 7; M[1][1] = 3; M[2][2] = -2; M[3][3] = -2;
    ProductDiagonal(M, 4, m); // This is indeeed an error in the code. The index needs to start with 0 not with 1
    CHECK(m==84);
    ProductDiagonal(M, 3, m); // This is indeeed an error in the code. The index needs to start with 0 not with 1
    CHECK(m==-42);
    CHECK(CheckTriangular(M, 4, m)==false);
    M[1][0] = -7*m;
    CHECK(CheckTriangular(M, 4, m)==true);
    CHECK(checkInverseModm(M, N, m)== false);
    CreateMatr(M, 4);
    CreateMatr(N, 4);
    m=13;
    M[0][0] = 13; M[0][1] = 5; M[0][2] = -8; M[0][3] = 2; M[1][1] = 1; M[2][2] = 1; M[3][3] = 1;
    N[0][0] = 1; N[1][0] = -5; N[2][0] = 8; N[3][0] = -2; N[1][1] = 13; N[2][2] = 13; N[3][3] = 13;
    CreateMatr(MT, 4);    
    CHECK(checkInverseModm(M, N, m)== true);
}

};