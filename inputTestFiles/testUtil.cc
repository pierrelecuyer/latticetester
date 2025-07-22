#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <iostream>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "../include/latticetester/Util.h"

NTL::RR r, s;
double d, e;

TEST_CASE("testing the swap function") {
    r = 39150490341.1029501;
    s = -13481480.190341;
    NTL::RR r_copy, s_copy;
    r_copy = r;
    s_copy = s;
    LatticeTester::swap9(r, s);
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
    CHECK(LatticeTester::mysqrt(d) == sqrt(d));
    CHECK(LatticeTester::mysqrt(-d) == -1);
    d = 27.9019512;
    // Should return the logarithm in base 2
    CHECK(LatticeTester::Lg(d) == log(d)/log(2));
}

TEST_CASE("testing the round function") {
    CHECK(LatticeTester::Round(1.3) == 1);
    CHECK(LatticeTester::Round(-1.5) == -2);
    CHECK(LatticeTester::Round(491409134.4914804814) == 491409134);
    CHECK(LatticeTester::Round(-491409134.4914804814) == -491409134);   
}

TEST_CASE("testing the factorial function") {
    CHECK(LatticeTester::Factorial(0) == 1);
    CHECK(LatticeTester::Factorial(1) == 1);
    CHECK(LatticeTester::Factorial(2) == 2);
    CHECK(LatticeTester::Factorial(3) == 6);
    CHECK(LatticeTester::Factorial(10) == 3628800); 
}

TEST_CASE("testing arithmetic functions") {
    NTL::ZZ a, b, q;
    a = 97;
    b = 3;
    LatticeTester::Quotient(a, b, q);
    CHECK(q == a/b);
    b = 7;
    LatticeTester::Modulo(a, b, q);
    CHECK(q == a % b);
    a = -818;
    LatticeTester::Modulo(a, b, q);
    CHECK(q == a % b);
    LatticeTester::ModuloPos(a, b, q);
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
    LatticeTester::ModuloVec(v, b);
    CHECK(v == w);
    a = 9;
    b = 7;
    LatticeTester::ModuloTowardZero(a, b, q);
    CHECK(q  == 2);
    a = 13;
    b = 7;
    LatticeTester::ModuloTowardZero(a, b, q);
    CHECK(q  == -1);
    r = 39150490341.1029501;
    s = -13481480.190341;
    NTL::RR qu, rem;
    LatticeTester::Divide(qu, rem, r, s);
    CHECK(qu == trunc(r/s));
    CHECK(rem == r - trunc(r/s)*s);
    LatticeTester::DivideRound(qu, r, s);
    CHECK(qu == -2904);
}

TEST_CASE("testing Euclidean algoritm functions") {
    CHECK(LatticeTester::gcd(9 ,7) == 1);
    CHECK(LatticeTester::gcd(4158, 15561) == 63);
    NTL::ZZ a, b, c, d, e, f, g;
    a = 4158; b = 15561;
    LatticeTester::Euclide(a, b, c, d, e, f, g);
    CHECK(c*a + d*b == 63);
    CHECK(e*a + f*b == 0);
}
