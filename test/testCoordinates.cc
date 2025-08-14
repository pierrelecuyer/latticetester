#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <iostream>
#include <cstdint>
#include <cmath>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <latticetester/Coordinates.h>
#include <latticetester/CoordinateSets.h>

namespace LatticeTester {

TEST_CASE("Test different types of coordinates") {
    // Default constructor: empty set
    Coordinates c;
    Coordinates c_test;

    // Construct from vector
    std::vector<std::size_t> vec = {3, 1, 4, 1, 5, 7, 9, 12};
    std::vector<std::size_t> vec_test = {1, 3, 4, 5, 7, 9, 12};
    c = Coordinates(vec);
    c_test = Coordinates(vec_test);
    CHECK( c == c_test);

    // Construct from iterator range (only first 4 elements)
    c = Coordinates (vec.begin(), vec.begin() + 4);
    vec_test = {1, 3, 4};
    c_test = Coordinates(vec_test);
    CHECK( c == c_test);

    // Add missing coordinates in c and add extra coordinate '1' which may not have an influence.
    c.insert(5); c.insert(7), c.insert(9); c.insert(12); c.insert(1);
    vec_test = {1, 3, 4, 5, 7, 9, 12};
    c_test = Coordinates(vec_test);
    CHECK( c == c_test);
}

};
