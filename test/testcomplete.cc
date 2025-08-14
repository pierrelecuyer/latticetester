#include <cstdlib>
#include <iostream>

int main() {
    std::cout << "##########################\nTesting Coordinates.h \n########################## \n \n";
    int result = system("./testCoordinates");
    std::cout << "##########################\nTesting Util.h \n########################## \n \n";
    result = system("./testUtil");
    std::cout << "##########################\nTesting Rank1Lattice.h \n########################## \n \n";
    result = system("./testRank1Lattice");
    return result;
}
