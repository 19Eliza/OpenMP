#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <fstream>
#include <iomanip>

using namespace std;

// =============================================================
//  КОНЦЕНТРИЧЕСКИЕ ЭЛЛИПСОИДЫ
// =============================================================
struct Ellipsoid {
    double a, b, c;  // полуоси
};

// Проверка нахождения точки внутри эллипсоида
bool checkInsideEllipsoid(double x, double y, double z, const Ellipsoid& E)
{
    double nx = x / E.a;
    double ny = y / E.b;
    double nz = z / E.c;
    return (nx*nx + ny*ny + nz*nz) <= 1.0;
}