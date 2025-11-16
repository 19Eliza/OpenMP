#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

// =============================================================
//  КОНЦЕНТРИЧЕСКИЕ ЭЛЛИПСОИДЫ
// =============================================================
struct Ellipsoid {
    double a, b, c;  // полуоси
};

// Проверка точки внутри эллипсоида
bool insideEllipsoid(double x, double y, double z, const Ellipsoid& E)
{
    double nx = x / E.a;
    double ny = y / E.b;
    double nz = z / E.c;
    return (nx*nx + ny*ny + nz*nz) <= 1.0;
}

// Определение слоя точки
// 0 = brain(E1), 1 = skull(E2), 2 = skin(E3), 3 = outside
int determineLayer(double x, double y, double z,
                   const Ellipsoid& E1,
                   const Ellipsoid& E2,
                   const Ellipsoid& E3)
{
    if (insideEllipsoid(x,y,z,E1)) return 0;
    if (insideEllipsoid(x,y,z,E2)) return 1;
    if (insideEllipsoid(x,y,z,E3)) return 2;
    return 3;
}

// =============================================================
//   MONTE-CARLO
// =============================================================
void stratifiedMonteCarloVoxel(
    double x0, double x1,
    double y0, double y1,
    double z0, double z1,
    const Ellipsoid& E1,
    const Ellipsoid& E2,
    const Ellipsoid& E3,
    int N,
    double fraction[4])
{
    int count[4] = {0,0,0,0};

    double dx = (x1 - x0) / N;
    double dy = (y1 - y0) / N;
    double dz = (z1 - z0) / N;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
            {
                double x = x0 + (i + 0.5)*dx;
                double y = y0 + (j + 0.5)*dy;
                double z = z0 + (k + 0.5)*dz;

                int L = determineLayer(x, y, z, E1, E2, E3);
                count[L]++;
            }

    int total = N*N*N;
    for (int l = 0; l < 4; l++)
        fraction[l] = double(count[l]) / total;
}


int main()
{

    Ellipsoid E1{30, 40, 50};   // brain
    Ellipsoid E2{35, 45, 55};   // skull
    Ellipsoid E3{40, 50, 60};   // skin  (внешний)

    // ----------------------------------------------
    // 
    // Минимальный осевой параллелепипед
    // ----------------------------------------------
    double xmin = -E3.a;
    double xmax =  E3.a;

    double ymin = -E3.b;
    double ymax =  E3.b;

    double zmin = -E3.c;
    double zmax =  E3.c;

    // ----------------------------------------------
    // Строим 3D-сетку
    // ----------------------------------------------
    const int Nx = 40;   // число вокселей по x
    const int Ny = 40;   // число вокселей по y
    const int Nz = 40;   // число вокселей по z

    double dx = (xmax - xmin) / Nx;
    double dy = (ymax - ymin) / Ny;
    double dz = (zmax - zmin) / Nz;

    const int Nsamples = 5;

    // ----------------------------------------------
    //  Перебираем воксели
    // ----------------------------------------------
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
            {
                double x0 = xmin + i*dx;
                double x1 = x0 + dx;

                double y0 = ymin + j*dy;
                double y1 = y0 + dy;

                double z0 = zmin + k*dz;
                double z1 = z0 + dz;

                double frac[4];
                stratifiedMonteCarloVoxel(x0,x1, y0,y1, z0,z1,
                                          E1, E2, E3,
                                          Nsamples,
                                          frac);


                if (frac[0] > 0 || frac[1] > 0 || frac[2] > 0)
                {
                    cout << "Voxel (" << i << "," << j << "," << k << "):\n";
                    cout << "   brain="  << frac[0] << "\n";
                    cout << "   skull="  << frac[1] << "\n";
                    cout << "   skin="   << frac[2] << "\n";
                    cout << "   outside="<< frac[3] << "\n";
                }
            }

    return 0;
}
