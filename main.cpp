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

using tripleEllipsoid = std::array<Ellipsoid, 3>;

// Тестовые случаи
struct TestCase{
    std::vector<tripleEllipsoid> testEllipsoid;// множество тестовых эллипсоидов
    int countPointX = 32;
    int countPointY = 32;
    int countPointZ = 32;
    std::unordered_set<int> Nsamples{2048}; // количество точек разбросанных в одном элементе

};

// Проверка нахождения точки внутри эллипсоида
bool checkInsideEllipsoid(double x, double y, double z, const Ellipsoid& E)
{
    double nx = x / E.a;
    double ny = y / E.b;
    double nz = z / E.c;
    return (nx*nx + ny*ny + nz*nz) <= 1.0;
}

// Определение слоя для точки
// 0 = brain(E1), 1 = skull(E2), 2 = skin(E3), 3 = outside
int determineLayer(double x, double y, double z,
                   const Ellipsoid& E1,
                   const Ellipsoid& E2,
                   const Ellipsoid& E3)
{
    if (checkInsideEllipsoid(x,y,z,E1)) return 0;
    if (checkInsideEllipsoid(x,y,z,E2)) return 1;
    if (checkInsideEllipsoid(x,y,z,E3)) return 2;
    return 3;
}

// случайный коэффициент от 0 до 1
inline double randomCoeff()
{
    return double(rand()) / RAND_MAX; 
}

// Генерация случайной точки внутри элемента
void randomPointInElement(
    double x0, double x1,
    double y0, double y1,
    double z0, double z1,
    double &x, double &y, double &z)
{
    x = x0 + randomCoeff() * (x1 - x0);
    y = y0 + randomCoeff() * (y1 - y0);
    z = z0 + randomCoeff() * (z1 - z0);
}

// =============================================================
//   MONTE-CARLO
// =============================================================
void MonteCarloMethod(
    double x0, double x1,
    double y0, double y1,
    double z0, double z1,
    const Ellipsoid& E1,
    const Ellipsoid& E2,
    const Ellipsoid& E3,
    int Nsamples,
    double fraction[4])
{
    int count[4] = {0,0,0,0};

    for (int t = 0; t < Nsamples; t++)
    {
        double x, y, z;
        randomPointInElement(x0, x1, y0, y1, z0, z1, x, y, z);

        int L = determineLayer(x, y, z, E1, E2, E3);
        count[L]++;
    }

    for (int i = 0; i < 4; i++)
        fraction[i] = double(count[i]) / Nsamples;
}


void writeHeader(ofstream& fout)
{
    fout << left
         << setw(5)  << "i"
         << setw(5)  << "j"
         << setw(5)  << "k"
         << setw(12) << "brain"
         << setw(12) << "skull"
         << setw(12) << "skin"
         << setw(12) << "outside"
         << "\n";

    fout << fixed << setprecision(6); 
}

// Функция форматированного вывода строки
void writeVoxelLine(ofstream& fout,
                    int i, int j, int k,
                    const double frac[4])
{
    fout << left
         << setw(5)  << i
         << setw(5)  << j
         << setw(5)  << k
         << setw(12) << frac[0]
         << setw(12) << frac[1]
         << setw(12) << frac[2]
         << setw(12) << frac[3]
         << "\n";
}


int main()
{

    Ellipsoid E1{30, 40, 50};   // brain
    Ellipsoid E2{35, 45, 55};   // skull
    Ellipsoid E3{40, 50, 60};   // skin 

    // ----------------------------------------------
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
    const int Nx = 32;   // число точек по x
    const int Ny = 32;   // число точек по y
    const int Nz = 32;   // число точек по z

    double dx = (xmax - xmin) / Nx;
    double dy = (ymax - ymin) / Ny;
    double dz = (zmax - zmin) / Nz;

    const int Nsamples = 2048;

    ofstream fout("result.txt");
    writeHeader(fout);
    fout.precision(8);

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
                MonteCarloMethod(x0,x1, y0,y1, z0,z1,
                                          E1, E2, E3,
                                          Nsamples,
                                          frac);


                // Сохраняем только если точка попала в какой нибудь слой
                if (frac[0] > 0 || frac[1] > 0 || frac[2] > 0)
                {
                    writeVoxelLine(fout, i, j, k, frac);
                }
            }

    fout.close();
    cout << "Date is saved to voxels.txt\n";

    return 0;
}
