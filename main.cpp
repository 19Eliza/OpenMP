#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <fstream>
#include <iomanip>

#include <omp.h>

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
    Ellipsoid E1{30.0, 40.0, 50.0}; // brain
    Ellipsoid E2{30.5, 40.5, 50.5}; // skull
    Ellipsoid E3{31.0, 41.0, 51.0}; // skin

    double xmin = -E3.a, xmax = E3.a;
    double ymin = -E3.b, ymax = E3.b;
    double zmin = -E3.c, zmax = E3.c;

    const int Nx = 32, Ny = 32, Nz = 32;
    double dx = (xmax - xmin) / Nx;
    double dy = (ymax - ymin) / Ny;
    double dz = (zmax - zmin) / Nz;

    const int Nsamples = 2048;

    ofstream fout("result.txt");
    writeHeader(fout);

    // ----------------------------------------------
    // Параллельный блок
    // ----------------------------------------------
    #pragma omp parallel
    {
        std::ostringstream local_buffer;  // каждый поток пишет в свой буфер

        cout << "Number of threads = " << omp_get_num_threads() << endl;

        #pragma omp for collapse(3) schedule(dynamic)
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
                    MonteCarloMethod(x0,x1, y0,y1, z0,z1, E1, E2, E3, Nsamples, frac);

                    int tid = omp_get_thread_num();

                    #pragma omp critical
                    cout << "Element (" << i << "," << j << "," << k << ") processed by thread " << tid << endl;

                    if (frac[0] > 0 || frac[1] > 0 || frac[2] > 0)
                    {
                        local_buffer << left
                                     << setw(5)  << i
                                     << setw(5)  << j
                                     << setw(5)  << k
                                     << setw(12) << frac[0]
                                     << setw(12) << frac[1]
                                     << setw(12) << frac[2]
                                     << setw(12) << frac[3]
                                     << "\n";
                    }
                }

        // Синхронизация: поток безопасно пишет свой буфер в файл
        #pragma omp critical
        {
            fout << local_buffer.str();
        }
    } // конец параллельного блока

    fout.close();
    cout << "Data is saved to result.txt\n";
    return 0;
}