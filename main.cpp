#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <fstream>
#include <iomanip>

#include"FormatFileOuter.h"
#include"Ellipsoid.h"
#include"TestCaseStruct.h"

#include <omp.h>

using namespace std;

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


int main()
{   
    for(int i=0;i<TestCase::testEllipsoid.size();++i){

        auto test = TestCase::testEllipsoid[i];

        Ellipsoid E1{test[0]}; // brain
        Ellipsoid E2{test[1]}; // skull
        Ellipsoid E3{test[2]}; // skin

        double xmin = -E3.a, xmax = E3.a;
        double ymin = -E3.b, ymax = E3.b;
        double zmin = -E3.c, zmax = E3.c;

        int Nx = TestCase::countPointX, Ny = TestCase::countPointY, Nz = TestCase::countPointZ;
        double dx = (xmax - xmin) / Nx;
        double dy = (ymax - ymin) / Ny;
        double dz = (zmax - zmin) / Nz;

        int Nsamples = TestCase::Nsamples;

        string fileName = "resultForTest"+to_string(i+1)+".txt";
        ofstream fout(fileName);
        writeHeader(fout);

        #pragma omp parallel
        {
            std::ostringstream local_buffer;  

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


                        if (frac[0] > 0 || frac[1] > 0 || frac[2] > 0)
                        {
                            #pragma omp critical
                            writeVoxelLine(fout,i,j,k,frac);
                        }
                    }

            #pragma omp critical
            {
                fout << local_buffer.str();
            }
        } 

        fout.close();
        cout << "Data is saved to" <<fileName<<endl;
    }
    return 0;
}