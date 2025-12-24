#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <fstream>
#include <iomanip>
#include <random>

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

// Случайный коэффициент от 0 до 1
inline double randomCoeff(std::mt19937 &gen)
{   
    static thread_local std::uniform_real_distribution<double> dist(0.0,1.0);
    return dist(gen);
}

// Генерация случайной точки внутри элемента
void randomPointInElement(
    double x0, double x1,
    double y0, double y1,
    double z0, double z1,
    double &x, double &y, double &z)
{   
    std::mt19937 gen(1234 + omp_get_thread_num());
    x = x0 + randomCoeff(gen) * (x1 - x0);
    y = y0 + randomCoeff(gen) * (y1 - y0);
    z = z0 + randomCoeff(gen) * (z1 - z0);
}

// =============================================================
//   Метод Монте-Карло. Приближённо оценивает объём каждого эллипсоида в элементе. 
//   В методе внутри каждого элемента разбрасываются равномерно случайные точки. 
//   Далее определяет, к какому слою принадлежит каждая точка.
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

        /// Минимальный ограничивающий параллелепипед
        double xmin = -E3.a, xmax = E3.a;
        double ymin = -E3.b, ymax = E3.b;
        double zmin = -E3.c, zmax = E3.c;

        int Nx = TestCase::countPointX, Ny = TestCase::countPointY, Nz = TestCase::countPointZ; /// количество точек разбиения по каждой оси

        /// шаги сетки 
        double dx = (xmax - xmin) / Nx; 
        double dy = (ymax - ymin) / Ny;
        double dz = (zmax - zmin) / Nz;

        int Nsamples = TestCase::Nsamples;

        // Файл вывода результатов
        string fileName = "resultForTest"+to_string(i+1)+".txt";
        ofstream fout(fileName);
        writeHeader(fout);

        double t_start = omp_get_wtime();


        #pragma omp parallel
        {
            std::ostringstream local_buffer;  

            #pragma omp for collapse(3) schedule(dynamic)
            for (int i = 0; i < Nx; i++) // номер элемента вдоль оси X 
                for (int j = 0; j < Ny; j++) // номер элемента вдоль оси Y
                    for (int k = 0; k < Nz; k++) // номер элемента вдоль оси Z
                    {
                        double x0 = xmin + i*dx;
                        double x1 = x0 + dx;
                        double y0 = ymin + j*dy;
                        double y1 = y0 + dy;
                        double z0 = zmin + k*dz;
                        double z1 = z0 + dz;

                        double frac[4];
                        MonteCarloMethod(x0,x1, y0,y1, z0,z1, E1, E2, E3, Nsamples, frac);

                        Вывод в файл с результаттами, если хотя бы одна доля в элементе > 0
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
        double t_end = omp_get_wtime();
        double elapsed = t_end - t_start;

        fout.close();
        int threads = omp_get_max_threads();
        cout << "Time for " << threads << " threads: " << elapsed << " seconds" << endl;
        cout << "Data is saved to" <<fileName<<endl;

    }
    return 0;
}