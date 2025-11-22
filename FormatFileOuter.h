#include<iostream>
#include <fstream>
#include <iomanip>

using namespace std;

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
