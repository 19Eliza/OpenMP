#pragma once
#include"Ellipsoid.h"


using tripleEllipsoid = std::array<Ellipsoid, 3>; /// {E1,E2,E3}, brain(E1), skull(E2), skin(E3)

// Тестовые случаи
struct TestCase{
    inline static const tripleEllipsoid test1{ Ellipsoid{30,40,50}, Ellipsoid{35,45,55}, Ellipsoid{40,50,60} };
    inline static const tripleEllipsoid test2{Ellipsoid {30.0, 40.0, 50.0}, Ellipsoid {30.5, 40.5, 50.5}, Ellipsoid {31.0, 41.0, 51.0} };
    inline static const std::vector<tripleEllipsoid> testEllipsoid{test1,test2};// множество тестовых эллипсоидов
    static const int countPointX = 32;
    static const int countPointY = 32;
    static const int countPointZ = 32;
    static const int Nsamples{2048}; // количество точек разбросанных в одном элементе

};