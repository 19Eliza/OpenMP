cmake_minimum_required(VERSION 3.22)
project(EllipsoidsTask)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 

find_package(OpenMP REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(EllipsoidsTask PUBLIC OpenMP::OpenMP_CXX)
endif()