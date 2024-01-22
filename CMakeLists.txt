cmake_minimum_required(VERSION 3.20)
project(Discrete_Fourier_Transform)

set(CMAKE_CXX_STANDARD 17)

add_executable(Discrete_Fourier_Transform
        main.cpp
        DFT.hpp
        DFT.cpp)

find_package(Python3 COMPONENTS Development REQUIRED)
target_link_libraries(Discrete_Fourier_Transform PRIVATE Python3::Python)

target_include_directories(Discrete_Fourier_Transform PRIVATE "lib_inc")