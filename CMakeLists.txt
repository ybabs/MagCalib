project(mag_calibration)
cmake_minimum_required(VERSION 3.11)

        
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)
set (CMAKE_PREFIX_PATH "/home/pc/matplotplusplus/")

find_package(Python3 COMPONENTS Interpreter Development NumPy)
find_package(OpenMP REQUIRED)
find_package(Matplot++ REQUIRED)

include_directories(${Python3_INCLUDE_DIRS} include)
include_directories("/usr/include/eigen3/")
# include_directories("/home/pc/matplotplusplus/" include)
# add_compile_options("$<$<CONFIG:RELEASE>:-W;-Wall;-O3;-pedantic>")
#add_compile_options("$<$<CONFIG:DEBUG>:-O0>")
add_compile_options("$<$<CONFIG:DEBUG>:-W;-Wall;-O0;-g;-pedantic;-fopenmp>")

add_executable(mag_calib src/main.cpp
                         src/mag.cpp)

target_link_libraries(mag_calib  ${Python3_LIBRARIES} OpenMP::OpenMP_CXX Matplot++::matplot)