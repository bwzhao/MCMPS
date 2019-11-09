cmake_minimum_required(VERSION 3.15)
project(MCMPS)

set(CMAKE_CXX_STANDARD 17)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

set(CMAKE_CXX_FLAGS "-Wall -larmadillo -DARMA_DONT_USE_WRAPPER -framework Accelerate ${CMAKE_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "-Wall -larmadillo ${CMAKE_CXX_FLAGS}")

add_executable(MCMPS main.cpp Config.h Class_MPSTensor.cpp Class_MPSTensor.h Class_PySpace.cpp Class_PySpace.h Class_Network.cpp Class_Network.h Class_Matrix.cpp Class_Matrix.h Class_DataMeasurement.h)