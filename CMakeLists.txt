cmake_minimum_required(VERSION 3.10)
project(CRF)

set(CMAKE_CXX_STANDARD 11)

add_executable(CRF main.cpp crflearning.h crflearning.cpp datasetmgr.h datasetmgr.cpp node.h node.cpp path.h path.cpp feature.h feature.cpp)