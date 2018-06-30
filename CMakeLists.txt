cmake_minimum_required(VERSION 3.10)
project(CRF)

set(CMAKE_CXX_STANDARD 11)

add_executable(CRF main.cpp crf.h crf.cpp datasetmgr.h datasetmgr.cpp node.h node.cpp path.h path.cpp featuremanager.h featuremanager.cpp decoder.h decoder.cpp crfthread.h crfthread.cpp lbfgs.h lbfgs.cpp)