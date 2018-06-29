//
// Created by  ngs on 27/04/2018.
//

#ifndef CLIONTEST_COMMON_H
#define CLIONTEST_COMMON_H

#include <iostream>
#include <vector>
#include <time.h>


const int LINE_MAX_SIZE = 18192;
const int LINE_FEATURE_MAX_SIZE = 18192;
const int TAG_MAX_SIZE = 1024;
const int DIGITAL_FLAG = 123456789;
const double SMOOTH_VALUE = 0.0001;
const std::string SPACE_OF_LINE = "\t";
const std::string OUT_FLAG = "OUT";
const std::string SPERATOR_FLAG = "+";
const std::string TAGER_BIO_B = "B";
const std::string TAGER_BIO_I = "I";
const std::string TAGER_BIO_O = "O";
const std::string NOT_FOUND_FLAG = "NOT_FOUND_FLAG";
const std::string NOT_FOUND = "NOT_FOUND";
const std::string PATH_START = "PATH_START";
const double EM_ITERATION_STOP = 0.1;
const double INITIAL_LOG_LIKEIHOOD = 10000;
const double ALPHA_START = 0.01;
const int RAND_MAX_NUM = 100;
const double DEFAULT_COST_VALUE = 0;
const double DEFAULT_COST_VALUE_DEC = 1;
const int MINUS_LOG_EPSILON = 50;

const double LEARNING_RATE = 0.05;
const double L2_FACTOR = 0.1;
const double CONVERGED_VALUE = -0.1;
//const double CONVERGED_VALUE = -9.54;
const int FEATURE_CODE_OFFSET = 500000;
const int FEATURE_NO_EXIST = -100;
const int START_NODE_FLAG = -100000;
const int STOP_NODE_FLAG = -200000;

#define TEST_MODE = 100;

const int START_NODE_ID = 100000;
const int STOP_NODE_ID =  200000;

//#define DEBUG_MODE = 200;
//#define DEBUG_MODE_ = 200;
//#define PERFORMANCE_CLOCK_TIME;
//#define PERFORMANCE_CLOCK_TIME_;

#define DEBUG_MODE__

struct HMMParameters{
    size_t  num_of_state_;
    std::vector<std::vector<double >> *ptr_e_;
    std::vector<std::vector<double >> *ptr_t_;
    std::vector<std::vector<double >> *ptr_count_uv_;
    std::vector<double> *ptr_count_u_;
    std::vector<std::vector<double>> *ptr_count_uk_;
    //for test only
    size_t num_of_x_;
};

#endif //CLIONTEST_COMMON_H
