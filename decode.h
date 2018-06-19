//
// Created by  ngs on 18/06/2018.
//
#include "crflearning.h"
#ifndef CRF_DECODE_H
#define CRF_DECODE_H
class Decode{
public:
    explicit Decode();
    void Inference();
private:
    CRFLearning *ptr_learning;
};
#endif //CRF_DECODE_H
