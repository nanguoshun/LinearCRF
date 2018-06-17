#include <iostream>
#include "crflearning.h"
int main(int argc, char **argv) {
    DatasetMgr *ptr_datamgr = new DatasetMgr(true);
    ptr_datamgr->OpenDataSet(argv[1], true);
    CRFLearning *ptr_crf = new CRFLearning(ptr_datamgr);
    ptr_crf->Learning();
    return 0;
}