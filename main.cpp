#include "crf.h"
#include "decoder.h"
int main(int argc, char **argv) {
    bool is_traing = true;
    bool is_decoding = true;
    bool is_calc_result = false;
    if(is_calc_result){
        Decoder *ptr_decoder = new Decoder();
        ptr_decoder->CalculateResult();
//        ptr_decoder->RewriteTrainandTestData("test.data","newtest.data");
    } else{
        if(!is_decoding){
            DatasetMgr *ptr_datamgr_training = new DatasetMgr(true);
            DatasetMgr *ptr_datamgr_test = new DatasetMgr(true);
            ptr_datamgr_training->OpenDataSet(argv[1], argv[2],true);
            ptr_datamgr_test->OpenDataSet(argv[1], argv[2],false);
            LinearCRF *ptr_crf = new LinearCRF(ptr_datamgr_training);
            ptr_crf->Training();
        }else{
            DatasetMgr *ptr_datamgr_test1 = new DatasetMgr(true);
            ptr_datamgr_test1->OpenDataSet(argv[1], argv[2],false);
            Decoder *ptr_decoder = new Decoder(ptr_datamgr_test1);
            ptr_decoder->Decoding();
        }
        std::cout << "pls check the file"<<std::endl;

    }
    return 0;
}