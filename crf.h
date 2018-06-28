//
// Created by  ngs on 09/06/2018.
//
#ifndef CRF_CRF_LEARNING_H
#define CRF_CRF_LEARNING_H
#include <iostream>
#include <vector>
#include "common.h"
#include "node.h"
#include "featuremanager.h"
#include "datasetmgr.h"
#include <math.h>
#include <thread>
#include "crfthread.h"

class LinearCRF{
public:
    explicit LinearCRF(DatasetMgr *ptr_datamgr_training);
    void AllocateSpace();
    void CreateTagObservMap();
    ~LinearCRF();
    void Init();
    double CalcLoglikelihoodFunction();
    void Training();
    void FromVectorToSet();
    void CRFRun();
    void MainThreadCalculation();
    void ThreadCalc(CRFThread *ptr_tast,int seq_no);
    void ThreadStart();
    void CalcGradient();
    void CalcFeatureExpectation(std::vector<std::string> seq, int seq_no);
    double CalcEmpiricalFi(std::vector<std::string> x_seq, std::vector<std::string> tag_seq, int seq_no);
    void CalcAllEmpiricalFi();
    void BuildLattice();
    void BuildNode(std::vector<std::string> seq, int seq_no);
    void BuildLPath(std::vector<std::string> seq, int seq_no);
    void BuildRPath(std::vector<std::string> seq, int seq_no);
    void DeleteLattice();
    void UpdateWeight();
    void ResetParameters();
    void SetPathFeature(std::pair<int,int> feature_pair, Path *ppath);
    void GenerateSeqFromVector(std::vector<std::string> *ptr_vector,std::vector<std::vector<std::string>> *ptr_seq_vector);
    void SaveModelToFile();
    inline double LogSumExp(double x, double y, bool isStart){
        // calc \alpha * cost;
        //      double value = cost * lalpha;
        if(isStart){
            return y; // init mode
        }
        const double vmin = std::min(x,y);
        const double vmax = std::max(x,y);
        if(vmax > vmin + MINUS_LOG_EPSILON)
        {
            return vmax;
        }else{
            return vmax + std::log(std::exp(vmin-vmax)+1.0);
        }
    }


private:
    DatasetMgr *ptr_datamgr_;
    //DatasetMgr *ptr_datamgr_test_;
    std::vector<std::vector<std::vector<Node *>>> node_matrix_;
    std::vector<Node> *ptr_start_node_;
    std::vector<Node> *ptr_stop_node_;
    //int tag_num_;
    std::vector<double> *ptr_Z_;
    std::vector<double> *ptr_gradient_;
    //data related
    std::vector<std::vector<std::string>> *ptr_seq_matrix_;
    std::vector<std::vector<std::string>> *ptr_tag_seq_;
    std::vector<std::set<std::string>> *ptr_tag_set_matrix_;
    int num_of_instance_;
    std::vector<std::string> *ptr_x_train_vector_;
    //std::vector<std::string> *ptr_x_test_vector_;
    std::set<std::string> *ptr_x_traing_set_;
    //std::set<std::string> *ptr_x_test_set_;
    std::vector<std::string> *ptr_tag_vector_;
    std::set<std::string> *ptr_tag_set_;
    //std::set<std::string> *ptr_test_tag_set_;
    std::map<std::string, int> *ptr_tag_map_;
    std::map<int, std::string> *ptr_tag_map_reverse_;
    std::map<std::string, int> *ptr_x_corpus_map_;
    std::map<int, std::string> *ptr_x_corpus_map_reverse_;
    //std::vector<std::string> *ptr_x_corpus_;
    //feature related
    Feature *ptr_feature_;
    std::vector<std::pair<int, int>> *ptr_feature_vector_;
    std::vector<double> *ptr_empirical_e_;
    std::vector<double> *ptr_e_;
    bool is_converged_;
    std::vector<std::vector<std::string>> *ptr_decoded_tag_;
    std::vector<double> *ptr_feature_bit_vector_;
    //
    bool is_training_;
    // for multiple threading
    int num_of_thread_;
    std::vector<std::thread>* ptr_thread_vector_;

};

#endif //CRF_CRF_LEARNING_H
