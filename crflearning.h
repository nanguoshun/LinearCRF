//
// Created by  ngs on 09/06/2018.
//
#ifndef CRF_CRF_LEARNING_H
#define CRF_CRF_LEARNING_H
#include <iostream>
#include <vector>
#include "common.h"
#include "node.h"
#include "feature.h"
#include "datasetmgr.h"

class CRFLearning{
public:
    explicit CRFLearning(DatasetMgr *ptr_datamgr);
    ~CRFLearning();
    void Init();
    double CalcLossFunction(std::vector<std::string> seq);
    void Learning();
    void CalcGradient();
    double CalcExpectedFi(std::vector<std::string> seq);
    double CalcEmpericialFi();
    void ForwardBackward(std::vector<std::string> seq);
    void BuildLattice(std::vector<std::string> seq);
    void DeleteLattice(std::vector<std::string> seq);
private:
    DatasetMgr *ptr_datamgr_;
    double loss_value_;
    std::vector<std::vector<Node *>> node_matrix_;
    Node *ptr_start_node_;
    Node *ptr_stop_node_;
    int tag_num_;
    double Z_;
    Feature *ptr_feature_;
    std::vector<double> *ptr_gradient_;
    //data related
    std::vector<std::string> *ptr_x_vector_;
    std::set<std::string> *ptr_x_set_;
    std::vector<std::string> *ptr_tag_vector_;
    std::set<std::string> *ptr_tag_set_;
    std::map<std::string, int> *ptr_tag_map_;

    std::map<std::string, int> *ptr_x_corpus_map_;
    std::vector<std::string> *ptr_x_corpus_;

};

#endif //CRF_CRF_LEARNING_H
