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

class CRFLearning{
public:
    explicit CRFLearning(DatasetMgr *ptr_datamgr);
    ~CRFLearning();
    void Init(std::vector<std::string> seq);
    double CalcLoglikelihoodFunction(std::vector<std::string> seq);
    void Learning();
    void CalcGradient(std::vector<std::string> seq);
    void CalcCost(std::vector<std::string> seq);
    void CalcFeatureExpectation(std::vector<std::string> seq);
    void CalcNodeFeatureExpectation(int index, std::pair<int, int> feature);
    void CalcEdgeFeatureExpectation(int index, std::pair<int, int> feature);
    double CalcEmpiricalFi(std::vector<std::string> seq);
    void ForwardBackward(std::vector<std::string> seq);
    void BuildLattice(std::vector<std::string> seq);
    void BuildNode(std::vector<std::string> seq);
    void BuildLPath(std::vector<std::string> seq);
    void BuildRPath(std::vector<std::string> seq);
    void DeleteLattice(std::vector<std::string> seq);
    void UpdateWeight();
    void Viterbi();
private:
    DatasetMgr *ptr_datamgr_;
    double loss_value_;
    std::vector<std::vector<Node *>> node_matrix_;
    Node *ptr_start_node_;
    Node *ptr_stop_node_;
    int tag_num_;
    double Z_;
    std::vector<double> *ptr_gradient_;
    //data related
    std::vector<std::string> *ptr_x_vector_;
    std::set<std::string> *ptr_x_set_;
    std::vector<std::string> *ptr_tag_vector_;
    std::set<std::string> *ptr_tag_set_;
    std::map<std::string, int> *ptr_tag_map_;
    std::map<std::string, int> *ptr_x_corpus_map_;
    std::vector<std::string> *ptr_x_corpus_;
    //feature related
    Feature *ptr_feature_;
    std::vector<std::pair<int, int>> *ptr_feature_vector_;
    std::vector<double> *ptr_empirical_e_;
    std::vector<double> *ptr_e_;
    bool is_initialized_;
    bool is_converged_;
};

#endif //CRF_CRF_LEARNING_H
