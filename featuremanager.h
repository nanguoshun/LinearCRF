//
// Created by  ngs on 13/06/2018.
//
#ifndef CRF_FEATURE_H
#define CRF_FEATURE_H

#include "node.h"
#include "path.h"
#include <map>
#include <string>
#include <math.h>
#include <unordered_map>
#include "common.h"
#include <boost/functional/hash.hpp>

class Feature {
public:
    explicit Feature(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *tag_vector,
                     std::map<std::string, int> *ptr_x_corpus_map, std::map<int,std::string> *ptr_tag_map_reverse, int num_of_seq);

    Feature(const char *feature_file_name);

    ~Feature();

    //void CalcCost(Node *ptr_node);

    //void CalcCost(Path *ptr_path);
/*
    inline void CalcCost(Node *ptrnode) {
        double cost = DEFAULT_COST_VALUE;
        int observation = ptrnode->GetX();
        int y = ptrnode->GetY();
        int index = GetFeatureIndex(std::make_pair(observation,y));
        if(FEATURE_NO_EXIST != index){
            double weight = (*ptr_w_vector_)[index];
//        cost = exp(weight);
            cost = weight;
        } else{
            // std::cout << "index error"<<std::endl;
        }
        ptrnode->SetCost(cost);
    }*/
    inline void CalcCost(Node *ptrnode) {
        double cost = DEFAULT_COST_VALUE;
        int index = ptrnode->GetFeatureIndex();
        if(FEATURE_NO_EXIST != index){
            double weight = (*ptr_w_vector_)[index];
//        cost = exp(weight);
            cost = weight;
        } else{
            // std::cout << "index error"<<std::endl;
        }
        ptrnode->SetCost(cost);
    }
/*
    inline void CalcCost(Path *ptrpath) {
        double cost = DEFAULT_COST_VALUE;
        int lnodeY = ptrpath->GetLNode()->GetY();
        int rnodeY = ptrpath->GetRNode()->GetY();
        int index = GetFeatureIndex(std::make_pair(lnodeY, rnodeY));
        if (FEATURE_NO_EXIST != index) {
            //for the pair like (START, y)
            double weight = (*ptr_w_vector_)[index];
//        cost = exp(weight);
            cost = weight;
        } else{
            // std::cout << "index error"<<std::endl;
        }
        ptrpath->SetCost(cost);
    }
*/
    inline void CalcCost(Path *ptrpath) {
        double cost = DEFAULT_COST_VALUE;
        int index = ptrpath->GetFeatureIndex();
        if (FEATURE_NO_EXIST != index) {
            //for the pair like (START, y)
            double weight = (*ptr_w_vector_)[index];
            cost = weight;
        } else{
            // std::cout << "index error"<<std::endl;
        }
        ptrpath->SetCost(cost);
    }

    //for path: a and b denote the left node and the right node, respectively.
    //for node: a and b denote the node and the observation string, respectively.
    //int GetFeatureIndex(std::pair<int, int> str_pair);
    inline int GetFeatureIndex(std::pair<int, int> str_pair) {
        if(ptr_feature_map_->find(str_pair) != ptr_feature_map_->end()){
            return ptr_feature_map_->find(str_pair)->second;
        } else{
            return FEATURE_NO_EXIST;
        }
    }

    std::vector<double> *GetWeightVector();

    int GetFeatureSize();

    void CreateFeatureMap(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                          std::map<std::string, int> *ptr_x_corpus_map,
                          std::map<int, std::string> *ptr_tag_map_reverse, int num_of_seq);
    void CreateFeatureMap(const char *feature_file_name);
    void InsertFeature(std::pair<int ,int > feature_pair, int *index, int offset,  std::ofstream* ptr_of, bool istest);

//    std::unordered_map<std::pair<int,int>, int> *GetFeatureMap();

//    std::unordered_map<int,std::pair<int, int>> *GetReverseFeatureMap();

    std::unordered_map<std::pair<int ,int >, int,boost::hash<std::pair<int, int>>> *GetFeatureMap();

    std::unordered_map<int, std::pair<int ,int >> *GetReverseFeatureMap();
    void SetWeightVector(int index, double weight);

    int GetFeatureSizeEdge();

    int GetFeatureSizeNode();

private:
    //weight vector.
    std::vector<double> *ptr_w_vector_;
    std::unordered_map<std::pair<int ,int >, int, boost::hash<std::pair<int, int>>> *ptr_feature_map_;
    std::unordered_map<int, std::pair<int ,int >> *ptr_reverse_feature_map_;
    //std::unordered_map<std::pair<int ,int >, int> *ptr_feature_map_;
    //std::map<int, std::pair<int ,int >> *ptr_reverse_feature_map_;
//    std::vector<double> *ptr_f_e_;
//    std::vector<double> *ptr_f_empirical_e_;
    int feature_size_;
    int feature_size_edge_;
    int feature_size_node_;
    int feature_size_start_stop_;

};


#endif //CRF_FEATURE_H
