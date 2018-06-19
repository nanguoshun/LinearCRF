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
#include "common.h"

class Feature {
public:
    explicit Feature(std::vector<std::string> *observ_vector, std::vector<std::string> *tag_vector,
                     std::map<std::string, int> *ptr_x_corpus_map, std::map<int,std::string> *ptr_tag_map_reverse);

    ~Feature();

    void CalcCost(Node *ptr_node);

    void CalcCost(Path *ptr_path);

    //for path: a and b denote the left node and the right node, respectively.
    //for node: a and b denote the node and the observation string, respectively.
    int GetFeatureIndex(std::pair<int, int> str_pair);

    std::vector<double> *GetWeightVector();

    int GetFeatureSize();

    void CreateFeatureMap(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                          std::map<std::string, int> *ptr_x_corpus_map, std::map<std::string, int> *ptr_tag_map);

    void CreateAllFeatureMap(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                             std::map<std::string, int> *ptr_x_corpus_map,
                             std::map<int, std::string> *ptr_tag_map_reverse);

    std::map<std::pair<int, int>, int> *GetFeatureMap();

    std::map<int, std::pair<int, int>> *GetReverseFeatureMap();

    void SetWeightVector(int index, double weight);

    int GetFeatureSizeEdge();

    int GetFeatureSizeNode();

private:
    //weight vector.
    std::vector<double> *ptr_w_vector_;
    std::map<std::pair<int, int>, int> *ptr_feature_map_;
    std::map<int, std::pair<int, int>> *ptr_reverse_feature_map_;
//    std::vector<double> *ptr_f_e_;
//    std::vector<double> *ptr_f_empirical_e_;
    int feature_size_;
    int feature_size_edge_;
    int feature_size_node_;
    int feature_size_start_stop_;

};


#endif //CRF_FEATURE_H
