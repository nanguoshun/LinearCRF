//
// Created by  ngs on 13/06/2018.
//
#include "feature.h"

Feature::Feature(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                 std::map<std::string, int> *ptr_x_corpus_map, std::map<std::string, int> *ptr_tag_map) {
    int offset = 0;
    //create transition features.
    for (int i = 0; i < ptr_tag_vector->size(); ++i) {
        int tag1 = ptr_tag_map->find((*ptr_tag_vector)[i])->second;
        int tag2 = ptr_tag_map->find((*ptr_tag_vector)[i + 1])->second;
        std::pair tag_pair = std::make_pair(tag1, tag2);
        //to avoid duplicate
        if (ptr_feature_map_->find(tag_pair) == ptr_feature_map_->end()) {
            ptr_feature_map_reverse_->insert(std::make_pair(i, tag_pair));
            ptr_feature_map_->insert(std::make_pair(tag_pair, i));
            offset++;
        }
    }
    //create emission features;
    for (int i = 0; i < ptr_tag_vector->size(); ++i) {
        int tag1 = ptr_tag_map->find((*ptr_tag_vector)[i])->second;
        int ob1 = ptr_x_corpus_map->find((*ptr_observ_vector)[i])->second;
        std::pair tag_obsv_pair = std::make_pair(tag1, ob1);
        if (ptr_feature_map_->find(tag_obsv_pair) == ptr_feature_map_->end()) {
            ptr_feature_map_reverse_->insert(std::make_pair((i + offset), tag_obsv_pair));
            ptr_feature_map_->insert(std::make_pair(tag_obsv_pair, (i + offset)));
        }
    }
}

void Feature::CalcCost(Node *ptrnode) {
    int y = ptrnode->GetY();
    int observation = ptrnode->GetX();
    int index = GetFeatureIndex(std::make_pair(y,observation));
    double cost = (*ptr_w_vector_)[index];
    ptrnode->SetCost(cost);
}

void Feature::CalcCost(Path *ptrpath) {
    int lnodeY = ptrpath->GetLNode()->GetY();
    int rnodeY = ptrpath->GetRNode()->GetY();
    int index = GetFeatureIndex(std::make_pair(lnodeY,rnodeY));
    double cost = (*ptr_w_vector_)[index];
    ptrpath->SetCost(cost);
}

int Feature::GetFeatureIndex(std::pair<int, int> str_pair) {
    return ptr_feature_map_->find(str_pair)->second;
}

std::vector<double>* Feature::GetWeightVector() {
    return ptr_w_vector_;
}

double Feature::CalcExpectation(int index) {
}

double Feature::GetFeatureExpectation(int index) {
    return (*ptr_f_e_)[index];
}

double Feature::GetEmpiricalFeatureExpectation(int index) {
    return (*ptr_f_empirical_e_)[index];
}

int Feature::GetFeatureSize() {
    return  feature_size_;
}