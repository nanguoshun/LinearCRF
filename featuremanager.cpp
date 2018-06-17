//
// Created by  ngs on 13/06/2018.
//
#include "featuremanager.h"

Feature::Feature(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                 std::map<std::string, int> *ptr_x_corpus_map, std::map<std::string, int> *ptr_tag_map) {
    ptr_f_e_ = new std::vector<double>;
    ptr_f_empirical_e_ = new std::vector<double>;
    ptr_feature_map_ = new std::map<std::pair<int, int>, int>;
    ptr_reverse_feature_map_ = new std::map<int, std::pair<int, int>>;
    CreateFeatureMap(ptr_observ_vector,ptr_tag_vector,ptr_x_corpus_map,ptr_tag_map);
    ptr_w_vector_ = new std::vector<double>(ptr_feature_map_->size());
    std::fill(ptr_w_vector_->begin(),ptr_w_vector_->end(),0);
}

Feature::~Feature() {
    delete ptr_w_vector_;
    delete ptr_feature_map_;
    delete ptr_reverse_feature_map_;
    delete ptr_f_e_;
    delete ptr_f_empirical_e_;
}

void Feature::CreateFeatureMap(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                               std::map<std::string, int> *ptr_x_corpus_map, std::map<std::string, int> *ptr_tag_map) {
    int index_offset = 0;
    //create transition features.
    for (int i = 0; i < ptr_tag_vector->size()-1; ++i) {
        int tag1 = ptr_tag_map->find((*ptr_tag_vector)[i])->second;
        int tag2 = ptr_tag_map->find((*ptr_tag_vector)[i + 1])->second;
        std::pair<int,int> tag_pair = std::make_pair(tag1, tag2);
        //to avoid duplicate
        if (ptr_feature_map_->find(tag_pair) == ptr_feature_map_->end()) {
            ptr_reverse_feature_map_->insert(std::make_pair(index_offset, tag_pair));
            ptr_feature_map_->insert(std::make_pair(tag_pair, index_offset));
            index_offset++;
        }
    }
    feature_size_edge_ = index_offset;
    index_offset = 0;
    //create emission features;
    for (int i = 0; i < ptr_tag_vector->size(); ++i) {
        int ob1 = ptr_x_corpus_map->find((*ptr_observ_vector)[i])->second + FEATURE_CODE_OFFSET;
        int tag1 = ptr_tag_map->find((*ptr_tag_vector)[i])->second;
        std::pair<int,int> tag_obsv_pair = std::make_pair(ob1,tag1);
        if (ptr_feature_map_->find(tag_obsv_pair) == ptr_feature_map_->end()) {
            ptr_reverse_feature_map_->insert(std::make_pair((index_offset + feature_size_edge_), tag_obsv_pair));
            ptr_feature_map_->insert(std::make_pair(tag_obsv_pair, (index_offset + feature_size_edge_)));
            index_offset++;
        }
    }
    feature_size_node_= index_offset;
    feature_size_ = feature_size_edge_ + feature_size_node_;
}

void Feature::CalcCost(Node *ptrnode) {
    double cost = 0;
    int observation = ptrnode->GetX();
    int y = ptrnode->GetY();
    int index = GetFeatureIndex(std::make_pair((observation+FEATURE_CODE_OFFSET),y));
    if(FEATURE_NO_EXIST != index){
        double weight = (*ptr_w_vector_)[index];
        cost = exp(weight);
    }
    ptrnode->SetCost(cost);
}

void Feature::CalcCost(Path *ptrpath) {
    double cost = 0;
    int lnodeY = ptrpath->GetLNode()->GetY();
    int rnodeY = ptrpath->GetRNode()->GetY();
    int index = GetFeatureIndex(std::make_pair(lnodeY, rnodeY));
    if (FEATURE_NO_EXIST != index) {
        //for the pair like (START, y)
        double weight = (*ptr_w_vector_)[index];
        cost = exp(weight);
    }
    if (lnodeY == START_NODE_FLAG) {
        cost = exp(1);
    }
    if (rnodeY == STOP_NODE_FLAG){
        cost = exp(1);
    }
    ptrpath->SetCost(cost);
}

int Feature::GetFeatureIndex(std::pair<int, int> str_pair) {
    if(ptr_feature_map_->find(str_pair) != ptr_feature_map_->end()){
        return ptr_feature_map_->find(str_pair)->second;
    } else{
        return FEATURE_NO_EXIST;
    }
}

std::vector<double>* Feature::GetWeightVector() {
    return ptr_w_vector_;
}

int Feature::GetFeatureSize() {
    return  feature_size_;
}

std::map<std::pair<int, int>, int>* Feature::GetFeatureMap() {
    return  ptr_feature_map_;
}

std::map<int, std::pair<int, int>> * Feature::GetReverseFeatureMap() {
    return  ptr_reverse_feature_map_;
}

void Feature::SetWeightVector(int index, double weight) {
    (*ptr_w_vector_)[index] = weight;
}

int Feature::GetFeatureSizeEdge() {
    return feature_size_edge_;
}

int Feature::GetFeatureSizeNode() {
    return feature_size_node_;
}