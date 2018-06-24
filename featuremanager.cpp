//
// Created by  ngs on 13/06/2018.
//
#include "featuremanager.h"
#include <iostream>
#include <tuple>
#include <fstream>
#include <sstream>

Feature::Feature(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                 std::map<std::string, int> *ptr_x_corpus_map, std::map<int,std::string> *ptr_tag_map_reverse, int num_of_seq) {
    ptr_feature_map_ = new std::map<std::pair<int, int>, int>;
    ptr_reverse_feature_map_ = new std::map<int, std::pair<int, int>>;
    CreateFeatureMap(ptr_observ_vector, ptr_tag_vector, ptr_x_corpus_map, ptr_tag_map_reverse, num_of_seq);
    ptr_w_vector_ = new std::vector<double>(ptr_feature_map_->size());
    std::fill(ptr_w_vector_->begin(),ptr_w_vector_->end(),0);
    //for test only
/*
    double array[12] = {0.65,0.2,0.5,0.25,0.45,0.15,0.35,0.4,0.6,0.3,0.1,0.15};
    for(int i=0; i<ptr_feature_map_->size(); ++i){
        (*ptr_w_vector_)[i] = array[i];
    }
*/
}

Feature::Feature(const char *feature_file_name) {
    ptr_feature_map_ = new std::map<std::pair<int, int>, int>;
    ptr_reverse_feature_map_ = new std::map<int, std::pair<int, int>>;
    CreateFeatureMap(feature_file_name);
    ptr_w_vector_ = new std::vector<double>(ptr_feature_map_->size());
    std::fill(ptr_w_vector_->begin(),ptr_w_vector_->end(),0);
}

Feature::~Feature() {
    delete ptr_w_vector_;
    delete ptr_feature_map_;
    delete ptr_reverse_feature_map_;
}

//create feature map from feature file during test phase
void Feature::CreateFeatureMap(const char *feature_file_name) {
    std::ifstream ifs(feature_file_name);
    std::string str;
    while(getline(ifs,str)){
        std::stringstream ss(str);
        int x;
        int index = 0;
        int data[3]; // data[0], data[1], data[2] are feature index, y_i_1(x_i) and y_i, respectively.
        while(ss >> x){
            data[index++] = x;
        }
        InsertFeature(std::make_pair(data[1],data[2]),&data[0],0,NULL,true);
    }
}
//create feature map from training dataset during training phase.
void Feature::CreateFeatureMap(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                               std::map<std::string, int> *ptr_x_corpus_map,
                               std::map<int, std::string> *ptr_tag_map_reverse, int num_of_seq) {
    std::ofstream *ptr_of = new std::ofstream("featuremap.txt");
    int index_offset = 0;
    for (int i = 0; i < ptr_tag_map_reverse->size(); ++i) {
        for (int j = 0; j < ptr_tag_map_reverse->size(); ++j) {
            InsertFeature(std::make_pair(i,j),&index_offset, 0, ptr_of,false);
        }
    }
    feature_size_edge_ = index_offset;
    index_offset = 0;
    int cout = 0;
    for (int i = 0; i < ptr_observ_vector->size(); ++i) {
        cout++;
        for (int j = 0; j < ptr_tag_map_reverse->size(); ++j) {
            int ob1 = ptr_x_corpus_map->find((*ptr_observ_vector)[i])->second + FEATURE_CODE_OFFSET;
            InsertFeature(std::make_pair(ob1, j), &index_offset, feature_size_edge_, ptr_of, false);
        }
    }
    feature_size_node_= index_offset;
    index_offset = 0;
    //for start_node
    int start_node_flg = START_NODE_FLAG;
    for (int seq_no = 0; seq_no < num_of_seq; ++seq_no) {
        for (int k = 0; k < ptr_tag_map_reverse->size() ; ++k) {
            InsertFeature(std::make_pair(start_node_flg,k),&index_offset,feature_size_edge_+feature_size_node_, ptr_of, false);
        }
        start_node_flg--;
    }
    //for stop_node
    int stop_node_flg = STOP_NODE_FLAG;
    for (int seq_no = 0; seq_no < num_of_seq; ++seq_no) {
        for (int k = 0; k < ptr_tag_map_reverse->size(); ++k) {
            InsertFeature(std::make_pair(k, stop_node_flg), &index_offset, feature_size_edge_ + feature_size_node_, ptr_of, false);
        }
        stop_node_flg--;
    }
    feature_size_start_stop_ = index_offset;
    feature_size_ = feature_size_edge_ + feature_size_node_ + feature_size_start_stop_;
    delete ptr_of;
}

void Feature::InsertFeature(std::pair<int, int> feature_pair, int *index, int offset, std::ofstream *ptr_of, bool istest) {
    if (ptr_feature_map_->find(feature_pair) == ptr_feature_map_->end()) {
        ptr_reverse_feature_map_->insert(std::make_pair((*index) + offset, feature_pair));
        ptr_feature_map_->insert(std::make_pair(feature_pair, (*index) + offset));
        //write file in training phase
        if(!istest){
            (*ptr_of) << std::to_string((*index) + offset) + " " + std::to_string(feature_pair.first) + " " +
                         std::to_string(feature_pair.second);
            (*ptr_of) << std::endl;
            (*index)++;
        }
    }
}

void Feature::CalcCost(Node *ptrnode) {
    double cost = 1;
    int observation = ptrnode->GetX();
    int y = ptrnode->GetY();
    int index = GetFeatureIndex(std::make_pair(observation,y));
    if(FEATURE_NO_EXIST != index){
        double weight = (*ptr_w_vector_)[index];
        cost = exp(weight);
    } else{
        std::cout << "index error"<<std::endl;
    }
    ptrnode->SetCost(cost);
}

void Feature::CalcCost(Path *ptrpath) {
    double cost = 1;
    int lnodeY = ptrpath->GetLNode()->GetY();
    int rnodeY = ptrpath->GetRNode()->GetY();
    int index = GetFeatureIndex(std::make_pair(lnodeY, rnodeY));
    if (FEATURE_NO_EXIST != index) {
        //for the pair like (START, y)
        double weight = (*ptr_w_vector_)[index];
        cost = exp(weight);
    } else{
        std::cout << "index error"<<std::endl;
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