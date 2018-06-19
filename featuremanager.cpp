//
// Created by  ngs on 13/06/2018.
//
#include "featuremanager.h"

Feature::Feature(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                 std::map<std::string, int> *ptr_x_corpus_map, std::map<int,std::string> *ptr_tag_map_reverse) {
    ptr_feature_map_ = new std::map<std::pair<int, int>, int>;
    ptr_reverse_feature_map_ = new std::map<int, std::pair<int, int>>;
    CreateAllFeatureMap(ptr_observ_vector,ptr_tag_vector,ptr_x_corpus_map,ptr_tag_map_reverse);
    ptr_w_vector_ = new std::vector<double>(ptr_feature_map_->size());
    std::fill(ptr_w_vector_->begin(),ptr_w_vector_->end(),0);
    //for test only
/*
    double array[12] = {0.65,0.2,0.5,0.25,0.45,0.15,0.35,0.4,0.6,0.3,0.1,0.15};
    for(int i=0; i<ptr_feature_map_->size(); ++i){
        (*ptr_w_vector_)[i] = array[i];
    }*/
}

Feature::~Feature() {
    delete ptr_w_vector_;
    delete ptr_feature_map_;
    delete ptr_reverse_feature_map_;
}

void Feature::CreateFeatureMap(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                               std::map<std::string, int> *ptr_x_corpus_map, std::map<std::string, int> *ptr_tag_map) {
    int index_offset = 0;
    //create transition features.
    for (int i = 0; i < ptr_tag_vector->size()-1; ++i) {
        int tag1 = ptr_tag_map->find((*ptr_tag_vector)[i])->second;
        int tag2 = ptr_tag_map->find((*ptr_tag_vector)[i + 1])->second;
        std::cout << "Feature index:"<<index_offset<<", the string and tag code are: "<<
                  (*ptr_tag_vector)[i]<<","<<(*ptr_tag_vector)[i+1]<<":"<<tag1<<","<<tag2<<std::endl;
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
        std::cout << "Feature index:"<<feature_size_edge_ + index_offset<<", the string and tag code are: "<<
                  (*ptr_observ_vector)[i]<<","<<(*ptr_tag_vector)[i]<<":"<<ob1<<","<<tag1<<std::endl;
        if (ptr_feature_map_->find(tag_obsv_pair) == ptr_feature_map_->end()) {
            ptr_reverse_feature_map_->insert(std::make_pair((index_offset + feature_size_edge_), tag_obsv_pair));
            ptr_feature_map_->insert(std::make_pair(tag_obsv_pair, (index_offset + feature_size_edge_)));
            index_offset++;
        }
    }
    feature_size_node_= index_offset;
    feature_size_ = feature_size_edge_ + feature_size_node_;
}

void Feature::CreateAllFeatureMap(std::vector<std::string> *ptr_observ_vector, std::vector<std::string> *ptr_tag_vector,
                                  std::map<std::string, int> *ptr_x_corpus_map,
                                  std::map<int,std::string> *ptr_tag_map_reverse) {
    int index_offset = 0;
    int size = ptr_tag_vector->size();
    for (int i = 0; i < ptr_tag_map_reverse->size(); ++i) {
        for (int j = 0; j < ptr_tag_map_reverse->size(); ++j) {
            //make pair referring to tag_map;
            std::pair<int, int> tag_pair = std::make_pair(i,j);
            ptr_reverse_feature_map_->insert(std::make_pair(index_offset,tag_pair));
            ptr_feature_map_->insert(std::make_pair(tag_pair,index_offset));
            std::cout << "Feature index:"<<index_offset<<", the tag and tag code are: "<<
                      ptr_tag_map_reverse->find(i)->second<<","<< ptr_tag_map_reverse->find(j)->second<<":"<<i<<","<<j<<std::endl;
            index_offset++;
        }
    }
    feature_size_edge_ = index_offset;
    index_offset = 0;
    for (int i = 0; i < ptr_observ_vector->size(); ++i) {
        for (int j = 0; j < ptr_tag_map_reverse->size(); ++j) {
            int ob1 = ptr_x_corpus_map->find((*ptr_observ_vector)[i])->second + FEATURE_CODE_OFFSET;
            std::pair<int, int> tag_obsv_pair = std::make_pair(ob1,j);
            if (ptr_feature_map_->find(tag_obsv_pair) == ptr_feature_map_->end()) {
                ptr_reverse_feature_map_->insert(std::make_pair((index_offset + feature_size_edge_), tag_obsv_pair));
                ptr_feature_map_->insert(std::make_pair(tag_obsv_pair, (index_offset + feature_size_edge_)));
                std::cout << "Feature index:"<<feature_size_edge_ + index_offset<<", the string and tag code are: "<<
                          (*ptr_observ_vector)[i]<<","<<ptr_tag_map_reverse->find(j)->second<<":"<<ob1<<","<<j<<std::endl;
                index_offset++;
            }
        }
    }
    feature_size_node_= index_offset;
    index_offset = 0;
    //for start_node
    for (int k = 0; k < ptr_tag_map_reverse->size() ; ++k) {
        int tag1 = START_NODE_FLAG;
        std::pair<int, int> tag_pair = std::make_pair(tag1,k);
        ptr_reverse_feature_map_->insert(std::make_pair(index_offset+feature_size_edge_+feature_size_node_,tag_pair));
        ptr_feature_map_->insert(std::make_pair(tag_pair,index_offset+feature_size_edge_+feature_size_node_));
        std::cout << "Feature index:"<<feature_size_edge_ + index_offset + feature_size_node_<<", the string and tag code are: "<<
                  "Start Node"<<","<<ptr_tag_map_reverse->find(k)->second<<":"<<tag1<<","<<k<<std::endl;
        index_offset++;
    }
    //for stop_node
    for (int k = 0; k < ptr_tag_map_reverse->size() ; ++k) {
        int tag2 = STOP_NODE_FLAG;
        std::pair<int, int> tag_pair = std::make_pair(k,tag2);
        ptr_reverse_feature_map_->insert(std::make_pair(index_offset+feature_size_edge_+feature_size_node_,tag_pair));
        ptr_feature_map_->insert(std::make_pair(tag_pair,index_offset+feature_size_edge_+feature_size_node_));
        std::cout << "Feature index:"<<feature_size_edge_ + index_offset+feature_size_node_<<", the string and tag code are: "<<
                  ptr_tag_map_reverse->find(k)->second<<","<<"Stop Node"<<":"<<k<<","<<tag2<<std::endl;

        index_offset++;
    }
    feature_size_start_stop_ = index_offset;
    feature_size_ = feature_size_edge_ + feature_size_node_ + feature_size_start_stop_;

}

void Feature::CalcCost(Node *ptrnode) {
    double cost = 1;
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
    double cost = 1;
    int lnodeY = ptrpath->GetLNode()->GetY();
    int rnodeY = ptrpath->GetRNode()->GetY();
    int index = GetFeatureIndex(std::make_pair(lnodeY, rnodeY));
    if (FEATURE_NO_EXIST != index) {
        //for the pair like (START, y)
        double weight = (*ptr_w_vector_)[index];
        cost = exp(weight);
    }
    /*
    if (lnodeY == START_NODE_FLAG) {
        cost = exp(0);
    }
    if (rnodeY == STOP_NODE_FLAG){
        cost = exp(0);
    }*/
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