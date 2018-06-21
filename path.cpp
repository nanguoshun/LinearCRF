//
// Created by  ngs on 13/06/2018.
//
#include "path.h"
#include "common.h"

Path::Path(Node *ptr_lnode, Node *ptr_rnode) {
    ptr_lnode_ = ptr_lnode;
    ptr_rnode_ = ptr_rnode;
    cost_ = 0;
    expectation_ = 0;
    feature_index_ = FEATURE_NO_EXIST;
    is_calculated_ = false;
    isSetFeature_ = false;
}

Path::~Path() {

}

void Path::CalcExpectation(double Z, std::vector<double> *ptr_expectation, double *ppath_expectation) {
    double alpha = ptr_lnode_->GetAlpha();
    double beta = ptr_rnode_->GetBeta();
    double cost = cost_ * ptr_lnode_->GetCost();
    expectation_  = (alpha * beta * cost)/Z;
    if(feature_index_!=FEATURE_NO_EXIST){
        (*ptr_expectation)[feature_index_] += expectation_;
       // std::cout << "The expection of feature " << feature_index_ <<" is: "<<(*ptr_expectation)[feature_index_]<<std::endl;
    } else{
        std::cout << "error" <<std::endl;
    }
    (*ppath_expectation) += expectation_;
}

void Path::AddNode(Node *ptr_lnode, Node *ptr_rnode) {

}

double Path::GetCost() {
    return cost_;
}

Node* Path::GetLNode() {
    return ptr_lnode_;
}

Node* Path::GetRNode() {
    return ptr_rnode_;
}

void Path::SetCost(double cost) {
    cost_ = cost;
}

int Path::GetFeatureIndex() {
    return feature_index_;
}

void Path::SetFeatureIndex(int index) {
    feature_index_ = index;
    isSetFeature_ = true;
}

void Path::SetPathID(std::pair<int, int> id) {
    path_id_ = id;
}

std::pair<int, int> Path::GetPathID() {
    return path_id_;
}

bool Path::isCalculated() {
    return is_calculated_;
}

void Path::SetCalculatedFlag(bool isCalc) {
    is_calculated_ = isCalc;
}

bool Path::isSetFeature() {
    return isSetFeature_;
}