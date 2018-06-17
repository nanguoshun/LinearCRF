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
}

Path::~Path() {

}

void Path::CalcExpectation(double Z, std::vector<double> *ptr_expectation) {
    double alpha = ptr_lnode_->GetAlpha();
    double beta = ptr_rnode_->GetBeta();
    expectation_  = (alpha * beta * cost_)/Z;
    if(feature_index_!=FEATURE_NO_EXIST){
        (*ptr_expectation)[feature_index_] += expectation_;
    }
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
}