//
// Created by  ngs on 12/06/2018.
//

#include "node.h"
#include "path.h"
#include <math.h>
#include "common.h"

Node::Node(int x, int y) {
    x_ = x;
    y_ = y;
    cost_ = 0;
    alpha_ = 0;
    beta_ = 0;
    expectation_ =0;
}

Node::~Node() {
    for (std::vector<Path *>::iterator it = lpath_.begin(); it!=lpath_.end();++it){
        delete (*it);
    }
    for (std::vector<Path *>::iterator it = rpath_.begin(); it!=rpath_.end();++it){
        delete (*it);
    }
}


// calc alpha recursively
void Node::CalcAlpha() {
    alpha_ = 0;
    for(std::vector<Path *>::iterator it = lpath_.begin(); it!=lpath_.end();++it){
        alpha_ = SumExp(alpha_, (*it)->GetCost() + cost_, (*it)->GetLNode()->GetAlpha(), (it == lpath_.begin()));
    }
}
// calc beta recursively
void Node::CalcBeta() {
    beta_ = 0;
    for(std::vector<Path *>::iterator it = rpath_.begin(); it!=rpath_.end();++it){
        beta_ = SumExp(beta_, (*it)->GetCost() + cost_, (*it)->GetRNode()->GetBeta(), (it == rpath_.begin()));
    }
}

double Node::GetAlpha() {
    return alpha_;
}

double Node::GetBeta() {
    return beta_;
}

void Node::SetAlpha(double alpha) {
    alpha_ = alpha;
}

void Node::SetBeta(double beta) {
    beta_ = beta;
}

int Node::GetX() {
    return x_;
}

int Node::GetY() {
    return y_;
}

void Node::SetCost(double cost) {
   cost_ = cost;
}

double Node::GetCost() {
    return cost_;
}

std::vector<Path *> Node::GetLPath() {
    return lpath_;
}

std::vector<Path *> Node::GetRPath() {
    return rpath_;
}

void Node::CalcExpectation(double Z, std::vector<double> *ptr_expectation){
    double value = alpha_  * beta_ * cost_;
    expectation_ = value / Z;
    if(feature_index_ != FEATURE_NO_EXIST){
        (*ptr_expectation)[feature_index_] += expectation_;
    }
    for(std::vector<Path *>::iterator it = lpath_.begin();it!=lpath_.end();++it){
        (*it)->CalcExpectation(Z, ptr_expectation);
    }
}

double Node::GetExpectation() {
    return expectation_;
}

void Node::AddPath(Path *ptr_path, bool isLpath) {
    if(isLpath){
        lpath_.push_back(ptr_path);
    } else{
        rpath_.push_back(ptr_path);
    }
}

void Node::SetFeatureIndex(int index) {
    feature_index_ = index;
}

int Node::GetFeatureIndex() {
    return  feature_index_;
}

double Node::GetBestCost() {
    return bestCost_;
}

void Node::SetBestCost(double cost) {
    bestCost_ = cost;
}

void Node::SetPreNode(Node *pNode) {
    ptr_pre_node_ = pNode;
}

Node* Node::GetPreNode() {
    return ptr_pre_node_;
}