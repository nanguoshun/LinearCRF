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
    expectation_ = 0;
    feature_index_ = FEATURE_NO_EXIST;
    bestCost_ = 0;
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
        alpha_ = LogSumExp(alpha_, (*it)->GetCost() + (*it)->GetLNode()->GetCost(), (*it)->GetLNode()->GetAlpha(), (it == lpath_.begin()));
    }
}
// calc beta recursively
void Node::CalcBeta() {
    beta_ = 0;
    for(std::vector<Path *>::iterator it = rpath_.begin(); it!=rpath_.end();++it){
        beta_ = LogSumExp(beta_, (*it)->GetCost() + (*it)->GetLNode()->GetCost(), (*it)->GetRNode()->GetBeta(), (it == rpath_.begin()));
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
    double expectation = 0;
    for(std::vector<Path *>::iterator it = rpath_.begin();it!=rpath_.end();++it) {
        (*it)->CalcExpectation(Z, ptr_expectation, &expectation);
    }
    if(feature_index_ != FEATURE_NO_EXIST){
        (*ptr_expectation)[feature_index_] += expectation;
         //std::cout << "The expection of feature " << feature_index_ <<" is: "<<(*ptr_expectation)[feature_index_]<<std::endl;
    }else{
        //std::cout << "error"<<std::endl;
    }
}

void Node::CalcLogExpectation(double Z, std::vector<double> *ptr_expectation) {
    double expectation = 0;
    for (std::vector<Path *>::iterator it = rpath_.begin(); it != rpath_.end(); ++it) {
        (*it)->CalcLogExpectation(Z, ptr_expectation, &expectation);
    }
    if(feature_index_ != FEATURE_NO_EXIST){
        (*ptr_expectation)[feature_index_] += expectation;//LogSumExp((*ptr_expectation)[feature_index_],expectation,((*ptr_expectation)[feature_index_]==0));
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

void Node::SetIsPrintAllPath(bool isprint) {
    isPrintAllPath_ = isprint;
}

void Node::SetNodeID(int id) {
    node_id_ = id;
}

int Node::GetNodeID() {
    return node_id_;
}

double Node::GetPreCost() {
    return pre_cost_;
}

void Node::SetPreCost(double cost) {
    pre_cost_ = cost;
}