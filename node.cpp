//
// Created by  ngs on 12/06/2018.
//

#include "node.h"
#include "path.h"
#include <math.h>

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
    lpath_.clear();
}


// calc alpha recursively
void Node::CalcAlpha() {
    alpha_ = 0;
    for(std::vector<Path *>::iterator it = lpath_.begin(); it!=lpath_.end();++it){
        alpha_ = SumExp(alpha_, (*it)->GetCost()+cost_, (*it)->GetLNode()->GetAlpha(), (it == lpath_.begin()));
    }
}
// calc beta recursively
void Node::CalcBeta() {
    beta_ = 0;
    for(std::vector<Path *>::iterator it = rpath_.begin(); it!=rpath_.end();++it){
        beta_ = SumExp(beta_, (*it)->GetCost()+cost_, (*it)->GetLNode()->GetAlpha(), (it == rpath_.begin()));
    }
}

double Node::GetAlpha() {
    return alpha_;
}

double Node::GetBeta() {
    return beta_;
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

void Node::CalcExpectation(double Z) {
    double value = alpha_  * beta_ * cost_;
    expectation_ = value / Z;
    for(std::vector<Path *>::iterator it = lpath_.begin();it!=lpath_.end();++it){
        (*it)->CalcExpectation(Z);
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