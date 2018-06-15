//
// Created by  ngs on 13/06/2018.
//

#ifndef CRF_PATH_H
#define CRF_PATH_H

#include "node.h"

class Path {
public:
    explicit Path(Node *ptr_lnode, Node *r_node);

    ~Path();

    Node *GetLNode();

    Node *GetRNode();

    double GetCost();

    void SetCost(double cost);

    //
    void CalcExpectation(double Z);

    //
    void AddNode(Node *ptr_lnode, Node *ptr_rnode);

private:

    Node *ptr_lnode_; //left node on an edge
    Node *ptr_rnode_; //right node on an edge
    //cost indicates the w_k * f_k(s', s, x) for a path (edge)
    double cost_;
    //
    double expectation_;
};

#endif //CRF_PATH_H
