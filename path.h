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
    inline double LogSumExp(double x, double y, bool isStart){
        // calc \alpha * cost;
        //      double value = cost * lalpha;
        if(isStart){
            return y; // init mode
        }
        const double vmin = std::min(x,y);
        const double vmax = std::max(x,y);
        if(vmax > vmin + MINUS_LOG_EPSILON)
        {
            return vmax;
        }else{
            return vmax + std::log(std::exp(vmin-vmax)+1.0);
        }
    }
    //
    void CalcExpectation(double Z,std::vector<double> *ptr_expectation,double *ppath_expectation);
    void CalcLogExpectation(double Z,std::vector<double> *ptr_expectation,double *ppath_expectation);

    //
    void AddNode(Node *ptr_lnode, Node *ptr_rnode);
    void SetFeatureIndex(int index);
    int GetFeatureIndex();
    void SetPathID(std::pair<int, int> id);
    std::pair<int, int> GetPathID();
    bool isCalculated();
    void SetCalculatedFlag(bool isCalc);
    bool isSetFeature();
private:

    Node *ptr_lnode_; //left node on an edge
    Node *ptr_rnode_; //right node on an edge
    //cost indicates the w_k * f_k(s', s, x) for a path (edge)
    double cost_;
    //
    double expectation_;
    int feature_index_;
    std::pair<int, int> path_id_;
    bool is_calculated_;
    bool  isSetFeature_;
};

#endif //CRF_PATH_H
