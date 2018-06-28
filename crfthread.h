//
// Created by  ngs on 28/06/2018.
//

#ifndef CRF_CRFTHREAD_H
#define CRF_CRFTHREAD_H

#include <pthread.h>
#include "node.h"
#include "path.h"
#include "common.h"
#include <vector>
#include <map>
#include <string>
#include <thread>

class CRFThread{
public:
    CRFThread(std::vector<Node> *ptr_stratnode,std::vector<Node> *ptr_stoonode,std::vector<std::vector<std::vector<Node *>>> node_matrix, std::vector<double> *ptr_Z);
    ~CRFThread();
    void CRFThreadRun(int x_size, int y_size, int seq_no, std::vector<double> *ptr_w_vector);
    void ForwardBackward(int x_size, int y_size, int seq_no);
    void CalcCost(int x_size, int y_size, int seq_no,std::vector<double> *ptr_w_vector);
    inline void CalcCost(Path *ptrpath, std::vector<double> *ptr_w_vector) {
        double cost = DEFAULT_COST_VALUE;
        int index = ptrpath->GetFeatureIndex();
        if (FEATURE_NO_EXIST != index) {
            //for the pair like (START, y)
            double weight = (*ptr_w_vector)[index];
            cost = weight;
        } else{
            // std::cout << "index error"<<std::endl;
        }
        ptrpath->SetCost(cost);
    }
    inline void CalcCost(Node *ptrnode, std::vector<double> *ptr_w_vector) {
        double cost = DEFAULT_COST_VALUE;
        int index = ptrnode->GetFeatureIndex();
        if(FEATURE_NO_EXIST != index){
            double weight = (*ptr_w_vector)[index];
//        cost = exp(weight);
            cost = weight;
        } else{
            // std::cout << "index error"<<std::endl;
        }
        ptrnode->SetCost(cost);
    }
private:
    std::vector<Node> *ptr_start_node_;
    std::vector<Node> *ptr_stop_node_;
    std::vector<std::vector<std::vector<Node *>>> node_matrix_;
    //data related
    std::vector<double> *ptr_Z_;
};

#endif //CRF_CRFTHREAD_H
