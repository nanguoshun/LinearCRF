//
// Created by  ngs on 28/06/2018.
//
#include "crfthread.h"

CRFThread::CRFThread(std::vector<Node> *ptr_start_node, std::vector<Node> *ptr_stop_node,
                     std::vector<std::vector<std::vector<Node *>>> node_matrix, std::vector<double> *ptr_Z) {

    ptr_start_node_ = ptr_start_node;
    ptr_stop_node_ = ptr_stop_node;
    node_matrix_ = node_matrix;
    ptr_Z_ = ptr_Z;
}

CRFThread::~CRFThread() {

}

void CRFThread::ForwardBackward(int x_size, int y_size, int seq_no) {
    (*ptr_start_node_)[seq_no].SetAlpha(0);
    (*ptr_stop_node_)[seq_no].SetBeta(0);
    int start_time, end_time;
    start_time = clock();
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            node_matrix_[seq_no][i][j]->CalcAlpha();
#ifdef DEBUG_MODE
            std::cout << "the alpha of node "<<i<<","<<j<<" is: "<<node_matrix_[i][j]->GetAlpha()<<std::endl;
#endif
        }
    }
    end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME_
    std::cout << "forward time is: "<< end_time - start_time <<std::endl;
#endif
    start_time = clock();
    for (int i = x_size - 1; i >= 0; --i) {
        for (int j = 0; j < y_size; ++j) {
            node_matrix_[seq_no][i][j]->CalcBeta();
#ifdef DEBUG_MODE
            std::cout << "the beta of node "<<i<<","<<j<<" is: "<<node_matrix_[i][j]->GetBeta()<<std::endl;
#endif
        }
    }
    end_time = clock();

#ifdef PERFORMANCE_CLOCK_TIME_
    std::cout << "backward time is: "<< end_time - start_time <<std::endl;
#endif
    start_time = clock();
    (*ptr_Z_)[seq_no] = 0;
    std::vector < Path * > stop_lpath = (*ptr_stop_node_)[seq_no].GetLPath();
    for (int j = 0; j < y_size; ++j) {
        Node *p_node = stop_lpath[j]->GetLNode();
//        (*ptr_Z_)[seq_no] = p_node->SumExp((*ptr_Z_)[seq_no],stop_lpath[j]->GetCost()*p_node->GetCost(),p_node->GetAlpha(),(j==0));
        (*ptr_Z_)[seq_no] = p_node->LogSumExp((*ptr_Z_)[seq_no], stop_lpath[j]->GetCost() + p_node->GetCost(),
                                              p_node->GetAlpha(), (j == 0));
    }
    end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME_
    std::cout << "Calc Z time is: "<< end_time - start_time <<std::endl;
#endif
#ifdef DEBUG_MODE_
    std::cout<<"the value of pnode is: "<<(*ptr_Z_)[seq_no]<<std::endl;
    std::vector<Path *> start_rpath= (*ptr_start_node_)[seq_no].GetRPath();
    double beata_Z = 0;
    for (int j = 0; j < (*ptr_tag_set_matrix_)[seq_no].size(); ++j) {
        Node *p_node = start_rpath[j]->GetRNode();
        beata_Z = p_node->LogSumExp(beata_Z, (*ptr_start_node_)[seq_no].GetCost() + start_rpath[j]->GetCost(),p_node->GetBeta(),(j==0));
    }
    std::cout << "Z_ derived by alpha is:" << beata_Z <<std::endl;
#endif
#ifdef DEBUG_MODE_
    std::cout << "Z_ derived by beta is"<<beata_Z<<std::endl;
    std::cout << "==========="<<std::endl;
    for(int i=0; i<seq.size(); ++i){
        double value = 0;
        for(int k=0; k<(*ptr_tag_set_matrix_)[seq_no].size(); k++){
            value = LogSumExp(value,node_matrix_[seq_no][i][k]->GetAlpha() + node_matrix_[seq_no][i][k]->GetBeta(),(value==0));
        }
        std::cout << i<<"th Z is: "<< value <<std::endl;
    }
#endif
}

//cost indicates w_k * \phi_k(s', s, x)
void CRFThread::CalcCost(int x_size, int y_size, int seq_no,std::vector<double> *ptr_w_weight) {
    (*ptr_start_node_)[seq_no].SetCost(DEFAULT_COST_VALUE);
    (*ptr_stop_node_)[seq_no].SetCost(DEFAULT_COST_VALUE);
    std::vector< Path*> ptr_start_rpath = (*ptr_start_node_)[seq_no].GetRPath();
    for(std::vector< Path*>::iterator it = ptr_start_rpath.begin(); it!= ptr_start_rpath.end(); ++it){
        CalcCost((*it),ptr_w_weight);
    }
    std::vector<Path*> ptr_stop_lpath = (*ptr_stop_node_)[seq_no].GetLPath();
    for(std::vector< Path*>::iterator it = ptr_stop_lpath.begin(); it!= ptr_stop_lpath.end(); ++it){
        CalcCost((*it),ptr_w_weight);
    }
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<y_size; ++j){
            //calc node cost
            CalcCost(node_matrix_[seq_no][i][j],ptr_w_weight);
            //calc left path cost of each node
            std::vector<Path *> lpath = node_matrix_[seq_no][i][j]->GetLPath();
            std::vector<Path *> rpath = node_matrix_[seq_no][i][j]->GetRPath();
            for(std::vector<Path *>::iterator it = lpath.begin(); it!=lpath.end(); ++it){
                CalcCost((*it),ptr_w_weight);
            }
#ifdef DEBUG_MODE
            if(i==seq.size()-1){
                std::cout << "this is the last row"<<std::endl;
            }
#endif
            for(std::vector<Path *>::iterator it = rpath.begin(); it!=rpath.end(); ++it){
                CalcCost((*it),ptr_w_weight);
            }
        }
    }
}

void CRFThread::CRFThreadRun(int x_size, int y_size, int seq_no, std::vector<double> *ptr_w_vector) {
    CalcCost(x_size,y_size,seq_no,ptr_w_vector);
    ForwardBackward(x_size,y_size,seq_no);
}
