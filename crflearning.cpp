//
// Created by  ngs on 09/06/2018.
//

#include "crflearning.h"

CRFLearning::CRFLearning(DatasetMgr *ptr_datamgr) {
    ptr_datamgr_ = ptr_datamgr;
    loss_value_ = 0;
    ptr_x_corpus_ = new std::vector<std::string>;
    ptr_x_corpus_map_ = new std::map<std::string, int>;
    ptr_x_set_ = ptr_datamgr_->GetTrainingXSet();
    ptr_tag_vector_ = ptr_datamgr_->GetTageVector();
    ptr_x_vector_  = ptr_datamgr_->GetTrainingXVector();
    ptr_tag_set_ = ptr_datamgr_->GetTagSet();
}

void CRFLearning::Init() {
    //to simplify the learning, we use the training x set as corpus.
    int index = 0;
    for (std::set<std::string>::iterator it = ptr_x_set_->begin(); it != ptr_x_set_->end(); ++it) {
        std::cout << (*it) << std::endl;
        ptr_x_corpus_map_->insert(std::make_pair((*it), index));
        ptr_x_corpus_->push_back((*it));
        index++;
    }
    index = 0
    for(std::set<std::string>::iterator it = ptr_tag_set_->begin(); it!= ptr_tag_set_->end(); ++it){
        ptr_tag_map_->insert(std::make_pair((*it),index));
        index++;
    }
    ptr_feature_ = new Feature(ptr_x_vector_,ptr_tag_vector_,ptr_x_corpus_map_,ptr_tag_map_);

}

CRFLearning::~CRFLearning() {
    delete ptr_x_corpus_;
    delete ptr_x_corpus_map_;

}

void CRFLearning::DeleteLattice(std::vector<std::string> seq) {
    for(int i=0; i<seq.size(); ++i){
        for(int j=0; j<tag_num_; ++j){
        }
    }
}

void CRFLearning::BuildLattice(std::vector<std::string> seq) {
    //create each row
    for(int i=0; i<seq.size(); ++i){
        std::vector<Node *> *ptr_node_vector = new std::vector<Node *>;
        node_matrix_.push_back(*ptr_node_vector);
    }
    //create node in each row
    int tag_num = ptr_tag_set_->size();
    for(int i=0; i<seq.size(); ++i){
        for(int j=0; j<tag_num; ++j){
            Node *pnode = new Node(i, j);
            node_matrix_[i].push_back(pnode);
        }
    }
    //create path for each node
    for (int i = 0; i <seq.size() ; ++i) {
        for(int j=0; j<tag_num;++j){
           if(i==0){
               Path *ppath = new Path(ptr_start_node_,node_matrix_[i][j]);
               node_matrix_[i][j]->AddPath(ppath, true);
           } else{
               for(int k=0; k<tag_num; ++k){
                   Path *ppath = new Path(node_matrix_[i-1][k],node_matrix_[i][j]);
                   node_matrix_[i][j]->AddPath(ppath, true);
               }
           }
        }
    }
}

//
void CRFLearning::Learning() {
    std::vector<std::string> seq = *ptr_x_vector_;
    BuildLattice(seq);
    //define loss function
    double loss_value = CalcLossFunction(seq);
    if(loss_value - loss_value_ > CONVERGED_VALUE){
        loss_value_ = loss_value;
    } else{
        std::cout << "Training completed"<<std::endl;
        return;
    }
    //training
    CalcGradient();
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    int size = ptr_feature_->GetFeatureSize();
    for(int k=0; k<size; ++k){
        (*ptr_weight)[k] = (*ptr_weight)[k] - LEARNING_RATE * (*ptr_gradient_)[k];
    }
}

double CRFLearning::CalcEmpericialFi() {
    //count the feature number

}

double CRFLearning::CalcExpectedFi(std::vector<std::string> seq) {
    int x_size = seq.size();
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcExpectation(Z_);
        }
    }
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<x_size; ++j){

        }
    }

}

// stochastic gradient descent;
void CRFLearning::CalcGradient() {
    int size = ptr_feature_->GetFeatureSize();
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    for(int k=0; k<size; ++k){
        double empirical_e_k = ptr_feature_->GetEmpiricalFeatureExpectation(k);
        double e_k = ptr_feature_->GetFeatureExpectation(k);
        double pre_w_k = (*ptr_weight)[k];
        double penalty = pre_w_k / L2_FACTOR;
        (*ptr_gradient_)[k] = empirical_e_k - e_k -penalty;
    }
}

double CRFLearning::CalcLossFunction(std::vector<std::string> seq) {
   std::vector<double> *p_w_vector = ptr_feature_->GetWeightVector();
   int feature_size = p_w_vector->size();
   double value = 0;
   for(int j = 1; j< seq.size(); ++j){
       int y_j_1 = ptr_tag_map_->find((*ptr_tag_vector_)[j-1])->second;
       int y_j = ptr_tag_map_->find((*ptr_tag_vector_)[j])->second;
       int index  = ptr_feature_->GetFeatureIndex(std::make_pair(y_j_1,y_j));
       value += (*p_w_vector)[index];
     }
    for(int j = 0; j< seq.size(); ++j) {
       int y_j = ptr_tag_map_->find((*ptr_tag_vector_)[j])->second;
       int x = ptr_x_corpus_map_->find(seq[j])->second;
       int index = ptr_feature_->GetFeatureIndex(std::make_pair(y_j,x));
       value += (*p_w_vector)[index];
    }
    return  value;
}
   

void CRFLearning::ForwardBackward(std::vector<std::string> seq) {
    int x_size = seq.size();
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcAlpha();
        }
    }
    for (int i = x_size-1; i>=0 ; --i) {
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcBeta();
        }
    }
    Z_ = 0;
    std::vector<Path *> stop_lpath_ = ptr_stop_node_->GetLPath();
    for (int j = 0; j<tag_num_; ++j) {
        Node *p_node = stop_lpath_[j]->GetLNode();
        double pathcost = stop_lpath_[j]->GetCost();
        Z_ = p_node->SumExp(Z_,pathcost,p_node->GetAlpha(),(j==0));
    }
}


