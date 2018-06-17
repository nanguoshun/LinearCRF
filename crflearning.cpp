//
// Created by  ngs on 09/06/2018.
//

#include "crflearning.h"
#include <numeric>
#include <cmath>

CRFLearning::CRFLearning(DatasetMgr *ptr_datamgr) {
    ptr_datamgr_ = ptr_datamgr;
    loss_value_ = -1000;
    ptr_x_corpus_ = new std::vector<std::string>;
    ptr_x_corpus_map_ = new std::map<std::string, int>;
    ptr_tag_map_ = new std::map<std::string, int>;
    ptr_x_set_ = ptr_datamgr_->GetTrainingXSet();
    ptr_tag_vector_ = ptr_datamgr_->GetTageVector();
    ptr_x_vector_  = ptr_datamgr_->GetTrainingXVector();
    ptr_tag_set_ = ptr_datamgr_->GetTagSet();
    is_initialized_ = false;
    is_converged_  = false;
}

void CRFLearning::Init(std::vector<std::string> seq) {
    //to simplify the learning, we use the training x set as corpus.
    int index = 0;
    for (std::set<std::string>::iterator it = ptr_x_set_->begin(); it != ptr_x_set_->end(); ++it) {
        std::cout << (*it) << std::endl;
        ptr_x_corpus_map_->insert(std::make_pair((*it), index));
        ptr_x_corpus_->push_back((*it));
        index++;
    }
    index = 0;
    for(std::set<std::string>::iterator it = ptr_tag_set_->begin(); it!= ptr_tag_set_->end(); ++it){
        ptr_tag_map_->insert(std::make_pair((*it),index));
        index++;
    }
    ptr_feature_ = new Feature(ptr_x_vector_,ptr_tag_vector_,ptr_x_corpus_map_,ptr_tag_map_);
    tag_num_ = ptr_tag_set_->size();
    ptr_feature_vector_ = new std::vector<std::pair<int, int>>;
    ptr_empirical_e_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_empirical_e_->begin(),ptr_empirical_e_->end(),0);
    CalcEmpiricalFi(seq);
    ptr_e_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_e_->begin(),ptr_e_->end(),0);
    ptr_gradient_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_gradient_->begin(),ptr_gradient_->end(),0);
}

CRFLearning::~CRFLearning() {
    delete ptr_x_corpus_;
    delete ptr_x_corpus_map_;
    delete ptr_tag_map_;
    delete ptr_feature_vector_;
    delete ptr_empirical_e_;
    delete ptr_e_;
    DeleteLattice(*ptr_x_vector_);
    //to be completed....
}

void CRFLearning::DeleteLattice(std::vector<std::string> seq) {
    for(int i=0; i<seq.size(); ++i){
        for(int j=0; j<tag_num_; ++j){
        }
    }
}

void CRFLearning::BuildLattice(std::vector<std::string> seq) {
    BuildNode(seq);
    BuildLPath(seq);
    BuildRPath(seq);
}

void CRFLearning::BuildNode(std::vector<std::string> seq) {
    ptr_start_node_ = new Node(START_NODE_FLAG, START_NODE_FLAG);
    ptr_stop_node_  = new Node(seq.size(),STOP_NODE_FLAG);
    //create each row
    for(int i=0; i<seq.size(); ++i){
        std::vector<Node *> *ptr_node_vector = new std::vector<Node *>;
        node_matrix_.push_back(*ptr_node_vector);
    }
    //create node in each row
    for(int i=0; i<seq.size(); ++i){
        for(std::set<std::string>::iterator it = ptr_tag_set_->begin(); it!=ptr_tag_set_->end(); ++it){
            std::cout << "the string value of x and y are: "<<(*ptr_x_vector_)[i]<<", "<<(*it)<<std::endl;
            int observ = ptr_x_corpus_map_->find((*ptr_x_vector_)[i])->second;
            int tag = ptr_tag_map_->find((*it))->second;
            Node *pnode = new Node(observ, tag);
            //set feature index for each node.
            int feature_index = FEATURE_NO_EXIST;
            if(ptr_feature_->GetFeatureMap()->find(std::make_pair((observ+FEATURE_CODE_OFFSET),tag)) != ptr_feature_->GetFeatureMap()->end()){
                feature_index = ptr_feature_->GetFeatureMap()->find(std::make_pair((observ+FEATURE_CODE_OFFSET),tag))->second;
            }
            pnode->SetFeatureIndex(feature_index);
            node_matrix_[i].push_back(pnode);
        }
    }
}
//build left path
void CRFLearning::BuildLPath(std::vector<std::string> seq) {
    //create the vector of left path for each node
    for (int i = 0; i <seq.size() ; ++i) {
        for(int j=0; j<tag_num_;++j){
            if(i==0){
                Path *ppath = new Path(ptr_start_node_,node_matrix_[i][j]);
                node_matrix_[i][j]->AddPath(ppath, true);
            } else{
                for(int k=0; k<tag_num_; ++k){
                    Path *ppath = new Path(node_matrix_[i-1][k],node_matrix_[i][j]);
                    //add left path vector
                    int tag1 = node_matrix_[i-1][k]->GetY();
                    int tag2 = node_matrix_[i][j]->GetY();
                    //set feature index for each path (edge)
                    int feature_index = FEATURE_NO_EXIST;
                    if(ptr_feature_->GetFeatureMap()->find(std::make_pair(tag1,tag2)) != ptr_feature_->GetFeatureMap()->end()){
                        feature_index = ptr_feature_->GetFeatureMap()->find(std::make_pair(tag1,tag2))->second;
                    }
                    ppath->SetFeatureIndex(feature_index);
                    node_matrix_[i][j]->AddPath(ppath, true);
                }
            }
        }
    }
    //create the left path vector for the stop node;
    for(int k=0; k<tag_num_; ++k){
        Path *ppath = new Path(node_matrix_[seq.size()-1][k],ptr_stop_node_);
        ptr_stop_node_->AddPath(ppath, true);
    }
}

//build right path
void CRFLearning::BuildRPath(std::vector<std::string> seq) {
    for (int i = seq.size()-1; i >=0 ; --i) {
        for(int j=0; j<tag_num_; ++j){
            if(i == seq.size()-1){
                Path *ppath = new Path(node_matrix_[i][j],ptr_stop_node_);
                node_matrix_[i][j]->AddPath(ppath, false);
            }else{
                for(int k=0; k<tag_num_; ++k){
                    Path *ppath = new Path(node_matrix_[i][j], node_matrix_[i+1][k]);
                    node_matrix_[i][j]->AddPath(ppath,false);
                }
            }
        }
    }
    //create the right path vector for the start node;
    for(int k=0; k<tag_num_; ++k){
        Path *ppath = new Path(ptr_start_node_,node_matrix_[0][k]);
        ptr_start_node_->AddPath(ppath, false);
    }
}

void CRFLearning::UpdateWeight() {
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    int size = ptr_feature_->GetFeatureSize();
    for(int k=0; k<size; ++k){
        double weight = (*ptr_weight)[k] - LEARNING_RATE * (*ptr_gradient_)[k];
        ptr_feature_->SetWeightVector(k,weight);
    }
}

//calc E'(f_i)
double CRFLearning::CalcEmpiricalFi(std::vector<std::string> seq) {
    //count the feature number
    int tag_size = ptr_tag_vector_->size();
    int x_size = seq.size();
    //calc count of trans feature
    for(int i=0; i<tag_size-1; ++i){
        int y_i_1 = ptr_tag_map_->find((*ptr_tag_vector_)[i])->second;
        int y_i =   ptr_tag_map_->find((*ptr_tag_vector_)[i+1])->second;
        int index = ptr_feature_->GetFeatureMap()->find(std::make_pair(y_i_1,y_i))->second;
        //count the num for each tag sequence
        (*ptr_empirical_e_)[index] += 1;
        std::cout << "the value of tag "<<y_i_1 << ", "<<y_i << " is: "<< (*ptr_empirical_e_)[index] <<",index is:"<<index<<std::endl;
    }
    //calc count of emission feature
    for(int i=0; i<x_size; ++i){
            int x = ptr_x_corpus_map_->find(seq[i])->second + FEATURE_CODE_OFFSET;
            int y = ptr_tag_map_->find((*ptr_tag_vector_)[i])->second;
            if(ptr_feature_->GetFeatureMap()->find(std::make_pair(x,y)) != ptr_feature_->GetFeatureMap()->end()){
                int index = ptr_feature_->GetFeatureMap()->find(std::make_pair(x,y))->second;
                //count the num for each
                (*ptr_empirical_e_)[index] += 1;
                //std::cout << "the observ and tag are "<< seq[i] <<"," << (*ptr_tag_vector_)[i] <<",index is:"<<index<<std::endl;
                std::cout << "the value of observ and tag "<<x<< ", "<<y << " is: "<< (*ptr_empirical_e_)[index] <<std::endl;
            }
    }
}

void CRFLearning::CalcFeatureExpectation(std::vector<std::string> seq) {
    int x_size = seq.size();
    //calc expectation for each node and its left path.
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcExpectation(Z_,ptr_e_);
            std::cout << "expectation of the node"<<i<< ", "<<j<< " is: "<< node_matrix_[i][j]->GetExpectation() <<std::endl;
            //std::vector<Path *> lpath = node_matrix_[i][j]->GetLPath();
            //for (std::vector<Path *>::iterator it = lpath.begin();  it!=lpath.end() ; ++it) {
             //   }
        }
    }
    //summation
    /*
    int feature_size = ptr_feature_->GetFeatureSize();
    for(int index = 0; index<feature_size; index++){
        std::pair<int,int > feature = ptr_feature_->GetReverseFeatureMap()->find(index)->second;
        if(index<ptr_feature_->GetFeatureSizeEdge()){
            CalcEdgeFeatureExpectation(index, feature);
        } else{
            CalcNodeFeatureExpectation(index, feature);
        }
    }
    */
    for(int i=0; i<ptr_feature_->GetFeatureSize(); ++i){
        std::cout << "feature expectation of the" <<i<<"th feature is"<<(*ptr_e_)[i]<<std::endl;
    }
}

void CRFLearning::CalcEdgeFeatureExpectation(int index, std::pair<int, int> feature) {
    //for trans feature
    int y_i_1 = feature.first;
    int y_i = feature.second;
    //summation over a tag instance.
    for(int i=0; i<ptr_tag_vector_->size()-1; ++i){
        int tag1 = ptr_tag_map_->find((*ptr_tag_vector_)[i])->second;
        int tag2 = ptr_tag_map_->find((*ptr_tag_vector_)[i+1])->second;
        //if it is exactly the edge that is bundled to the feature. here, an edge indicates two tag.
        if(y_i_1==tag1 && y_i == tag2){
            (*ptr_e_)[index] +=  node_matrix_[y_i_1][y_i]->GetExpectation();
            std::cout << "the expectation of edge "<<y_i_1<< ", "<<y_i << " is: "<< (*ptr_e_)[index] <<std::endl;
        }
    }
}

void CRFLearning::CalcNodeFeatureExpectation(int index, std::pair<int, int> feature) {
//for emission feature
    for (int i = 0; i < ptr_tag_vector_->size(); ++i) {
        std::cout << "the observ and tag is: " << (*ptr_x_vector_)[i] + SPERATOR_FLAG + (*ptr_tag_vector_)[i]
                  << std::endl;
    }
    int x = feature.first - FEATURE_CODE_OFFSET;
    int y = feature.second;
    for (int i = 0; i < ptr_x_vector_->size(); ++i) {
        int observ = ptr_x_corpus_map_->find((*ptr_x_vector_)[i])->second;
        int tag = ptr_tag_map_->find((*ptr_tag_vector_)[i])->second;
        if (observ == x && tag == y) {
            (*ptr_e_)[index] += node_matrix_[i][y]->GetExpectation();
            std::cout << "the expectation of oberv and tag " << x << ", " << y << " is: " << (*ptr_e_)[index]
                      << std::endl;
        }
    }
}

//cost indicates w_k * \phi_k(s', s, x)
void CRFLearning::CalcCost(std::vector<std::string> seq) {
    ptr_stop_node_->SetCost(exp(1));
    for(int i=0; i<seq.size(); ++i){
        for(int j=0; j<tag_num_; ++j){
            //calc node cost
            ptr_feature_->CalcCost(node_matrix_[i][j]);
            //calc left path cost of each node
            std::vector<Path *> lpath = node_matrix_[i][j]->GetLPath();
            std::vector<Path *> rpath = node_matrix_[i][j]->GetRPath();
            for(std::vector<Path *>::iterator it = lpath.begin(); it!=lpath.end(); ++it){
                ptr_feature_->CalcCost(*it);
            }
            if(i==seq.size()-1){
                std::cout << "this is the last row"<<std::endl;
            }
            for(std::vector<Path *>::iterator it = rpath.begin(); it!=rpath.end(); ++it){
                ptr_feature_->CalcCost(*it);
            }
        }
    }
}

// stochastic gradient descent;
void CRFLearning::CalcGradient(std::vector<std::string> seq) {
    CalcCost(seq);
    ForwardBackward(seq);
    CalcFeatureExpectation(seq);
    int size = ptr_feature_->GetFeatureSize();
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    for(int k=0; k<size; ++k){
        double empirical_e_k = (*ptr_empirical_e_)[k];
        double e_k = (*ptr_e_)[k];
        double pre_w_k = (*ptr_weight)[k];
        double penalty = pre_w_k / (L2_FACTOR * L2_FACTOR);
        (*ptr_gradient_)[k] = empirical_e_k - e_k -penalty;
        std::cout << "the gradient of the "<<k<<"th feature is: "<<(*ptr_gradient_)[k]<<std::endl;
    }
}

//calc the loss function.
double CRFLearning::CalcLoglikelihoodFunction(std::vector<std::string> seq) {
    int feature_size = ptr_feature_->GetFeatureSize();
    double sum_numerator = std::accumulate(ptr_feature_->GetWeightVector()->begin(), ptr_feature_->GetWeightVector()->end(), 0.0f);
    double sum_denominator = log(Z_);
    double l2 = 0;
    for (int i = 0; i < feature_size; i++) {
        double w = (*(ptr_feature_->GetWeightVector()))[i];
        l2 += (w * w) / (2 * L2_FACTOR * L2_FACTOR);
    }
    double loglikelihood = sum_numerator - sum_denominator - l2;
    return loglikelihood;
}
   
//calc alpha and beta and Z for all nodes in the graph.
void CRFLearning::ForwardBackward(std::vector<std::string> seq) {
    int x_size = seq.size();
    ptr_start_node_->SetAlpha(1);
    ptr_stop_node_->SetBeta(1);
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcAlpha();
            std::cout << "the alpha of node "<<i<<","<<j<<" is: "<<node_matrix_[i][j]->GetAlpha()<<std::endl;
        }
    }
    for (int i = x_size-1; i>=0 ; --i) {
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcBeta();
            std::cout << "the beta of node "<<i<<","<<j<<" is: "<<node_matrix_[i][j]->GetBeta()<<std::endl;
        }
    }
    Z_ = 0;
    std::vector<Path *> stop_lpath= ptr_stop_node_->GetLPath();
    for (int j = 0; j<tag_num_; ++j) {
        Node *p_node = stop_lpath[j]->GetLNode();
        Z_ = p_node->SumExp(Z_,ptr_stop_node_->GetCost(),p_node->GetAlpha(),(j==0));
    }

}

void CRFLearning::Viterbi() {

}

void CRFLearning::Learning() {
    std::vector<std::string> seq = *ptr_x_vector_;
    if(!is_initialized_){
        Init(seq);
        is_initialized_ = true;
    }
    BuildLattice(seq);
    //training
    while(!is_converged_){
        CalcGradient(seq);
        UpdateWeight();
        //calc loss function
        double loss_value = CalcLoglikelihoodFunction(seq);
        if(std::abs(loss_value) - std::abs(loss_value_) < CONVERGED_VALUE){
            loss_value_ = loss_value;
        } else{
            //std::cout << "Training completed"<<std::endl;
            //is_converged_ = true;
            //return;
        }
    }
}