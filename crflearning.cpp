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
    ptr_x_corpus_map_reverse_ = new std::map<int, std::string>;
    ptr_tag_map_ = new std::map<std::string, int>;
    ptr_tag_map_reverse_ = new std::map<int, std::string>;
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
        ptr_x_corpus_map_reverse_->insert(std::make_pair(index,(*it)));
        ptr_x_corpus_->push_back((*it));
        index++;
    }
    index = 0;
    for(std::set<std::string>::iterator it = ptr_tag_set_->begin(); it!= ptr_tag_set_->end(); ++it){
        ptr_tag_map_->insert(std::make_pair((*it),index));
        ptr_tag_map_reverse_->insert(std::make_pair(index,(*it)));
        index++;
    }
    ptr_feature_ = new Feature(ptr_x_vector_,ptr_tag_vector_,ptr_x_corpus_map_,ptr_tag_map_reverse_);
    tag_num_ = ptr_tag_set_->size();
    ptr_feature_vector_ = new std::vector<std::pair<int, int>>;
    ptr_empirical_e_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_empirical_e_->begin(),ptr_empirical_e_->end(),0);

    ptr_feature_bit_vector_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_feature_bit_vector_->begin(),ptr_feature_bit_vector_->end(),0);


    ptr_e_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_e_->begin(),ptr_e_->end(),0);
    ptr_gradient_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_gradient_->begin(),ptr_gradient_->end(),0);
    ptr_decoded_tag_ = new std::vector<std::string>(seq.size());

    CalcEmpiricalFi(seq);

}

CRFLearning::~CRFLearning() {
    delete ptr_x_corpus_;
    delete ptr_x_corpus_map_;
    delete ptr_tag_map_;
    delete ptr_feature_vector_;
    delete ptr_empirical_e_;
    delete ptr_e_;
    delete ptr_gradient_;
    delete ptr_decoded_tag_;
    DeleteLattice(*ptr_x_vector_);
}

void CRFLearning::DeleteLattice(std::vector<std::string> seq) {
    delete ptr_start_node_;
    delete ptr_stop_node_;
    for(int i=0; i<seq.size(); ++i){
        for(int j=0; j<tag_num_; ++j){
            delete node_matrix_[i][j];
        }
    }
    for(int i=0; i<seq.size(); ++i){
        std::vector<Node *> ptr_node_vector = node_matrix_[i];
        delete &ptr_node_vector;
    }
}

void CRFLearning::ResetParameters() {
    std::fill(ptr_e_->begin(),ptr_e_->end(),0.0);
}

void CRFLearning::BuildLattice(std::vector<std::string> seq) {
    BuildNode(seq);
    BuildLPath(seq);
    BuildRPath(seq);
}

void CRFLearning::BuildNode(std::vector<std::string> seq) {
    ptr_start_node_ = new Node(START_NODE_FLAG, START_NODE_FLAG);
    ptr_start_node_->SetNodeID(START_NODE_ID);
    ptr_stop_node_  = new Node(seq.size(),STOP_NODE_FLAG);
    ptr_stop_node_->SetNodeID(STOP_NODE_ID);
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
            Node *pnode = new Node(observ + FEATURE_CODE_OFFSET, tag);
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

void CRFLearning::SetPathFeature(std::pair<int, int> feature_pair, Path *ppath) {
        int feature_index = FEATURE_NO_EXIST;
        if(ptr_feature_->GetFeatureMap()->find(feature_pair) != ptr_feature_->GetFeatureMap()->end()){
            feature_index = ptr_feature_->GetFeatureMap()->find(feature_pair)->second;
        }
        ppath->SetFeatureIndex(feature_index);
}

//build left path
void CRFLearning::BuildLPath(std::vector<std::string> seq) {
    //create the vector of left path for each node
    for (int i = 0; i <seq.size() ; ++i) {
        for(int j=0; j<tag_num_;++j){
            if(i==0){
                Path *ppath = new Path(ptr_start_node_,node_matrix_[i][j]);
                //set feature index for each path (edge)
                SetPathFeature(std::make_pair(START_NODE_FLAG,node_matrix_[i][j]->GetY()),ppath);
                node_matrix_[i][j]->AddPath(ppath, true);
            } else{
                for(int k=0; k<tag_num_; ++k){
                    Path *ppath = new Path(node_matrix_[i-1][k],node_matrix_[i][j]);
                    //set feature index for each path (edge)
                    SetPathFeature(std::make_pair(node_matrix_[i-1][k]->GetY(),node_matrix_[i][j]->GetY()),ppath);
                    node_matrix_[i][j]->AddPath(ppath, true);
                }
            }
        }
    }
    //create the left path vector for the stop node;
    for(int k=0; k<tag_num_; ++k){
        Path *ppath = new Path(node_matrix_[seq.size()-1][k],ptr_stop_node_);
        SetPathFeature(std::make_pair(node_matrix_[seq.size()-1][k]->GetY(),STOP_NODE_FLAG),ppath);
        ptr_stop_node_->AddPath(ppath, true);
    }
}

//build right path
void CRFLearning::BuildRPath(std::vector<std::string> seq) {
    for (int i = seq.size()-1; i >=0 ; --i) {
        for(int j=0; j<tag_num_; ++j){
            if(i == seq.size()-1){
                Path *ppath = new Path(node_matrix_[i][j],ptr_stop_node_);
                SetPathFeature(std::make_pair(node_matrix_[i][j]->GetY(),STOP_NODE_FLAG),ppath);
                node_matrix_[i][j]->AddPath(ppath, false);
            }else{
                for(int k=0; k<tag_num_; ++k){
                    Path *ppath = new Path(node_matrix_[i][j], node_matrix_[i+1][k]);
                    SetPathFeature(std::make_pair(node_matrix_[i][j]->GetY(),node_matrix_[i+1][k]->GetY()),ppath);
                    node_matrix_[i][j]->AddPath(ppath,false);
                }
            }
        }
    }
    //create the right path vector for the start node;
    for(int k=0; k<tag_num_; ++k){
        Path *ppath = new Path(ptr_start_node_,node_matrix_[0][k]);
        SetPathFeature(std::make_pair(START_NODE_FLAG,node_matrix_[0][k]->GetY()),ppath);
        ptr_start_node_->AddPath(ppath, false);
    }
}

//calc E'(f_i)
double CRFLearning::CalcEmpiricalFi(std::vector<std::string> seq) {
    //count the feature number
    int tag_size = ptr_tag_vector_->size();
    int x_size = seq.size();
    //calc count of trans feature
    //start node;
    int index = ptr_feature_->GetFeatureIndex(std::make_pair(START_NODE_FLAG, ptr_tag_map_->find((*ptr_tag_vector_)[0])->second));
    (*ptr_empirical_e_)[index] = 1;
    (*ptr_feature_bit_vector_)[index] = 1;

    for(int i = 0; i < tag_size; ++i){
        int y_i = 0;
        int y_i_1 = ptr_tag_map_->find((*ptr_tag_vector_)[i])->second;
        if(i == tag_size-1){
             y_i = STOP_NODE_FLAG;
        } else{
             y_i = ptr_tag_map_->find((*ptr_tag_vector_)[i+1])->second;
        }
        index = ptr_feature_->GetFeatureIndex(std::make_pair(y_i_1,y_i));
        //count the num for each tag sequence
        (*ptr_empirical_e_)[index] += 1;
        (*ptr_feature_bit_vector_)[index] += 1;
#ifdef DEBUG_MODE
        std::cout << "the value of tag "<<y_i_1 << ", "<<y_i << " is: "<< (*ptr_empirical_e_)[index] <<",index is:"<<index<<std::endl;
#endif
    }
    //calc count of emission feature
    for(int i=0; i<x_size; ++i){
            int x = ptr_x_corpus_map_->find(seq[i])->second + FEATURE_CODE_OFFSET;
            int y = ptr_tag_map_->find((*ptr_tag_vector_)[i])->second;
            if(ptr_feature_->GetFeatureMap()->find(std::make_pair(x,y)) != ptr_feature_->GetFeatureMap()->end()){
                 index = ptr_feature_->GetFeatureIndex(std::make_pair(x,y));
                //count the num for each
                (*ptr_empirical_e_)[index] += 1;
                (*ptr_feature_bit_vector_)[index] += 1;
                //std::cout << "the observ and tag are "<< seq[i] <<"," << (*ptr_tag_vector_)[i] <<",index is:"<<index<<std::endl;
#ifdef DEBUG_MODE
                std::cout << "the value of observ and tag "<<x<< ", "<<y << " is: "<< (*ptr_empirical_e_)[index] <<std::endl;
#endif
            }
    }
    //stop node
    index = ptr_feature_->GetFeatureIndex(std::make_pair(ptr_tag_map_->find((*ptr_tag_vector_)[seq.size()-1])->second,STOP_NODE_FLAG));
    (*ptr_empirical_e_)[index] = 1;
    (*ptr_feature_bit_vector_)[index] = 1;
}

void CRFLearning::CalcFeatureExpectation(std::vector<std::string> seq) {
    int x_size = seq.size();
    //cal expectation of the start node's right paths.
    std::vector<Path *> ppath = ptr_start_node_->GetRPath();
    double expectation = 0;
    for(std::vector<Path *>::iterator it = ppath.begin();it!=ppath.end();++it) {
        (*it)->CalcExpectation(Z_, ptr_e_, &expectation);
    }
    //calc expectation for each node and its right path.
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcExpectation(Z_,ptr_e_);
#ifdef DEBUG_MODE
            std::cout << "expectation of the node"<<i<< ", "<<j<< " is: "<< node_matrix_[i][j]->GetExpectation() <<std::endl;
#endif
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
#ifdef DEBUG_MODE
    for(int i=0; i<ptr_feature_->GetFeatureSize(); ++i){
        std::cout << "feature expectation of the" <<i<<"th feature is"<<(*ptr_e_)[i]<<std::endl;
    }
#endif
}

//cost indicates w_k * \phi_k(s', s, x)
void CRFLearning::CalcCost(std::vector<std::string> seq) {
    ptr_start_node_->SetCost(1);
    ptr_stop_node_->SetCost(1);
    std::vector< Path*> ptr_start_rpath = ptr_start_node_->GetRPath();
    for(std::vector< Path*>::iterator it = ptr_start_rpath.begin(); it!= ptr_start_rpath.end(); ++it){
        ptr_feature_->CalcCost((*it));
    }
    std::vector<Path*> ptr_stop_lpath = ptr_stop_node_->GetLPath();
    for(std::vector< Path*>::iterator it = ptr_stop_lpath.begin(); it!= ptr_stop_lpath.end(); ++it){
        ptr_feature_->CalcCost((*it));
    }

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
#ifdef DEBUG_MODE
            if(i==seq.size()-1){
                std::cout << "this is the last row"<<std::endl;
            }
#endif
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
        double penalty = 0; //pre_w_k / (L2_FACTOR * L2_FACTOR);
        (*ptr_gradient_)[k] = empirical_e_k - e_k - penalty;
        std::cout << "the gradient of the "<<k<<"th feature is: "<<(*ptr_gradient_)[k]<<std::endl;
#ifdef DEBUG_MODE_
        std::cout << "the gradient of the "<<k<<"th feature is: "<<(*ptr_gradient_)[k]<<std::endl;
#endif
    }
}

void CRFLearning::UpdateWeight() {
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    int size = ptr_feature_->GetFeatureSize();
    for(int k=0; k<size; ++k){
        double weight = (*ptr_weight)[k] + LEARNING_RATE * (*ptr_gradient_)[k];
        ptr_feature_->SetWeightVector(k,weight);
    }
    /*
    double value =0;
    for(int k=0; k<size; ++k){
        value += ((*ptr_gradient_)[k] * (*ptr_gradient_)[k]);
    }
    value = sqrt(value);
    std::cout << "the gradient norm is: "<<value<<std::endl;*/
}

//calc the loss function.
double CRFLearning::CalcLoglikelihoodFunction(std::vector<std::string> seq) {
    int feature_size = ptr_feature_->GetFeatureSize();
    std::vector<double > *ptr_weight = ptr_feature_->GetWeightVector();
    double sum_numerator = 0; // (*ptr_weight)[2] + (*ptr_weight)[5] + (*ptr_weight)[6] + (*ptr_weight)[9] + (*ptr_weight)[10];
    for(int i=0; i<feature_size; ++i){
        sum_numerator += (*ptr_weight)[i] * (*ptr_feature_bit_vector_)[i];
    }
    //std::cout << (*ptr_weight)[2]<<","<<(*ptr_weight)[5]<<","<<(*ptr_weight)[6]<<","<<(*ptr_weight)[9]<<","<<(*ptr_weight)[10]<<std::endl;
    double sum_denominator = log(Z_);
    double l2 = 0;
    for (int i = 0; i < feature_size; ++i) {
        double w = (*(ptr_feature_->GetWeightVector()))[i];
        l2 += (w * w) / (2 * L2_FACTOR * L2_FACTOR);
    }
    l2 = 0;
    return  sum_numerator - sum_denominator - l2;
}
   
//calc alpha and beta and Z for all nodes in the graph.
void CRFLearning::ForwardBackward(std::vector<std::string> seq) {
    int x_size = seq.size();
    ptr_start_node_->SetAlpha(1);
    ptr_stop_node_->SetBeta(1);
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcAlpha();
#ifdef DEBUG_MODE
            std::cout << "the alpha of node "<<i<<","<<j<<" is: "<<node_matrix_[i][j]->GetAlpha()<<std::endl;
#endif
        }
    }
    for (int i = x_size-1; i>=0 ; --i) {
        for(int j=0; j<tag_num_; ++j){
            node_matrix_[i][j]->CalcBeta();
#ifdef DEBUG_MODE
            std::cout << "the beta of node "<<i<<","<<j<<" is: "<<node_matrix_[i][j]->GetBeta()<<std::endl;
#endif
        }
    }
    Z_ = 0;
    std::vector<Path *> stop_lpath= ptr_stop_node_->GetLPath();
    for (int j = 0; j<tag_num_; ++j) {
        Node *p_node = stop_lpath[j]->GetLNode();
        Z_ = p_node->SumExp(Z_,stop_lpath[j]->GetCost()*p_node->GetCost(),p_node->GetAlpha(),(j==0));
    }
    //std::cout << "Z_ derived by alpha is"<<Z_<<std::endl;

#ifdef DEBUG_MODE_
    std::cout << "Z_ derived by alpha is"<<Z_<<std::endl;
#endif
    std::vector<Path *> start_rpath= ptr_start_node_->GetRPath();
    double beata_Z = 0;
    for (int j = 0; j < tag_num_; ++j) {
        Node *p_node = start_rpath[j]->GetRNode();
        beata_Z = p_node->SumExp(beata_Z, ptr_start_node_->GetCost() * start_rpath[j]->GetCost(),p_node->GetBeta(),(j==0));
    }
    //std::cout<<"Z by beta is: "<< beata_Z <<std::endl;

#ifdef DEBUG_MODE_
    std::cout << "Z_ derived by beta is"<<beata_Z<<std::endl;
    std::cout << "==========="<<std::endl;
#endif
    for(int i=0; i<seq.size(); ++i){
        double value = 0;
        for(int k=0; k<tag_num_; k++){
            value += node_matrix_[i][k]->GetAlpha() * node_matrix_[i][k]->GetBeta();
        }
        //std::cout << i<<"th Z is: "<< value <<std::endl;
        Z_ = value;
    }

}

void CRFLearning::Viterbi(std::vector<std::string> seq) {
    for(int i=0; i<seq.size(); ++i){
        for(int j=0; j<tag_num_; ++j){
            double best_cost = 1e10;
            Node *p_best_node;
            std::vector<Path *> ppath = node_matrix_[i][j]->GetLPath();
            for(std::vector<Path *>::iterator it = ppath.begin(); it!=ppath.end(); ++it){
                Node *pnode = (*it)->GetLNode();
                double cost = pnode->GetBestCost() * ((*it)->GetCost() + node_matrix_[i][j]->GetCost());
                if(cost > best_cost){
                    best_cost = cost;
                    p_best_node = pnode;
                }
            }
            node_matrix_[i][j]->SetBestCost(best_cost);
            node_matrix_[i][j]->SetPreNode(p_best_node);
        }
    }
    //for the last node;
    double best_cost = 1e10;
    Node *p_best_node;
    for(int j=0; j<tag_num_; ++j){
       double cost = node_matrix_[seq.size()-1][j]->GetBestCost();
       if(cost > best_cost){
            best_cost = cost;
            p_best_node = node_matrix_[seq.size()-1][j];
       }
    }
    ptr_stop_node_->SetPreNode(p_best_node);
    //backtracking
    ViterbiBackTracking(seq);
}

void CRFLearning::ViterbiBackTracking(std::vector<std::string> seq) {
    Node *ptr_node = ptr_stop_node_->GetPreNode();
    std::string str = ptr_tag_map_reverse_->find(ptr_node->GetY())->second;
    std::cout << "tag inferred by viterbi is:"<<str<<std::endl;
    (*ptr_decoded_tag_)[seq.size() - 1] = str;
    for (int i = seq.size() - 2; i >= 0; --i) {
        ptr_node = ptr_node->GetPreNode();
        str = ptr_tag_map_reverse_->find(ptr_node->GetY())->second;
        (*ptr_decoded_tag_)[i] = str;
        std::cout << "tag inferred by viterbi is:"<<str<<std::endl;
    }
}

void CRFLearning::PrintPath(Node *pNode) {
    if(pNode->GetNodeID() == START_NODE_ID){
        return;
    }
    std::vector<Path* > pPath = pNode->GetLPath();
    for(std::vector<Path *>::iterator it = pPath.begin(); it != pPath.end(); ++it){
        std::string lstr = ptr_tag_map_reverse_->find((*it)->GetLNode()->GetY())->second;
        if((*it)->GetLNode()->GetNodeID()==START_NODE_ID){
            lstr = "START_NODE";
        }
        std::string rstr;
        if(pNode->GetNodeID() == STOP_NODE_ID){
            rstr = "STOP_NODE";
        }else{
            rstr = ptr_tag_map_reverse_->find(pNode->GetY())->second;
        }
        std::cout <<"the feature(t) id is"<<(*it)->GetFeatureIndex()<<" and pair is "<<lstr<<","<<rstr<<std::endl;
        //std::cout <<"the feature(e) id is"<<(*it)->GetLNode()->GetFeatureIndex()<<std::endl;
        PrintPath((*it)->GetLNode());
    }
}

void CRFLearning::Learning() {
    std::vector<std::string> seq = *ptr_x_vector_;
    if(!is_initialized_){
        Init(seq);
        is_initialized_ = true;
    }
    BuildLattice(seq);
    //for test only;
#ifdef DEBUG_MODE
    PrintPath(ptr_stop_node_);
#endif
    //training
    while(!is_converged_){
        ResetParameters();
        CalcGradient(seq);
        UpdateWeight();
        //calc loss function
        double loss_value = CalcLoglikelihoodFunction(seq);
        std::cout << "loglikelihood is:"<<loss_value<<std::endl;
        if(std::abs(loss_value) - std::abs(loss_value_) > CONVERGED_VALUE){
            loss_value_ = loss_value;
        } else{
            std::cout << "Training completed"<<std::endl;
            is_converged_ = true;
        }
    }
    Viterbi(seq);
}

