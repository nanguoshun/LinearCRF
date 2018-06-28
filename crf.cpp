//
// Created by  ngs on 09/06/2018.
//

#include "crf.h"
#include <numeric>
#include <cmath>
#include <string>
#include <pthread.h>

LinearCRF::LinearCRF(DatasetMgr *ptr_datamgr_training) {
    AllocateSpace();
    //is_training_ = isTraining;
    ptr_datamgr_ = ptr_datamgr_training;
    //ptr_datamgr_test_ = ptr_datamgr_test;
    ptr_x_traing_set_ = ptr_datamgr_->GetTrainingXSet();
    //ptr_x_test_set_ = ptr_datamgr_test_->GetTrainingXSet();
    ptr_tag_vector_ = ptr_datamgr_->GetTageVector();
    ptr_x_train_vector_  = ptr_datamgr_->GetTrainingXVector();
    //ptr_x_test_vector_ = ptr_datamgr_test_->GetTrainingXVector();
    ptr_tag_set_ = ptr_datamgr_->GetTagSet();
    //ptr_test_tag_set_ = ptr_datamgr_test_->GetTagSet();
    GenerateSeqFromVector(ptr_x_train_vector_,ptr_seq_matrix_);
    GenerateSeqFromVector(ptr_tag_vector_,ptr_tag_seq_);
    num_of_instance_ = ptr_datamgr_->GetNumOfTrainingSeqs();
    //tag_num_ = ptr_tag_set_->size();
    CreateTagObservMap();
    Init();
}

void LinearCRF::AllocateSpace() {
    //ptr_x_corpus_ = new std::vector<std::string>;
    ptr_x_corpus_map_ = new std::map<std::string, int>;
    ptr_x_corpus_map_reverse_ = new std::map<int, std::string>;
    ptr_tag_map_ = new std::map<std::string, int>;
    ptr_tag_map_reverse_ = new std::map<int, std::string>;
    ptr_seq_matrix_ = new std::vector<std::vector<std::string>>;
    ptr_tag_seq_ = new std::vector<std::vector<std::string>>;
    ptr_start_node_ = new std::vector<Node>;
    ptr_stop_node_ = new std::vector<Node>;
}

void LinearCRF::Init() {
    //to simplify the learning, we use the training x set as corpus.
    ptr_feature_ = new Feature(ptr_x_train_vector_, ptr_tag_vector_, ptr_x_corpus_map_, ptr_tag_map_reverse_,
                               num_of_instance_);
    ptr_feature_vector_ = new std::vector<std::pair<int, int>>;
    ptr_empirical_e_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_empirical_e_->begin(),ptr_empirical_e_->end(),0);
    ptr_feature_bit_vector_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_feature_bit_vector_->begin(),ptr_feature_bit_vector_->end(),0);
    ptr_e_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_e_->begin(),ptr_e_->end(),0);
    ptr_gradient_ = new std::vector<double>(ptr_feature_->GetFeatureSize());
    std::fill(ptr_gradient_->begin(),ptr_gradient_->end(),0);
    ptr_decoded_tag_ = new std::vector<std::vector<std::string>>();
    ptr_tag_set_matrix_ = new std::vector<std::set<std::string>>;
    for(int seq_no = 0; seq_no<num_of_instance_; ++seq_no){
        std::vector<std::string> *pvector = new std::vector<std::string>((*ptr_seq_matrix_)[seq_no].size());
        ptr_decoded_tag_->push_back(*pvector);
        std::set<std::string> *pset = new std::set<std::string>;
        ptr_tag_set_matrix_->push_back(*pset);
    }
    FromVectorToSet();
    ptr_Z_ = new std::vector<double>(num_of_instance_);
    is_converged_  = false;
    num_of_thread_ =  std::thread::hardware_concurrency();
    ptr_thread_vector_ = new std::vector<std::thread>;

}

void LinearCRF::CreateTagObservMap() {
    std::ofstream ofs_x_map("x_map.txt");
    std::ofstream ofs_tag_map("tag_map.txt");
    int index = 0;
    for (std::set<std::string>::iterator it = ptr_x_traing_set_->begin(); it != ptr_x_traing_set_->end(); ++it) {
        //std::cout << (*it) << std::endl;
        ptr_x_corpus_map_->insert(std::make_pair((*it), index));
        ptr_x_corpus_map_reverse_->insert(std::make_pair(index,(*it)));
        ofs_x_map << std::to_string(index) + " " + (*it);
        ofs_x_map << std::endl;
        index++;
    }
    index = 0;
    for(std::set<std::string>::iterator it = ptr_tag_set_->begin(); it!= ptr_tag_set_->end(); ++it){
        ptr_tag_map_->insert(std::make_pair((*it),index));
        ptr_tag_map_reverse_->insert(std::make_pair(index,(*it)));
        ofs_tag_map << std::to_string(index) + " "+(*it);
        ofs_tag_map << std::endl;
        index++;
    }
}

LinearCRF::~LinearCRF() {
    //delete ptr_x_corpus_;
    delete ptr_x_corpus_map_;
    delete ptr_x_corpus_map_reverse_;
    delete ptr_tag_map_;
    delete ptr_tag_map_reverse_;
    delete ptr_feature_vector_;
    delete ptr_empirical_e_;
    delete ptr_e_;
    delete ptr_gradient_;
    delete ptr_decoded_tag_;
    delete ptr_seq_matrix_;
    delete ptr_tag_seq_;
    delete ptr_Z_;
    delete ptr_thread_vector_;
    DeleteLattice();
}

void LinearCRF::DeleteLattice() {
    delete ptr_start_node_;
    delete ptr_stop_node_;
    for(int seq_no = 0; seq_no<num_of_instance_; ++seq_no){
        for(int i=0; i<(*ptr_seq_matrix_)[seq_no].size(); ++i){
            for(int j=0; j<(*ptr_tag_set_matrix_)[seq_no].size(); ++j){
                delete node_matrix_[seq_no][i][j];
            }
        }
    }
    for(int seq_no = 0; seq_no<num_of_instance_; ++seq_no) {
        for (int i = 0; i < (*ptr_seq_matrix_)[seq_no].size(); ++i) {
            std::vector<Node *> ptr_node_vector = node_matrix_[seq_no][i];
            delete &ptr_node_vector;
        }
    }
    for(int seq_no = 0; seq_no<num_of_instance_; ++seq_no){
        std::set<std::string> *ptr_tagset = &((*ptr_tag_set_matrix_)[seq_no]);
        std::vector<std::string> *ptr_decoded_tag = &((*ptr_decoded_tag_)[seq_no]);
        delete ptr_tagset;
        delete ptr_decoded_tag;
    }
    delete ptr_tag_set_matrix_;
    delete ptr_decoded_tag_;

}

void LinearCRF::ResetParameters() {
    std::fill(ptr_e_->begin(),ptr_e_->end(),0.0);
}

void LinearCRF::BuildLattice() {
#ifdef DEBUG_MODE_
    int m=1;
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        std::cout << seq_no<<"th sequence ============"<<std::endl;
        for(int k=0; k<seq.size(); ++k){
            std::cout <<m<<"line: "<<seq[k]<<std::endl;
            m++;
        }
        m+=2;
    }
#endif
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        Node *p_startnode = new Node(START_NODE_FLAG-seq_no, START_NODE_FLAG-seq_no);
        ptr_start_node_->push_back(*p_startnode);
        Node *p_stopnode = new Node(STOP_NODE_FLAG-seq_no, STOP_NODE_FLAG-seq_no);
        ptr_stop_node_->push_back(*p_stopnode);
    }
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no){
        std::vector<std::vector<Node *>> *ptr_matrix = new std::vector<std::vector<Node *>>;
        node_matrix_.push_back(*ptr_matrix);
    }
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no){
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        BuildNode(seq, seq_no);
    }
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        BuildLPath(seq,seq_no);
    }
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        BuildRPath(seq,seq_no);
    }
}


void LinearCRF::FromVectorToSet() {
    for(int seq_no=0; seq_no < num_of_instance_; ++seq_no){
        std::set<std::string> *ptr_tag_set = &((*ptr_tag_set_matrix_)[seq_no]);
        std::vector<std::string> tag_vector = (*ptr_tag_seq_)[seq_no];
        for(int i = 0; i<tag_vector.size(); ++i){
            if(ptr_tag_set->find(tag_vector[i]) == ptr_tag_set->end()){
                ptr_tag_set->insert(tag_vector[i]);
            }
        }
    }
}

void LinearCRF::BuildNode(std::vector<std::string> seq, int seq_no) {
    //create each row
    for (int i = 0; i < seq.size(); ++i) {
        std::vector<Node *> *ptr_node_vector = new std::vector<Node *>;
        node_matrix_[seq_no].push_back(*ptr_node_vector);
    }
    //create node in each row
    for (int i = 0; i < seq.size(); ++i) {
        for (std::set<std::string>::iterator it = (*ptr_tag_set_matrix_)[seq_no].begin(); it != (*ptr_tag_set_matrix_)[seq_no].end(); ++it) {
//        for (std::set<std::string>::iterator it = ptr_tag_set_->begin(); it != ptr_tag_set_->end(); ++it) {
            std::cout << "the string value of x and y are: " << seq[i] << ", " << (*it) << std::endl;
            int observ = ptr_x_corpus_map_->find(seq[i])->second;
            int tag = ptr_tag_map_->find((*it))->second;
            //std::cout << "tag Y is"<<tag<<std::endl;
            Node *pnode = new Node(observ + FEATURE_CODE_OFFSET, tag);
            //set feature index for each node.
            int feature_index = FEATURE_NO_EXIST;
            if (ptr_feature_->GetFeatureMap()->find(std::make_pair((observ + FEATURE_CODE_OFFSET), tag)) !=
                ptr_feature_->GetFeatureMap()->end()) {
                feature_index = ptr_feature_->GetFeatureMap()->find(
                        std::make_pair((observ + FEATURE_CODE_OFFSET), tag))->second;
            }
            pnode->SetFeatureIndex(feature_index);
            node_matrix_[seq_no][i].push_back(pnode);
        }
    }
}

void LinearCRF::SetPathFeature(std::pair<int, int> feature_pair, Path *ppath) {
        int feature_index = FEATURE_NO_EXIST;
        if(ptr_feature_->GetFeatureMap()->find(feature_pair) != ptr_feature_->GetFeatureMap()->end()){
            feature_index = ptr_feature_->GetFeatureMap()->find(feature_pair)->second;
        }
        ppath->SetFeatureIndex(feature_index);
}

//build left path
void LinearCRF::BuildLPath(std::vector<std::string> seq, int seq_no) {
    //create the vector of left path for each node
    for (int i = 0; i <seq.size() ; ++i) {
        for(int j=0; j<(*ptr_tag_set_matrix_)[seq_no].size();++j){
            if(i==0){
                Path *ppath = new Path(&((*ptr_start_node_)[seq_no]),node_matrix_[seq_no][i][j]);
                //set feature index for each path (edge)
                SetPathFeature(std::make_pair(START_NODE_FLAG-seq_no,node_matrix_[seq_no][i][j]->GetY()),ppath);
                node_matrix_[seq_no][i][j]->AddPath(ppath, true);
            } else{
                for(int k=0; k<(*ptr_tag_set_matrix_)[seq_no].size(); ++k){
                    Path *ppath = new Path(node_matrix_[seq_no][i-1][k],node_matrix_[seq_no][i][j]);
                    //set feature index for each path (edge)
                    SetPathFeature(std::make_pair(node_matrix_[seq_no][i-1][k]->GetY(),node_matrix_[seq_no][i][j]->GetY()),ppath);
                    node_matrix_[seq_no][i][j]->AddPath(ppath, true);
                }
            }
        }
    }
    //create the left path vector for the stop node;
    for(int k=0; k<(*ptr_tag_set_matrix_)[seq_no].size(); ++k){
        Path *ppath = new Path(node_matrix_[seq_no][seq.size()-1][k],&((*ptr_stop_node_)[seq_no]));
        SetPathFeature(std::make_pair(node_matrix_[seq_no][seq.size()-1][k]->GetY(),STOP_NODE_FLAG-seq_no),ppath);
        (*ptr_stop_node_)[seq_no].AddPath(ppath, true);
    }
}

//build right path
void LinearCRF::BuildRPath(std::vector<std::string> seq, int seq_no) {
    for (int i = seq.size()-1; i >=0 ; --i) {
        for(int j=0; j<(*ptr_tag_set_matrix_)[seq_no].size(); ++j){
            if(i == seq.size()-1){
                Path *ppath = new Path(node_matrix_[seq_no][i][j],&((*ptr_stop_node_)[seq_no]));
                SetPathFeature(std::make_pair(node_matrix_[seq_no][i][j]->GetY(),STOP_NODE_FLAG-seq_no),ppath);
                node_matrix_[seq_no][i][j]->AddPath(ppath, false);
            }else{
                for(int k=0; k<(*ptr_tag_set_matrix_)[seq_no].size(); ++k){
                    Path *ppath = new Path(node_matrix_[seq_no][i][j], node_matrix_[seq_no][i+1][k]);
                    SetPathFeature(std::make_pair(node_matrix_[seq_no][i][j]->GetY(),node_matrix_[seq_no][i+1][k]->GetY()),ppath);
                    node_matrix_[seq_no][i][j]->AddPath(ppath,false);
                }
            }
        }
    }
    //create the right path vector for the start node;
    for(int k=0; k<(*ptr_tag_set_matrix_)[seq_no].size(); ++k){
        Path *ppath = new Path(&((*ptr_start_node_)[seq_no]),node_matrix_[seq_no][0][k]);
        SetPathFeature(std::make_pair(START_NODE_FLAG-seq_no,node_matrix_[seq_no][0][k]->GetY()),ppath);
        (*ptr_start_node_)[seq_no].AddPath(ppath, false);
    }
}

void LinearCRF::CalcAllEmpiricalFi() {
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> x_seq = (*ptr_seq_matrix_)[seq_no];
        std::vector<std::string> tag_seq = (*ptr_tag_seq_)[seq_no];
        CalcEmpiricalFi(x_seq,tag_seq,seq_no);
    }
}

//calc E'(f_i)
double LinearCRF::CalcEmpiricalFi(std::vector<std::string> x_seq, std::vector<std::string> tag_seq, int seq_no) {
    //count the feature number
    int tag_size = tag_seq.size();
    //calc count of trans feature
    //start node;
    int index = ptr_feature_->GetFeatureIndex(std::make_pair((START_NODE_FLAG-seq_no), ptr_tag_map_->find(tag_seq[0])->second));
    (*ptr_empirical_e_)[index] = 1;
    (*ptr_feature_bit_vector_)[index] = 1;
    for(int i = 0; i < tag_size; ++i){
        int y_i = 0;
        int y_i_1 = ptr_tag_map_->find(tag_seq[i])->second;
        if(i == tag_size-1){
             y_i = STOP_NODE_FLAG-seq_no;
        } else{
             y_i = ptr_tag_map_->find(tag_seq[i+1])->second;
        }
        index = ptr_feature_->GetFeatureIndex(std::make_pair(y_i_1,y_i));
        //count the num for each tag sequence
        (*ptr_empirical_e_)[index] += 1;
        (*ptr_feature_bit_vector_)[index] += 1;
#ifdef DEBUG_MODE
        std::cout << "the value of tag "<<y_i_1 << ", "<<y_i << " is: "<< (*ptr_empirical_e_)[index] <<",index is:"<<index<<std::endl;
#endif
    }
    int x_size = x_seq.size();
    for(int i=0; i<x_size; ++i){
            int x = ptr_x_corpus_map_->find(x_seq[i])->second + FEATURE_CODE_OFFSET;
            int y = ptr_tag_map_->find(tag_seq[i])->second;
            if(ptr_feature_->GetFeatureMap()->find(std::make_pair(x,y)) != ptr_feature_->GetFeatureMap()->end()){
                 index = ptr_feature_->GetFeatureIndex(std::make_pair(x,y));
                //count the num for each
                (*ptr_empirical_e_)[index] += 1;
                (*ptr_feature_bit_vector_)[index] += 1;
#ifdef DEBUG_MODE
                std::cout << "the value of observ and tag "<<x<< ", "<<y << " is: "<< (*ptr_empirical_e_)[index] <<std::endl;
#endif
            }
    }
    //stop node
    index = ptr_feature_->GetFeatureIndex(std::make_pair(ptr_tag_map_->find(tag_seq[x_seq.size()-1])->second,STOP_NODE_FLAG-seq_no));
    (*ptr_empirical_e_)[index] = 1;
    (*ptr_feature_bit_vector_)[index] = 1;
}

void LinearCRF::CalcFeatureExpectation(std::vector<std::string> seq, int seq_no) {
    //cal expectation of the start node's right paths.
    std::vector<Path *> ppath = (*ptr_start_node_)[seq_no].GetRPath();
    double expectation = 0;
    for (std::vector<Path *>::iterator it = ppath.begin(); it != ppath.end(); ++it) {
        (*it)->CalcLogExpectation((*ptr_Z_)[seq_no], ptr_e_, &expectation);
    }
    int x_size = seq.size();
    //calc expectation for each node and its right path.
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<(*ptr_tag_set_matrix_)[seq_no].size(); ++j){
            node_matrix_[seq_no][i][j]->CalcLogExpectation((*ptr_Z_)[seq_no],ptr_e_);
#ifdef DEBUG_MODE
            std::cout << "expectation of the node"<<i<< ", "<<j<< " is: "<< node_matrix_[i][j]->GetExpectation() <<std::endl;
#endif
        }
    }
#ifdef DEBUG_MODE
    for(int i=0; i<ptr_feature_->GetFeatureSize(); ++i){
        std::cout << "feature expectation of the" <<i<<"th feature is"<<(*ptr_e_)[i]<<std::endl;
    }
#endif
}

//cost indicates w_k * \phi_k(s', s, x)
void LinearCRF::CalcCost(std::vector<std::string> seq, int seq_no) {
    (*ptr_start_node_)[seq_no].SetCost(DEFAULT_COST_VALUE);
    (*ptr_stop_node_)[seq_no].SetCost(DEFAULT_COST_VALUE);
    std::vector< Path*> ptr_start_rpath = (*ptr_start_node_)[seq_no].GetRPath();
    for(std::vector< Path*>::iterator it = ptr_start_rpath.begin(); it!= ptr_start_rpath.end(); ++it){
        ptr_feature_->CalcCost((*it));
    }
    std::vector<Path*> ptr_stop_lpath = (*ptr_stop_node_)[seq_no].GetLPath();
    for(std::vector< Path*>::iterator it = ptr_stop_lpath.begin(); it!= ptr_stop_lpath.end(); ++it){
        ptr_feature_->CalcCost((*it));
    }
    for(int i=0; i<seq.size(); ++i){
        for(int j=0; j<(*ptr_tag_set_matrix_)[seq_no].size(); ++j){
            //calc node cost
            ptr_feature_->CalcCost(node_matrix_[seq_no][i][j]);
            //calc left path cost of each node
            std::vector<Path *> lpath = node_matrix_[seq_no][i][j]->GetLPath();
            std::vector<Path *> rpath = node_matrix_[seq_no][i][j]->GetRPath();
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

void LinearCRF::MainThreadCalculation() {
    int num_of_batch = num_of_instance_ / num_of_thread_;
    int num_of_rest = num_of_instance_ % num_of_thread_;
    CRFThread  *ptr_task = new CRFThread(ptr_start_node_,ptr_stop_node_,node_matrix_,ptr_Z_);
    for(int batch_no = 0; batch_no < num_of_batch; ++batch_no){
        for(int thread_no = 0; thread_no < num_of_thread_; ++thread_no) {
            int seq_no = batch_no * num_of_thread_ + thread_no;
            ThreadCalc(ptr_task,seq_no);
        }
        ThreadStart();
        ptr_thread_vector_->clear();
    }
    for(int seq_no = num_of_instance_ - num_of_rest; seq_no<num_of_instance_; ++seq_no){
        ThreadCalc(ptr_task,seq_no);
    }
    ThreadStart();
    ptr_thread_vector_->clear();
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        CalcFeatureExpectation(seq,seq_no);
    }
}

void LinearCRF::ThreadCalc(CRFThread *ptr_task,int seq_no) {
    int x_size = (*ptr_seq_matrix_)[seq_no].size();
    int y_size = (*ptr_tag_set_matrix_)[seq_no].size();
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    ptr_thread_vector_->push_back(std::thread(&CRFThread::CRFThreadRun,ptr_task,x_size,y_size,seq_no,ptr_weight));

}

void LinearCRF::ThreadStart() {
    for(auto it = ptr_thread_vector_->begin(); it!=ptr_thread_vector_->end(); ++it){
        (*it).join();
    }
}

void LinearCRF::MainCalculation() {
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        clock_t start_time, end_time;
        start_time = clock();
        CalcCost(seq, seq_no);
        end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME__
        std::cout << "calc cost time is: "<<end_time-start_time<<std::endl;
#endif
        start_time = clock();
        ForwardBackward(seq,seq_no);
        end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME__
        std::cout << "forward and backward time is: "<<end_time-start_time<<std::endl;
#endif
        start_time = clock();
        CalcFeatureExpectation(seq,seq_no);
        end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME__
        std::cout << "calc expectation time is: "<<end_time-start_time<<std::endl;
#endif
    }
}

// stochastic gradient descent;
void LinearCRF::CalcGradient() {
    clock_t start_time, end_time;
    start_time = clock();
//    MainCalculation();
    MainThreadCalculation();
    end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME
    std::cout << "main calc time is: "<<end_time-start_time<<std::endl;
#endif
    int size = ptr_feature_->GetFeatureSize();
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    for(int k=0; k<size; ++k){
        double empirical_e_k = (*ptr_empirical_e_)[k];
        double e_k = (*ptr_e_)[k];
        double pre_w_k = (*ptr_weight)[k];
        double penalty = 0; //pre_w_k / (L2_FACTOR * L2_FACTOR);
        (*ptr_gradient_)[k] = empirical_e_k - e_k - penalty;
        //std::cout << "the gradient of the "<<k<<"th feature is: "<<(*ptr_gradient_)[k]<<std::endl;
#ifdef DEBUG_MODE_
        std::cout << "the gradient of the "<<k<<"th feature is: "<<(*ptr_gradient_)[k]<<std::endl;
#endif
    }
}

void LinearCRF::UpdateWeight() {
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    int size = ptr_feature_->GetFeatureSize();
    for(int k=0; k<size; ++k){
        double weight = (*ptr_weight)[k] + LEARNING_RATE * (*ptr_gradient_)[k];
        ptr_feature_->SetWeightVector(k,weight);
    }
}

//calc the loglikelihood function.
double LinearCRF::CalcLoglikelihoodFunction() {
    int feature_size = ptr_feature_->GetFeatureSize();
    std::vector<double > *ptr_weight = ptr_feature_->GetWeightVector();
    double sum_numerator = 0; // (*ptr_weight)[2] + (*ptr_weight)[5] + (*ptr_weight)[6] + (*ptr_weight)[9] + (*ptr_weight)[10];
    for(int i=0; i<feature_size; ++i){
        sum_numerator += (*ptr_weight)[i] * (*ptr_feature_bit_vector_)[i];
    }
    //std::cout << (*ptr_weight)[2]<<","<<(*ptr_weight)[5]<<","<<(*ptr_weight)[6]<<","<<(*ptr_weight)[9]<<","<<(*ptr_weight)[10]<<std::endl;
    double sum_denominator = 0;
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        sum_denominator += (*ptr_Z_)[seq_no];
//        sum_denominator += (*ptr_Z_)[seq_no];
    }
    double l2 = 0;
    /*
    for (int i = 0; i < feature_size; ++i) {
        double w = (*(ptr_feature_->GetWeightVector()))[i];
        l2 += (w * w) / (2 * L2_FACTOR * L2_FACTOR);
    }
    l2 = 0;
     */
    return  sum_numerator - sum_denominator - l2;
}
   
//calc alpha and beta and Z for all nodes in the graph.

void LinearCRF::ForwardBackward(std::vector<std::string> seq, int seq_no) {
    int x_size = seq.size();
    (*ptr_start_node_)[seq_no].SetAlpha(0);
    (*ptr_stop_node_)[seq_no].SetBeta(0);
    int start_time, end_time;
    start_time = clock();
    for(int i=0; i<x_size; ++i){
        for(int j=0; j<(*ptr_tag_set_matrix_)[seq_no].size(); ++j){
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
    for (int i = x_size-1; i>=0 ; --i) {
        for(int j=0; j<(*ptr_tag_set_matrix_)[seq_no].size(); ++j){
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
    std::vector<Path *> stop_lpath= (*ptr_stop_node_)[seq_no].GetLPath();
    for (int j = 0; j<(*ptr_tag_set_matrix_)[seq_no].size(); ++j) {
        Node *p_node = stop_lpath[j]->GetLNode();
//        (*ptr_Z_)[seq_no] = p_node->SumExp((*ptr_Z_)[seq_no],stop_lpath[j]->GetCost()*p_node->GetCost(),p_node->GetAlpha(),(j==0));
        (*ptr_Z_)[seq_no] = p_node->LogSumExp((*ptr_Z_)[seq_no],stop_lpath[j]->GetCost()+p_node->GetCost(),p_node->GetAlpha(),(j==0));
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

void LinearCRF::Viterbi(std::vector<std::string> seq, int seq_no) {
    (*ptr_start_node_)[seq_no].SetBestCost(1);
    for(int i=0; i<seq.size(); ++i){
        for(int j=0; j<(*ptr_tag_set_matrix_)[seq_no].size(); ++j){
            SelectBestNode(node_matrix_[seq_no][i][j]);
        }
    }
    //for the last node;
    std::vector<Path *> ppath = (*ptr_stop_node_)[seq_no].GetLPath();
    for(std::vector<Path *>::iterator it = ppath.begin(); it!=ppath.end(); ++it){
        SelectBestNode(&((*ptr_stop_node_)[seq_no]));
    }
    //backtracking
    ViterbiBackTracking(seq, seq_no);
}

void LinearCRF::SelectBestNode(Node *pNode) {
    double best_cost = 0;
    Node *p_best_node;
    std::vector<Path *> ppath = pNode->GetLPath();
    for(std::vector<Path *>::iterator it = ppath.begin(); it!=ppath.end(); ++it){
        double cost =  (*it)->GetLNode()->GetBestCost() * ((*it)->GetCost() *  (*it)->GetLNode()->GetCost());
        if(cost > best_cost){
            best_cost = cost;
            p_best_node = (*it)->GetLNode();
        }
    }
    pNode->SetBestCost(best_cost);
    pNode->SetPreNode(p_best_node);
    //int x = p_best_node->GetX();
    //int y = p_best_node->GetY();
    //std::cout << "x and y are "<<x<<","<<y<<std::endl;
}

void LinearCRF::ViterbiBackTracking(std::vector<std::string> seq, int seq_no) {
    Node *ptr_node = (*ptr_stop_node_)[seq_no].GetPreNode();
    std::string str = ptr_tag_map_reverse_->find(ptr_node->GetY())->second;
    //std::cout << "tag inferred by viterbi is:"<<str<<std::endl;
    (*ptr_decoded_tag_)[seq_no][seq.size() - 1] = str;
    for (int i = seq.size() - 2; i >= 0; --i) {
        ptr_node = ptr_node->GetPreNode();
        str = ptr_tag_map_reverse_->find(ptr_node->GetY())->second;
        (*ptr_decoded_tag_)[seq_no][i] = str;
    }
}

void LinearCRF::WriteDecodingTagtoFile() {
    std::ofstream of("encodingfile.txt");
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        for(int i=0; i<seq.size(); i++){
            of << seq[i] + " " + (*ptr_decoded_tag_)[seq_no][i];
            of << std::endl;
        }
        of << ". . O";
        of << std::endl;
        of << std::endl;
    }
}

void LinearCRF::PrintPath(Node *pNode) {
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

void LinearCRF::GenerateSeqFromVector(std::vector<std::string> *ptr_vector,
                               std::vector<std::vector<std::string>> *ptr_seq_vector) {
    std::vector<std::string> seq;
    for (std::vector<std::string>::iterator it = ptr_vector->begin(); it != ptr_vector->end(); it++) {
        //std::cout<<*it<<std::endl;
        if (*it == SPERATOR_FLAG) {
            ptr_seq_vector->push_back(seq);
            seq.clear();
            continue;
        } else {
            seq.push_back(*it);
            //do not forget the last seq which doesn't contain a SPEARATOR_FLAG at the end.
            if (it == (ptr_vector->end() - 1)) {
                ptr_seq_vector->push_back(seq);
            }
        }
    }
}

void LinearCRF::CRFRun() {
    //training
    int k = 0;
    while(!is_converged_){
        clock_t start_time, end_time;
        start_time = clock();
        ResetParameters();
        end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME
        std::cout << "reset time is: "<<end_time-start_time<<std::endl;
#endif
        start_time = clock();
        CalcGradient();
        end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME
        std::cout << "calc gradient time is: "<<end_time-start_time<<std::endl;
#endif
        start_time = clock();
        UpdateWeight();
        end_time = clock();
#ifdef PERFORMANCE_CLOCK_TIME
        std::cout << "update weight time is: "<<end_time-start_time<<std::endl;
#endif
        //calc loglikelihood function
        double loss_value = CalcLoglikelihoodFunction();
        std::cout << "loglikelihood is:"<<loss_value<<std::endl;
        if(loss_value > CONVERGED_VALUE && k>100){
            std::cout << "Training completed"<<std::endl;
            is_converged_ = true;
        }
        k++;
    }
}

void LinearCRF::SaveModelToFile() {
    std::ofstream of("modelfile.txt");
    int size = ptr_feature_->GetFeatureSize();
    for (int i = 0; i < size; ++i) {
        double weight = (*ptr_e_)[i];
        of << std::to_string(weight);
        of << std::endl;
    }
}

void LinearCRF::Training() {
    CalcAllEmpiricalFi();
    BuildLattice();
    //for test only;
#ifdef DEBUG_MODE
    PrintPath(ptr_stop_node_);
#endif
    CRFRun();
    SaveModelToFile();
    //run viterbi
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        Viterbi(seq, seq_no);
    }
    WriteDecodingTagtoFile();
}

void LinearCRF::Decoding() {
    BuildLattice();
    std::vector<double> *ptr_weight = ptr_feature_->GetWeightVector();
    std::ifstream is("modelfile.txt");
    std::string weightstr;
    std::string::size_type sz;
    int weight_index = 0;
    while(getline(is,weightstr)){
        double weight = std::stod(weightstr,&sz);
        ptr_feature_->SetWeightVector(weight_index, weight);
        weight_index++;
    }
    //calc cost
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        CalcCost(seq, seq_no);
    }
    //run viterbi
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        Viterbi(seq, seq_no);
    }
    WriteDecodingTagtoFile();
}

