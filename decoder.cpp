//
// Created by  ngs on 24/06/2018.
//
#include <sstream>
#include "decoder.h"
#include <thread>

Decoder::Decoder(DatasetMgr *ptr_datamgr) {
    ptr_datamgr_ = ptr_datamgr;
    ptr_seq_matrix_ = new std::vector<std::vector<std::string>>;
    ptr_tag_seq_ = new std::vector<std::vector<std::string>>;
    ptr_feature_ = new Feature("featuremap.txt");
    ptr_x_corpus_map_ = new std::map<std::string, int>;
    ptr_x_corpus_map_reverse_ = new std::map<int, std::string>;
    ptr_tag_map_ = new std::map<std::string, int>;
    ptr_tag_map_reverse_ = new std::map<int, std::string>;
    num_of_instance_ = ptr_datamgr_->GetNumOfTrainingSeqs();
    ptr_tag_set_ = ptr_datamgr_->GetTagSet();
    tag_num_ = ptr_tag_set_->size();
    ptr_tag_vector_ = ptr_datamgr_->GetTageVector();
    ptr_x_test_vector_ = ptr_datamgr_->GetTrainingXVector();
    ptr_x_test_set_ = ptr_datamgr_->GetTrainingXSet();
    GenerateSeqFromVector(ptr_x_test_vector_,ptr_seq_matrix_);
    GenerateSeqFromVector(ptr_tag_vector_,ptr_tag_seq_);
    CreateTagObservMap("x_map.txt","tag_map.txt");
    ptr_start_node_ = new std::vector<Node>;
    ptr_stop_node_ = new std::vector<Node>;
    ptr_decoded_tag_ = new std::vector<std::vector<std::string>>();
    ptr_tag_set_matrix_ = new std::vector<std::set<std::string>>;
    for(int seq_no = 0; seq_no<num_of_instance_; ++seq_no){
        std::cout << (*ptr_seq_matrix_)[seq_no].size() <<std::endl;
        std::vector<std::string> *pvector = new std::vector<std::string>((*ptr_seq_matrix_)[seq_no].size());
        ptr_decoded_tag_->push_back(*pvector);
        std::set<std::string> *pset = new std::set<std::string>;
        ptr_tag_set_matrix_->push_back(*pset);
    }
    FromVectorToSet();
    ptr_thread_vector_ = new std::vector<std::thread>;
    num_of_thread_ =  std::thread::hardware_concurrency();

}

Decoder::Decoder() {

}

Decoder::~Decoder() {
    delete ptr_seq_matrix_;
    delete ptr_tag_seq_;
    delete ptr_feature_;
    delete ptr_x_corpus_map_;
    delete ptr_x_corpus_map_reverse_;
    delete ptr_tag_map_;
    delete ptr_tag_map_reverse_;
    delete ptr_start_node_;
    delete ptr_stop_node_;
    DeleteLattice();
}

void Decoder::DeleteLattice() {
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


void Decoder::ReadModel() {

    std::ifstream is("file.txt");
    std::string weightstr;
    std::string::size_type sz;
    int weight_index = 0;
    while(getline(is,weightstr)){
        double weight = std::stod(weightstr,&sz);
        ptr_feature_->SetWeightVector(weight_index, weight);
        weight_index++;
    }
}

void Decoder::CreateTagObservMap(const char *x_map_file, const char* tag_map_file){
    int index = ReadMapfile(x_map_file,ptr_x_corpus_map_,ptr_x_corpus_map_reverse_);
    for (std::set<std::string>::iterator it = ptr_x_test_set_->begin(); it != ptr_x_test_set_->end(); ++it) {
        InsertMap(index, (*it), ptr_x_corpus_map_, ptr_x_corpus_map_reverse_, &index);
    }
    index = ReadMapfile(tag_map_file,ptr_tag_map_, ptr_tag_map_reverse_);
    for(std::set<std::string>::iterator it = ptr_tag_set_->begin(); it!= ptr_tag_set_->end(); ++it){
        InsertMap(index, (*it), ptr_tag_map_, ptr_tag_map_reverse_, &index);
    }
}

int Decoder::ReadMapfile(const char *ptr_x_map_file, std::map<std::string, int> *ptr_map,
                          std::map<int, std::string> *ptr_map_reverse_) {
    std::ifstream ifs(ptr_x_map_file);
    std::string str;
    int map_size = 0;
    while(getline(ifs,str)){
        std::stringstream ss(str);
        int index = 0;
        std::string map_str;
        for(int i=0; i<2; ++i){
            if(i==0){
                ss >> index;
            } else{
                ss >> map_str;
            }
        }
        InsertMap(index, map_str, ptr_map, ptr_map_reverse_, &map_size);
    }
    return map_size;
}

void Decoder::InsertMap(int index, std::string str, std::map<std::string, int> *pmap,
                            std::map<int, std::string> *pmap_reverse, int *map_size) {
    if(pmap->find(str) == pmap->end()){
        pmap->insert(std::make_pair(str, index));
        pmap_reverse->insert(std::make_pair(index,str));
        (*map_size)++;
    }
}

void Decoder::BuildLattice() {
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
        std::cout << seq.size() <<std::endl;
        BuildRPath(seq,seq_no);
    }
}



void Decoder::BuildNode(std::vector<std::string> seq, int seq_no) {
    //create each row
    for (int i = 0; i < seq.size(); ++i) {
        std::vector<Node *> *ptr_node_vector = new std::vector<Node *>;
        node_matrix_[seq_no].push_back(*ptr_node_vector);
    }
    //create node in each row
    for (int i = 0; i < seq.size(); ++i) {
        for (std::set<std::string>::iterator it = (*ptr_tag_set_matrix_)[seq_no].begin(); it != (*ptr_tag_set_matrix_)[seq_no].end(); ++it) {
            //std::cout << "the string value of x and y are: " << seq[i] << ", " << (*it) << std::endl;
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


void Decoder::SetPathFeature(std::pair<int, int> feature_pair, Path *ppath) {
    int feature_index = FEATURE_NO_EXIST;
    if(ptr_feature_->GetFeatureMap()->find(feature_pair) != ptr_feature_->GetFeatureMap()->end()){
        feature_index = ptr_feature_->GetFeatureMap()->find(feature_pair)->second;
    }
    ppath->SetFeatureIndex(feature_index);
}

//build left path
void Decoder::BuildLPath(std::vector<std::string> seq, int seq_no) {
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
void Decoder::BuildRPath(std::vector<std::string> seq, int seq_no) {
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
    std::cout << seq.size() <<std::endl;

}

void Decoder::GenerateSeqFromVector(std::vector<std::string> *ptr_vector,
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
            std::cout << "x is: "<<(*it)<<std::endl;
            //do not forget the last seq which doesn't contain a SPEARATOR_FLAG at the end.
            if (it == (ptr_vector->end() - 1)) {
                ptr_seq_vector->push_back(seq);
            }
        }
    }
}

void Decoder::CalcCost() {
    int num_of_batch = num_of_instance_ / num_of_thread_;
    int num_of_rest = num_of_instance_ % num_of_thread_;
    CRFThread  *ptr_task = new CRFThread(ptr_start_node_,ptr_stop_node_,node_matrix_);
    for(int batch_no = 0; batch_no < num_of_batch; ++batch_no){
        for(int thread_no = 0; thread_no < num_of_thread_; ++thread_no) {
            int seq_no = batch_no * num_of_thread_ + thread_no;
            ThreadCalc(ptr_task,seq_no);
        }
        ThreadStart();
        ptr_thread_vector_->clear();
    }
    //calc the rest of the
    for(int seq_no = num_of_instance_ - num_of_rest; seq_no<num_of_instance_; ++seq_no){
        ThreadCalc(ptr_task,seq_no);
    }
    ThreadStart();
    ptr_thread_vector_->clear();
}

void Decoder::ThreadCalc(CRFThread *ptr_task,int seq_no) {
    int x_size = (*ptr_seq_matrix_)[seq_no].size();
    int y_size = (*ptr_tag_set_matrix_)[seq_no].size();
    ptr_thread_vector_->push_back(std::thread(&CRFThread::CRFThreadViterbi,ptr_task,x_size,y_size,seq_no,ptr_feature_->GetWeightVector()));

}

void Decoder::ThreadStart() {
    for(auto it = ptr_thread_vector_->begin(); it!=ptr_thread_vector_->end(); ++it){
        (*it).join();
    }
}


void Decoder::CalcCost(std::vector<std::string> seq, int seq_no) {
    (*ptr_start_node_)[seq_no].SetCost(DEFAULT_COST_VALUE_DEC);
    (*ptr_stop_node_)[seq_no].SetCost(DEFAULT_COST_VALUE_DEC);
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

void Decoder::Viterbi(std::vector<std::string> seq, int seq_no) {
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

void Decoder::SelectBestNode(Node *pNode) {
    double best_cost = 0;
    Node *p_best_node;
    std::vector<Path *> ppath = pNode->GetLPath();
    for(std::vector<Path *>::iterator it = ppath.begin(); it!=ppath.end(); ++it){
        double cost =  (*it)->GetLNode()->GetBestCost() * ((*it)->GetCost() * (*it)->GetLNode()->GetCost());
        if(cost >= best_cost){
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

void Decoder::ViterbiBackTracking(std::vector<std::string> seq, int seq_no) {
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

void Decoder::WriteDecodingTagtoFile() {
    std::ofstream of("encodingfile.txt");
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        for(int i=0; i<seq.size(); ++i){
            of << seq[i] + " " + (*ptr_tag_seq_)[seq_no][i] + " " + (*ptr_decoded_tag_)[seq_no][i];
            of << std::endl;
        }
        of << ". . O";
        of << std::endl;
        of << std::endl;
    }
}

void Decoder::RewriteTrainandTestData(const char *origfile, const char *newfile) {
//    std::ifstream ifs("train.txt");
//    std::ofstream ofs("newtrain.txt");
    std::ifstream ifs(origfile);
    std::ofstream ofs(newfile);
    std::string str;
    while (getline(ifs,str)){
        std::stringstream ss(str);
        std::string x;
        std::string tag;
        std::string bio;
        ss >> x;
        ss >> tag;
        ss >> bio;
//        ofs << x + " "+bio+" "+tag<<std::endl;
//        ofs << x + " "+tag+" "+tag<<std::endl;
        ofs <<x+" "+tag<<std::endl;
//        ofs << x +" "+tag<<std::endl;
    }
}

void Decoder::CalculateResult() {
//    std::ifstream ifs("encodingfile.txt");
    std::ifstream ifs("test_info_test");
    std::string str;
    double correct_prediction = 0;
    double all_prediction = 0;
    while(getline(ifs,str)){
        std::string predicted_str;
        std::string ground_truth;
        std::stringstream ss(str);
        std::string tmp;
        ss >> tmp;
        ss >> ground_truth;
        ss >> predicted_str;
        if("." != tmp){
            all_prediction ++;
            if(ground_truth == predicted_str){
                correct_prediction++;
            }
        }
    }
    double correct_ratio = correct_prediction / all_prediction;
    std::cout << "correct ratio is "<<correct_ratio<<std::endl;
}

void Decoder::FromVectorToSet() {
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

void Decoder::Decoding() {
    BuildLattice();
    ReadModel();
    //calc cost
    //CalcCost();
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::cout << (*ptr_seq_matrix_)[seq_no].size() <<std::endl;
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        CalcCost(seq, seq_no);
    }
    //run viterbi
    for(int seq_no = 0; seq_no < num_of_instance_; ++seq_no) {
        std::vector<std::string> seq = (*ptr_seq_matrix_)[seq_no];
        Viterbi(seq, seq_no);
        std::cout << "run viterbi for "<< seq_no << "th sentences"<<std::endl;
    }
    WriteDecodingTagtoFile();
}