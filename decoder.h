//
// Created by  ngs on 24/06/2018.
//

#ifndef CRF_DECODE_H
#define CRF_DECODE_H

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include "datasetmgr.h"
#include "common.h"
#include "node.h"
#include "path.h"
#include "featuremanager.h"

class Decoder{
public:
    explicit Decoder(DatasetMgr *ptr_datamgr);
    Decoder();
    ~Decoder();

    void DeleteLattice();
    void CreateTagObservMap(const char *x_map_file, const char* tag_map_file);

    int ReadMapfile(const char *ptr_x_map_file, std::map<std::string, int> *ptr_map,
                     std::map<int, std::string> *ptr_map_reverse_);
    void InsertMap(int index, std::string str, std::map<std::string, int> *pmap,
                   std::map<int, std::string> *pmap_reverse, int *map_size);
    void Decoding();
    void ReadModel();
    void ReadCoupus();
    void CalcCost(std::vector<std::string> seq, int seq_no);
    void Viterbi();
    void SaveResult();

    void BuildLattice();
    void BuildNode(std::vector<std::string> seq, int seq_no);
    void BuildLPath(std::vector<std::string> seq, int seq_no);
    void BuildRPath(std::vector<std::string> seq, int seq_no);

    void SetPathFeature(std::pair<int,int> feature_pair, Path *ppath);
    void GenerateSeqFromVector(std::vector<std::string> *ptr_vector,std::vector<std::vector<std::string>> *ptr_seq_vector);

    void Viterbi(std::vector<std::string> seq, int seq_no);
    void SelectBestNode(Node *pNode);
    void ViterbiBackTracking(std::vector<std::string> seq, int seq_no);
    void WriteDecodingTagtoFile();
    void CalculateResult();
    void RewriteTrainandTestData(const char *origfile, const char *newfile);
    void FromVectorToSet();

private:
    //
    DatasetMgr *ptr_datamgr_;
    //lattice related
    std::vector<std::vector<std::vector<Node *>>> node_matrix_;
    std::vector<Node> *ptr_start_node_;
    std::vector<Node> *ptr_stop_node_;

    Feature* ptr_feature_;

    //data related
    int num_of_instance_;
    int tag_num_;
    std::vector<std::vector<std::string>> *ptr_seq_matrix_;
    std::vector<std::vector<std::string>> *ptr_tag_seq_;
    std::vector<std::string> *ptr_x_test_vector_;
    std::set<std::string> *ptr_x_test_set_;
    std::vector<std::string> *ptr_tag_vector_;
    std::set<std::string> *ptr_tag_set_;
    std::vector<std::vector<std::string>> *ptr_decoded_tag_;
    //corpus
    std::map<std::string, int> *ptr_tag_map_;
    std::map<int, std::string> *ptr_tag_map_reverse_;
    std::map<std::string, int> *ptr_x_corpus_map_;
    std::map<int, std::string> *ptr_x_corpus_map_reverse_;
    std::vector<std::set<std::string>> *ptr_tag_set_matrix_;


};

#endif //CRF_DECODE_H
