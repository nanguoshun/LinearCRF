//
// Created by  ngs on 12/06/2018.
//

#ifndef CRF_NODE_H
#define CRF_NODE_H

#include <vector>
#include <string>

class Path;

class Node{
public:
    explicit Node(int x, int y);
    ~Node();
    void CalcAlpha();
    void CalcBeta();
    void CalcExpectation(double Z);
    void SetCost(double cost);
    /**
     *
     * @param iteration_alpha : for each \alpha_{i-1} * cost
     * @param cost : w_k * f_k(s', s, x)
     * @param lalpha: \alpha_{i-1}
     * @param isStart: check the based case
     * @return
     */
    inline double SumExp(double iteration_alpha, double cost, double lalpha, bool isStart){
        // calc \alpha * cost;
        double value = cost * lalpha;
        if(isStart){
            return value; // init mode
        }else{
            return  (iteration_alpha + value);
        }
    }
    double GetAlpha();
    double GetBeta();
    int GetX();
    int GetY();
    double GetCost();
    std::vector<Path *> GetLPath();
    std::vector<Path *> GetRPath();
    double GetExpectation();
    void AddPath(Path *ptr_path, bool isLpath);
private:
    int x_;
    int y_;
    double cost_;
    double alpha_;
    double beta_;
    double expectation_;
    double bestCost_; //for viterbi;
    Node *pre_; // for viterbi;
    std::vector<Path *> lpath_; //edges in the left
    std::vector<Path *> rpath_; //edges in the right
};

#endif //CRF_NODE_H
