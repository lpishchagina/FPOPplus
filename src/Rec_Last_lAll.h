#ifndef REC_LAST_LALL_H
#define REC_LAST_LALL_H

#include "pRectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Rec_Last_lAll {
private:
  unsigned int dim;
  unsigned int tau;
  pRectangle* rect;
  double** cumSumData;
  double** cumSumData2;
  double* vectOfCosts;
  std::list<pSphere> lDiskPass;

public:
  Rec_Last_lAll(): dim(0), tau(0), rect(0), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL) { }
  Rec_Last_lAll(unsigned  int p): dim(p), tau(0), rect(new pRectangle(p)), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL) { }
  Rec_Last_lAll(unsigned int p, unsigned int t): dim(p), tau(t), rect(new pRectangle(p)), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL) { }
  Rec_Last_lAll(const Rec_Last_lAll & candidate);
  ~Rec_Last_lAll();
  
  unsigned int getTau()const;
  std::list<pSphere> getlDiskPass() const;
  
  double calculRadius2(Cost & cost, unsigned int i, unsigned int j);
  void cleanOfCandidate();
  bool isEmptyOfCandidate();
  void initialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts,  std::vector <unsigned int> & lDiskIndexPass);
  void updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<Rec_Last_lAll>::iterator> &vectlinktocands, unsigned int& realNbExclus);
};
#endif //REC_LAST_LALL_H
