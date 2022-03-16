#ifndef RECSPH_LAST_LALL_H
#define RECSPH_LAST_LALL_H

#include "pRectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class RecSph_Last_lAll {
private:
  unsigned int dim;
  unsigned int tau;
  pRectangle* rect;
  double** cumSumData;
  double** cumSumData2;
  double* vectOfCosts;
  std::list<pSphere> lDiskPass;


public:
  RecSph_Last_lAll(): dim(0), tau(0), rect(0), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL) { }
  RecSph_Last_lAll(unsigned  int p): dim(p), tau(0), rect(new pRectangle(p)),  cumSumData(NULL),cumSumData2(NULL), vectOfCosts(NULL) { }
  RecSph_Last_lAll(unsigned int p, unsigned int t): dim(p), tau(t), rect(new pRectangle(p)), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL) { }
  RecSph_Last_lAll(const RecSph_Last_lAll & candidate);
  ~RecSph_Last_lAll();

  std::list<pSphere> getlDiskPass() const;
  unsigned int getTau()const;
  
  double calculRadius2(Cost & cost, unsigned int i, unsigned int j);
  void cleanOfCandidate();
  bool isEmptyOfCandidate();
  void initialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & vDiskIndexPass);
  void updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<RecSph_Last_lAll>::iterator> &vectlinktocands, unsigned int& realNbExclus);
};
#endif //RECSPH_LAST_LALL_H
