#ifndef SPH_ALLWITHI_LALL_H
#define SPH_ALLWITHI_LALL_H
#include "pSphere.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Sph_AllwithI_lAll {
private:
  unsigned int dim;
  unsigned int tau;
  double** cumSumData;
  double** cumSumData2;
  double* vectOfCosts;
  bool flEmpty;
  std::list<pSphere> lDiskPass;
  
public:
  Sph_AllwithI_lAll(): dim(0), tau(0), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL),  flEmpty(false) { }
  Sph_AllwithI_lAll(unsigned  int p): dim(p), tau(0), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL), flEmpty(false) { }
  Sph_AllwithI_lAll(unsigned int p, unsigned int t): dim(p), tau(t), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL), flEmpty(false) { }
  Sph_AllwithI_lAll(const Sph_AllwithI_lAll & candidate);
  ~Sph_AllwithI_lAll();

  unsigned int getTau()const;
  
  double calculRadius2(Cost & cost, unsigned int i, unsigned int j);
  void cleanOfCandidate();
  bool isEmptyOfCandidate();
  void initialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & vDiskIndexPass);
  void updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<Sph_AllwithI_lAll>::iterator> &vectlinktocands, unsigned int& realNbExclus);
};
#endif //SPH_ALLWITHI_LALL_H
//------------------------------------------------------------------------------
