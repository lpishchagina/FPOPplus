#ifndef REC_ALL_LALL_H
#define REC_ALL_LALL_H

#include "pRectangle.h"
#include "Cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Rec_All_lAll {
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  std::list<pSphere> DiskListBefore;

public:
  Rec_All_lAll(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) {}
  Rec_All_lAll(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)),  CumSumData(NULL),CumSumData2(NULL), VectOfCosts(NULL) {}
  Rec_All_lAll(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) {}
  Rec_All_lAll(const Rec_All_lAll & candidate);
  ~Rec_All_lAll();

  double Dist(double* a, double*b);
  std::list<pSphere> GetDiskListBefore() const;
  unsigned int GetTau()const;
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & DiskIndexBefore);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_lAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //REC_ALL_LALL_H
