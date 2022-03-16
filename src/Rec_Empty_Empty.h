#ifndef Rec_Empty_Empty_H
#define Rec_Empty_Empty_H

#include "pRectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Rec_Empty_Empty{
private:
  unsigned int dim;
  unsigned int tau;
  double** cumSumData;
  double** cumSumData2;
  double* vectOfCosts;
  bool flEmpty;

public:
  Rec_Empty_Empty(): dim(0), tau(0), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL), flEmpty(false) { }
  Rec_Empty_Empty(unsigned  int p): dim(p), tau(0),  cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL),  flEmpty(false) { }
  Rec_Empty_Empty(unsigned int p, unsigned int t): dim(p), tau(t), cumSumData(NULL), cumSumData2(NULL), vectOfCosts(NULL),  flEmpty(false) { }
  Rec_Empty_Empty(const Rec_Empty_Empty & candidate);
  ~Rec_Empty_Empty();

  unsigned int getTau()const;

  double calculRadius2(Cost & cost, unsigned int i, unsigned int j);
  void cleanOfCandidate();
  bool isEmptyOfCandidate();
  void initialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & vDiskIndexPass);
  void updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<Rec_Empty_Empty>::iterator> &vectlinktocands, unsigned int& realNbExclus);
};
#endif //Rec_Empty_Empty_H
