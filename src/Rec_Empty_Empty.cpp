#include "Rec_Empty_Empty.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Rec_Empty_Empty::Rec_Empty_Empty(const Rec_Empty_Empty & candidate) {
  dim = candidate.dim;
  tau = candidate.tau;
  cumSumData = candidate.cumSumData;
  cumSumData2 = candidate.cumSumData;
  vectOfCosts = candidate.vectOfCosts;
}

Rec_Empty_Empty::~Rec_Empty_Empty() { cumSumData = NULL; cumSumData2 = NULL;  vectOfCosts = NULL; }

unsigned int Rec_Empty_Empty::getTau()const { return tau; }

double Rec_Empty_Empty::calculRadius2(Cost & cost, unsigned int i, unsigned int j){
  cost.InitialCost(dim, i, j, cumSumData, cumSumData2, vectOfCosts);
  double radius2 = (vectOfCosts[j + 1] - vectOfCosts[i] - cost.get_coef_Var())/cost.get_coef();
  return radius2;
}

void Rec_Empty_Empty::cleanOfCandidate() { cumSumData = NULL;  cumSumData2 = NULL; vectOfCosts = NULL; flEmpty = false; }

bool Rec_Empty_Empty::isEmptyOfCandidate() { return flEmpty; }

void Rec_Empty_Empty::initialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & vDiskIndexPass) {
  tau = t;
  cumSumData = cumsumdata;
  cumSumData2 = cumsumdata2;
  vectOfCosts = vectofcosts;
  flEmpty = false;
}

void Rec_Empty_Empty::updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<Rec_Empty_Empty>::iterator> &vectlinktocands, unsigned int& realNbExclus) {
  flEmpty = false;
  //pelt
  Cost cost = Cost(dim);
  unsigned int lastT = vectlinktocands[vectlinktocands.size()-1] -> getTau();
  if (calculRadius2(cost,  tau, lastT) < 0) {//pelt
    flEmpty = true;
  }
}
