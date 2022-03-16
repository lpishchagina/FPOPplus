#include "Rec_Last_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Rec_Last_lAll::Rec_Last_lAll(const Rec_Last_lAll & candidate) {
  dim = candidate.dim;
  tau = candidate.tau;
  rect = new pRectangle(dim);
  cumSumData = candidate.cumSumData;
  cumSumData2 = candidate.cumSumData2;
  vectOfCosts = candidate.vectOfCosts;
  lDiskPass.clear();
  lDiskPass = candidate.lDiskPass;
}

Rec_Last_lAll::~Rec_Last_lAll() { delete rect;  cumSumData = NULL;  cumSumData2 = NULL;  vectOfCosts = NULL; }

unsigned int Rec_Last_lAll::getTau()const { return tau; }
std::list<pSphere> Rec_Last_lAll::getlDiskPass()const { return lDiskPass; }


double Rec_Last_lAll::calculRadius2(Cost & cost, unsigned int i, unsigned int j){
  cost.InitialCost(dim, i, j, cumSumData, cumSumData2, vectOfCosts);
  double radius2 = (vectOfCosts[j + 1] - vectOfCosts[i] - cost.get_coef_Var())/cost.get_coef();
  return radius2;
}

void Rec_Last_lAll::cleanOfCandidate() { cumSumData = NULL;  cumSumData2 = NULL;  vectOfCosts = NULL; lDiskPass.clear();}bool Rec_Last_lAll::isEmptyOfCandidate() { return rect -> IsEmpty_rect(); }

void Rec_Last_lAll::initialOfCandidate(unsigned int t, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts,  std::vector <unsigned int> & vDiskIndexPass) {
  tau = t;
  cumSumData = cumsumdata;
  cumSumData2 = cumsumdata2;
  vectOfCosts = vectofcosts;
  lDiskPass.clear();
  //create lDiskPass
  if (vDiskIndexPass.size() != 0) {
    Cost cost = Cost(dim);
    pSphere diskPass = pSphere(dim);
    double radius2;
    for (std::vector<unsigned int>::iterator it = vDiskIndexPass.begin() ; it != vDiskIndexPass.end(); ++it) {
      radius2 = calculRadius2(cost, (*it), (tau-1));
        
      diskPass.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
      lDiskPass.push_back(diskPass);
    }
  }
}

void Rec_Last_lAll::updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<Rec_Last_lAll>::iterator> &vectlinktocands, unsigned int& realNbExclus) {
  realNbExclus = 0;
  Cost cost = Cost(dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> getTau();
  double radius2 = calculRadius2(cost,tau, j); 
  if (radius2 < 0) { rect -> DoEmpty_rect(); return;}   //pelt
  pSphere diskI = pSphere(dim);
  diskI.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
  //rect last intersection :
  rect -> Intersection_disk(diskI);
  if (rect -> IsEmpty_rect()) { return; }
  //rect exclusion :
  if (lDiskPass.size() != 0) {
    std::list<pSphere>::iterator iter = lDiskPass.begin();
    while ( (iter != lDiskPass.end()) && (!rect -> IsEmpty_rect())) {
      if (rect -> EmptyIntersection(*iter)) {
        iter = lDiskPass.erase(iter);
      } else {
        rect -> Exclusion_disk(*iter);
        ++iter;
        ++realNbExclus;
      }
    }
  }
}
