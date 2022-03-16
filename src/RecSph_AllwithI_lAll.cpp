#include "RecSph_AllwithI_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

RecSph_AllwithI_lAll::RecSph_AllwithI_lAll(const RecSph_AllwithI_lAll & candidate) {
  dim = candidate.dim;
  tau = candidate.tau;
  rect = new pRectangle(dim);
  cumSumData = candidate.cumSumData;
  cumSumData2 = candidate.cumSumData2;
  vectOfCosts = candidate.vectOfCosts;
  lDiskPass.clear();
  lDiskPass = candidate.lDiskPass;
}

RecSph_AllwithI_lAll::~RecSph_AllwithI_lAll() { delete rect;  cumSumData = NULL; cumSumData2 = NULL; vectOfCosts = NULL; }

std::list<pSphere> RecSph_AllwithI_lAll::getlDiskPass()const { return lDiskPass; }
unsigned int RecSph_AllwithI_lAll::getTau()const { return tau; }

void RecSph_AllwithI_lAll::cleanOfCandidate() { cumSumData = NULL; cumSumData2 = NULL; vectOfCosts = NULL;  lDiskPass.clear(); }
bool RecSph_AllwithI_lAll::isEmptyOfCandidate() { return rect -> IsEmpty_rect(); }

double RecSph_AllwithI_lAll::calculRadius2(Cost & cost, unsigned int i, unsigned int j){
  cost.InitialCost(dim, i, j, cumSumData, cumSumData2, vectOfCosts);
  double radius2 = (vectOfCosts[j + 1] - vectOfCosts[i] - cost.get_coef_Var())/cost.get_coef();
  return radius2;
}

void RecSph_AllwithI_lAll::initialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & vDiskIndexPass) {
  tau = t;
  cumSumData = cumsumdata;
  cumSumData2 = cumsumdata2;
  vectOfCosts = vectofcosts;
  lDiskPass.clear();
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

void RecSph_AllwithI_lAll::updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<RecSph_AllwithI_lAll>::iterator> &vectlinktocands, unsigned int &realNbExclus) {
  realNbExclus = 0;
  Cost cost= Cost(dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> getTau();
  double radius2 = calculRadius2(cost, tau, j);
  if (radius2 < 0) { rect -> DoEmpty_rect(); return;}   //pelt
  pSphere diskI = pSphere(dim);
  diskI.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
  std::list<pSphere>::iterator iter;
  
  //sphere intersection:
  if (indexToLinkOfUpdCand  != (vectlinktocands.size()-1)) {
    pSphere diskIold = pSphere(dim);
    for (unsigned int i = indexToLinkOfUpdCand; i < (vectlinktocands.size()-1); i++) {
      j = vectlinktocands[i] -> getTau();
      radius2 = calculRadius2(cost, tau, j);
      if (radius2 < 0) { rect -> DoEmpty_rect();; return;}   //pelt
      diskIold.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
      if (!diskI.isIntersection(diskIold)) {
        rect -> DoEmpty_rect();
        return;
      }
    }
  }
  //sphere inclusion :
  if (lDiskPass.size() != 0) {
    iter = lDiskPass.begin();
    while( iter != lDiskPass.end()) {
      if ( diskI.isIntersection((*iter))) {
        realNbExclus++;
        if (diskI.isInclusion((*iter))) {
          rect -> DoEmpty_rect();
          return;
        } else { 
          ++iter; 
        } 
      } else {
        iter = lDiskPass.erase(iter);
      }
    }
  }
  //rect intersection :
  for (unsigned int i = indexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    j = vectlinktocands[i] -> getTau();
    radius2 = calculRadius2(cost, tau, j);
    if (radius2 < 0) { rect -> DoEmpty_rect(); return; }//pelt
    diskI.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
    rect -> Intersection_disk(diskI);
    if (rect -> IsEmpty_rect()) { return;}
  }
  //rect exclusion :
  if (lDiskPass.size() != 0) {
    iter = lDiskPass.begin();
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
