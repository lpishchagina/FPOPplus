#include "Sph_AllwithI_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Sph_AllwithI_lAll::Sph_AllwithI_lAll(const Sph_AllwithI_lAll & candidate) {
  dim = candidate.dim;
  tau = candidate.tau;
  cumSumData = candidate.cumSumData;
  cumSumData2 = candidate.cumSumData2;
  vectOfCosts = candidate.vectOfCosts;
  flEmpty = candidate.flEmpty;
  lDiskPass.clear();
  lDiskPass = candidate.lDiskPass;
}

Sph_AllwithI_lAll::~Sph_AllwithI_lAll() { cumSumData = NULL; cumSumData2 = NULL;  vectOfCosts = NULL; }

unsigned int Sph_AllwithI_lAll::getTau()const { return tau; }

double Sph_AllwithI_lAll::calculRadius2(Cost & cost, unsigned int i, unsigned int j){
  cost.InitialCost(dim, i, j, cumSumData, cumSumData2, vectOfCosts);
  double radius2 = (vectOfCosts[j + 1] - vectOfCosts[i] - cost.get_coef_Var())/cost.get_coef();
  return radius2;
}

void Sph_AllwithI_lAll::cleanOfCandidate() { cumSumData = NULL; cumSumData2 = NULL; vectOfCosts = NULL; }

bool Sph_AllwithI_lAll::isEmptyOfCandidate() { return flEmpty; }

void Sph_AllwithI_lAll::initialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & vDiskIndexPass) {
  tau = t;
  cumSumData = cumsumdata;
  cumSumData2 = cumsumdata2;
  vectOfCosts = vectofcosts;
  flEmpty = false;
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

void Sph_AllwithI_lAll::updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<Sph_AllwithI_lAll>::iterator> &vectlinktocands, unsigned int& realNbExclus) {
  realNbExclus = 0;
  Cost cost = Cost(dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> getTau();
  double radius2 = calculRadius2(cost, tau, j);
  if (radius2 < 0) { flEmpty = true; return;}   //pelt
  pSphere diskI = pSphere(dim);
  diskI.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
  //sphere intersection :
  if (indexToLinkOfUpdCand  != (vectlinktocands.size()-1)) {
    pSphere diskIold = pSphere(dim);
    for (unsigned int i = indexToLinkOfUpdCand; i < (vectlinktocands.size()-1); i++) {
      j = vectlinktocands[i] -> getTau();
      radius2 = calculRadius2(cost, tau, j);
      if (radius2 < 0) { flEmpty = true; return;}   //pelt
      diskIold.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
      if (!diskI.isIntersection(diskIold)) {
        flEmpty = true; 
        return;
      }
    }
  }
  //sphere inclusion :
  if (lDiskPass.size() != 0) {
    std::list<pSphere>::iterator iter = lDiskPass.begin();
    while( iter != lDiskPass.end()) {
      if ( diskI.isIntersection((*iter))) {
        realNbExclus++;
        if (diskI.isInclusion((*iter))) {
          flEmpty = true; 
          return;
        } else { 
          ++iter; 
        } 
      } else {
        iter = lDiskPass.erase(iter);
      }
    }
  }
}
