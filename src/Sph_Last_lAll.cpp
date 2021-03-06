#include "Sph_Last_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Sph_Last_lAll::Sph_Last_lAll(const Sph_Last_lAll & candidate) {
  dim = candidate.dim;
  tau = candidate.tau;
  cumSumData = candidate.cumSumData;
  cumSumData2 = candidate.cumSumData2;
  vectOfCosts = candidate.vectOfCosts;
  flEmpty = candidate.flEmpty;
  lDiskPass.clear();
  lDiskPass = candidate.lDiskPass;
}

Sph_Last_lAll::~Sph_Last_lAll() { cumSumData = NULL; cumSumData2 = NULL;  vectOfCosts = NULL; }

unsigned int Sph_Last_lAll::getTau()const { return tau; }

void Sph_Last_lAll::cleanOfCandidate() { cumSumData = NULL; cumSumData2 = NULL; vectOfCosts = NULL; }

bool Sph_Last_lAll::isEmptyOfCandidate() { return flEmpty; }

double Sph_Last_lAll::calculRadius2(Cost & cost, unsigned int i, unsigned int j){
  cost.InitialCost(dim, i, j, cumSumData, cumSumData2, vectOfCosts);
  double radius2 = (vectOfCosts[j + 1] - vectOfCosts[i] - cost.get_coef_Var())/cost.get_coef();
  return radius2;
}


void Sph_Last_lAll::initialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & vDiskIndexPass) {
  tau = t;
  cumSumData = cumsumdata;
  cumSumData2 = cumsumdata2;
  vectOfCosts = vectofcosts;
  flEmpty = false;
  lDiskPass.clear();

  if (vDiskIndexPass.size() != 0) {
    Cost cost = Cost(dim);
    pSphere diskPass = pSphere(dim);
    double radius2;
    for (std::vector<unsigned int>::iterator it = vDiskIndexPass.begin(); it != vDiskIndexPass.end(); ++it) {
      radius2 = calculRadius2(cost, (*it), tau-1);
      diskPass.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
      lDiskPass.push_back(diskPass);
    }
  }
}

void Sph_Last_lAll::updateOfCandidate(unsigned int indexToLinkOfUpdCand, std::vector<std::list<Sph_Last_lAll>::iterator> &vectlinktocands, unsigned int& realNbExclus) {
  realNbExclus = 0;
  Cost cost = Cost(dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> getTau();
  double radius2 = calculRadius2(cost, tau, j);
  if (radius2 < 0) { flEmpty = true; return;}   //pelt

  pSphere diskI = pSphere(dim);
  diskI.InitialpSphere(dim, cost.get_mu(), sqrt(radius2));
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
