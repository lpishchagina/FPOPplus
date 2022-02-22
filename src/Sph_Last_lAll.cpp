#include "Sph_Last_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Sph_Last_lAll::Sph_Last_lAll(const Sph_Last_lAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  fl_empty = candidate.fl_empty;
  DiskListBefore.clear();
  DiskListBefore = candidate.DiskListBefore;
}

Sph_Last_lAll::~Sph_Last_lAll() { CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Sph_Last_lAll::GetTau()const { return Tau; }

double Sph_Last_lAll::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) { dist = dist + (a[k] - b[k])*(a[k] - b[k]); }
  return sqrt(dist);
}

void Sph_Last_lAll::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Sph_Last_lAll::EmptyOfCandidate() { return fl_empty; }

void Sph_Last_lAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts, std::vector <unsigned int> & DiskIndexBefore) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
  DiskListBefore.clear();

  if (DiskIndexBefore.size() != 0) {
    Cost cost = Cost(Dim);
    pSphere DiskBefore = pSphere(Dim);
    double Radius2;
    for (std::vector<unsigned int>::iterator it = DiskIndexBefore.begin() ; it != DiskIndexBefore.end(); ++it) {
      cost.InitialCost(Dim, (*it), Tau-1, CumSumData, CumSumData2, VectOfCosts);
      Radius2 = (VectOfCosts[Tau] - VectOfCosts[(*it)] - cost.get_coef_Var()) / cost.get_coef();
      DiskBefore.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
      DiskListBefore.push_back(DiskBefore);
    }
  }
}

void Sph_Last_lAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Last_lAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  Cost CostNew = Cost(Dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  CostNew.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - CostNew.get_coef_Var())/CostNew.get_coef();

  if (Radius2 < 0) { fl_empty = true; return;}   //pelt

  pSphere DiskNew = pSphere(Dim);
  DiskNew.InitialpSphere(Dim, CostNew.get_mu(), sqrt(Radius2));
  double dist;
  //inclusion :
  if (DiskListBefore.size() != 0) {//pelt+
    std::list<pSphere>::iterator iter = DiskListBefore.begin();
    while( iter != DiskListBefore.end()) {
      dist = Dist(DiskNew.get_center(),(*iter).get_center());
      if (dist < ((*iter).get_radius() + DiskNew.get_radius())){
        RealNbExclus++;
        if (dist <= ((*iter).get_radius() - DiskNew.get_radius())){
          fl_empty = true;
          return;
        } else {
          ++iter;
        }
      } else {
        iter = DiskListBefore.erase(iter);
      }
    }
  }
}