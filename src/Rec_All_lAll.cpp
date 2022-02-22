#include "Rec_All_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Rec_All_lAll::Rec_All_lAll(const Rec_All_lAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  DiskListBefore.clear();
  DiskListBefore = candidate.DiskListBefore;
}

Rec_All_lAll::~Rec_All_lAll() { delete Rect;  CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

std::list<pSphere> Rec_All_lAll::GetDiskListBefore()const { return DiskListBefore; }
unsigned int Rec_All_lAll::GetTau()const { return Tau; }

void Rec_All_lAll::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL;  DiskListBefore.clear(); }
bool Rec_All_lAll::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Rec_All_lAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts,  std::vector <unsigned int>& DiskIndexBefore) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
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

double Rec_All_lAll::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Rec_All_lAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_lAll>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  RealNbExclus = 0;

  unsigned int j;
  double Radius2;
  Cost CostNew = Cost(Dim);
  pSphere DiskNew = pSphere(Dim);
  //intersection :
  for (unsigned int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    j = vectlinktocands[i] -> GetTau();
    CostNew.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - CostNew.get_coef_Var())/CostNew.get_coef();
    if (Radius2 < 0) { Rect -> DoEmpty_rect(); return; }     //pelt

    DiskNew.InitialpSphere(Dim, CostNew.get_mu(), sqrt(Radius2));
    Rect -> Intersection_disk(DiskNew);//fpop intersection
    if (Rect -> IsEmpty_rect()) { return; }
  }
  //exclusion :
  if (DiskListBefore.size() != 0) {
    std::list<pSphere>::iterator iter = DiskListBefore.begin();
    while ( (iter != DiskListBefore.end()) && (!Rect -> IsEmpty_rect())) {
      if (Rect -> EmptyIntersection(*iter)) {
        iter = DiskListBefore.erase(iter);
      } else {
        Rect -> Exclusion_disk(*iter);//fpop exclusion
        ++iter;
        ++RealNbExclus;
      }
    }
  }
}

