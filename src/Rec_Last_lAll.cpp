#include "Rec_Last_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Rec_Last_lAll::Rec_Last_lAll(const Rec_Last_lAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  DiskListBefore.clear();
  DiskListBefore = candidate.DiskListBefore;
}

Rec_Last_lAll::~Rec_Last_lAll() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Rec_Last_lAll::GetTau()const { return Tau; }
std::list<pSphere> Rec_Last_lAll::GetDiskListBefore()const { return DiskListBefore; }
void Rec_Last_lAll::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; DiskListBefore.clear();}
bool Rec_Last_lAll::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }
void Rec_Last_lAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts,  std::vector <unsigned int> & DiskIndexBefore) {
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

double Rec_Last_lAll::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) { dist = dist + (a[k] - b[k])*(a[k] - b[k]); }
  return sqrt(dist);
}

void Rec_Last_lAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Last_lAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  Cost CostNew = Cost(Dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  CostNew.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - CostNew.get_coef_Var()) / CostNew.get_coef();

  if (Radius2 < 0) { Rect -> DoEmpty_rect(); return;}   //pelt

  pSphere DiskNew = pSphere(Dim);
  DiskNew.InitialpSphere(Dim, CostNew.get_mu(), sqrt(Radius2));
  //intersection :
  Rect -> Intersection_disk(DiskNew);
  if (Rect -> IsEmpty_rect()) { return; }
  //exclusion :
  if (DiskListBefore.size() != 0) {
    std::list<pSphere>::iterator iter = DiskListBefore.begin();
    while ( (iter != DiskListBefore.end()) && (!Rect -> IsEmpty_rect())) {
      if (Rect -> EmptyIntersection(*iter)) {
        iter = DiskListBefore.erase(iter);
      } else {
        Rect -> Exclusion_disk(*iter);
        ++iter;
        ++RealNbExclus;
      }
    }
  }
}
