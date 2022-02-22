#include "FPOP.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//converting parameters("approximation", "intersection", "exclusion") to a numeric value.
int NmbOfapproxFPOP(std::string approximation, std::string intersection,  std::string exclusion) {
  int type_approx = 0;

  if ((intersection == "empty") && (exclusion == "empty")) { type_approx = 1; }//PELT

  if (approximation == "rectangle") {
    if  (intersection == "last")  { type_approx = 3; }//LastFPOP
    if  (intersection == "all")  { type_approx = 5; }//FPOP
  }

  if (approximation == "sphere") { type_approx = 2;}

  if (approximation == "sphere_rectangle") {
    if  (intersection == "last")  { type_approx = 4; }
    if  (intersection == "all")  { type_approx = 6; }
  }

  return type_approx;
}

//Comparison of two FPOP-methods
bool TestOfComparisonTwoFPOP(Rcpp::NumericMatrix data, double penalty, unsigned int type_approx2, double UnpenalizedCost1, unsigned int* LastChpt1) {
  bool res = false;
  if (type_approx2 == 1) {
    FPOP<Rec_Empty_Empty> X = FPOP<Rec_Empty_Empty>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else  if (type_approx2 == 2) {
    FPOP<Sph_Last_lAll> X = FPOP<Sph_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 3) {
    FPOP<Rec_Last_lAll> X = FPOP<Rec_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  }  else if (type_approx2 == 4) {
    FPOP<RecSph_Last_lAll> X = FPOP<RecSph_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 5) {
    FPOP<Rec_All_lAll> X = FPOP<Rec_All_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 6) {
    FPOP<RecSph_All_lAll> X = FPOP<RecSph_All_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  }
  return res;
}


//' @title FPOPplus
//'
//' @description FPOP method (using the rectangle approximation of the sets) for  the multiple changepoint detection.
//' @param data is a matrix of data (p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param approximation is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection is the type of intersection : 'empty', 'all', 'last' or 'random' (by default, 'last').
//' @param exclusion is the type of intersection : 'empty', 'all'or 'random'(by default, 'all').
//' @param NbOfCands is the logical parameter (if NbOfCands = TRUE, than the file "NbOfCands.txt" contains the number of change candidates for each iteration.
//'
//' @return a list of  elements  = (changes, means, UnpenalizedCost, NumberOfCandidates).
//'
//' \describe{
//' \item{\code{changes}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
//' \item{\code{means}}{is the list of successive means for the p-variate time series.}
//' \item{\code{UnpenalizedCost}}{is a number equal to the global cost.}
//' \item{\code{NumberOfCandidates}}{is a number of candidates at each iteration (vector).}
//' }
//'
//' @examples
//' N <- 100
//' Dim <- 2
//' Penalty <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
//' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere', intersection = 'last', exclusion = 'all',NbOfCands = TRUE)
//' FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'last', exclusion = 'all',NbOfCands = TRUE)
//' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere_rectangle', intersection = 'last', exclusion = 'all',NbOfCands = TRUE)
//' FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'all', exclusion = 'all',NbOfCands = TRUE)
//' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere_rectangle', intersection = 'all', exclusion = 'all', NbOfCands = TRUE)
//' FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'empty', exclusion = 'empty', NbOfCands = TRUE)

// [[Rcpp::export]]
List FPOPplus(Rcpp::NumericMatrix data, double penalty, std::string approximation = "rectangle", std::string intersection = "all",  std::string exclusion = "all", bool NbOfCands = false) {
  List res;
  int type_approx = NmbOfapproxFPOP(approximation, intersection, exclusion);
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {throw std::range_error("Penalty should be a non-negative number!");}
  if(type_approx == 0){throw std::range_error("This combination of parameters 'intersection' and 'exclusion' is not available. ");}
  //----------------------------------------------------------------------------
  if (type_approx == 1) {
   FPOP<Rec_Empty_Empty> X = FPOP<Rec_Empty_Empty>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else  if (type_approx == 2) {
    FPOP<Sph_Last_lAll> X = FPOP<Sph_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else  if (type_approx == 3) {
    FPOP<Rec_Last_lAll> X = FPOP<Rec_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else  if (type_approx == 4) {
    FPOP<RecSph_Last_lAll> X = FPOP<RecSph_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  }  else  if (type_approx == 5) {
    FPOP<Rec_All_lAll> X = FPOP<Rec_All_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else  if (type_approx == 6) {
    FPOP<RecSph_All_lAll> X = FPOP<RecSph_All_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  }
  return res;
}

//'@title TestTwoFPOPplus
//'
//' @description Сomparing the parameters ("UnpenalizedCost", "LastChpt") for two different FPOP methods (using the rectangle approximation of the sets) .
//' @param data is a matrix of data (p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param approximation1 is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection1 is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'(by default, 'last').
//' @param exclusion1 is the type of intersection : 'empty', 'all', 'random' or 'sphere'(by default, 'all').
//' @param approximation2 is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection2 is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'(by default, 'empty').
//' @param exclusion2 is the type of intersection : 'empty', 'all', 'random' or 'sphere'(by default, 'empty').
//'
//' @return TRUE or FALSE
//'
//' \describe{
//' \item{\code{TRUE}}{ 'TRUE' if parameters are the same.}
//' \item{\code{FALSE}}{'TRUE' if parameters are different.}
//' }
//'
//' @examples
//' N <- 1000
//' Dim <- 2
//' Penality <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
//' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'all', exclusion1 = 'all')
//' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'last', exclusion1 = 'all')
// [[Rcpp::export]]
bool TestTwoFPOPplus(Rcpp::NumericMatrix data, double penalty, std::string approximation1 = "rectangle", std::string intersection1 = "last",  std::string exclusion1 = "all", std::string approximation2 = "rectangle", std::string intersection2 = "empty",  std::string exclusion2 = "empty") {
  int type_approx1 = NmbOfapproxFPOP(approximation1, intersection1, exclusion1);
  int type_approx2 = NmbOfapproxFPOP(approximation2, intersection2, exclusion2);

  //----------stop--------------------------------------------------------------
  if (penalty < 0) {
    throw std::range_error("Penalty should be a non-negative number!");
    return false;
  }
  if((type_approx1 == 0)||(type_approx2 == 0)) {
    throw std::range_error("These combinations of parameters 'intersection',and 'exclusion' are  not available. ");
    return false;
  }
  if(type_approx1 == type_approx2) {
    throw std::range_error("These combinations have same parameters 'intersection'and 'exclusion'. ");
    return true;
  }
  //----------------------------------------------------------------------------
  bool res = false;
  if (type_approx1 == 1) {
    FPOP<Rec_Empty_Empty> X1 = FPOP<Rec_Empty_Empty>(data, penalty);
     X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else  if (type_approx1 == 2) {
    FPOP<Sph_Last_lAll> X1 = FPOP<Sph_Last_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 3) {
    FPOP<Rec_Last_lAll> X1 = FPOP<Rec_Last_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 4) {
    FPOP<RecSph_Last_lAll> X1 = FPOP<RecSph_Last_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 5) {
    FPOP<Rec_All_lAll> X1 = FPOP<Rec_All_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 6) {
    FPOP<RecSph_All_lAll> X1 = FPOP<RecSph_All_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  }
  return res;
}

