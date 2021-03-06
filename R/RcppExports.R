# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title FPOPplus
#'
#' @description FPOP method (using the rectangle approximation of the sets) for  the multiple changepoint detection.
#' @param data is a matrix of data (p-rows x n-columns).
#' @param penalty is a value of penalty (a non-negative real number).
#' @param approximation is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
#' @param intersection is the type of intersection : 'empty', 'all', 'last' or 'random' (by default, 'last').
#' @param exclusion is the type of intersection : 'empty', 'all'or 'random'(by default, 'all').
#' @param NbOfCands is the logical parameter (if NbOfCands = TRUE, than the file "NbOfCands.txt" contains the number of change candidates for each iteration.
#'
#' @return a list of  elements  = (changes, means, UnpenalizedCost, NumberOfCandidates).
#'
#' \describe{
#' \item{\code{changes}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
#' \item{\code{means}}{is the list of successive means for the p-variate time series.}
#' \item{\code{UnpenalizedCost}}{is a number equal to the global cost.}
#' \item{\code{NumberOfCandidates}}{is a number of candidates at each iteration (vector).}
#' }
#'
#' @examples
#' N <- 1000
#' Dim <- 2
#' Penalty <- 2*Dim*log(N)
#' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere', intersection = 'last', exclusion = 'all',NbOfCands = TRUE)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'last', exclusion = 'all',NbOfCands = TRUE)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere_rectangle', intersection = 'last', exclusion = 'all',NbOfCands = TRUE)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'all', exclusion = 'all',NbOfCands = TRUE)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere_rectangle', intersection = 'all', exclusion = 'all', NbOfCands = TRUE)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'empty', exclusion = 'empty', NbOfCands = TRUE)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere_rectangle', intersection = 'allInt', exclusion = 'all', NbOfCands = TRUE)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere_rectangle', intersection = 'lastInt', exclusion = 'all', NbOfCands = TRUE)
#' FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere', intersection = 'lastInt', exclusion = 'all', NbOfCands = TRUE)
NULL

FPOPplus <- function(data, penalty, approximation = "rectangle", intersection = "all", exclusion = "all", NbOfCands = FALSE) {
    .Call(`_FPOPplus_FPOPplus`, data, penalty, approximation, intersection, exclusion, NbOfCands)
}

#'@title TestTwoFPOPplus
#'
#' @description Сomparing the parameters ("UnpenalizedCost", "LastChpt") for two different FPOP methods (using the rectangle approximation of the sets) .
#' @param data is a matrix of data (p-rows x n-columns).
#' @param penalty is a value of penalty (a non-negative real number).
#' @param approximation1 is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
#' @param intersection1 is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'(by default, 'last').
#' @param exclusion1 is the type of intersection : 'empty', 'all', 'random' or 'sphere'(by default, 'all').
#' @param approximation2 is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
#' @param intersection2 is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'(by default, 'empty').
#' @param exclusion2 is the type of intersection : 'empty', 'all', 'random' or 'sphere'(by default, 'empty').
#'
#' @return TRUE or FALSE
#'
#' \describe{
#' \item{\code{TRUE}}{ 'TRUE' if parameters are the same.}
#' \item{\code{FALSE}}{'TRUE' if parameters are different.}
#' }
#'
#' @examples
#' N <- 1000
#' Dim <- 2
#' Penality <- 2*Dim*log(N)
#' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
#' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'all', exclusion1 = 'all')
#' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'last', exclusion1 = 'all')
TestTwoFPOPplus <- function(data, penalty, approximation1 = "rectangle", intersection1 = "last", exclusion1 = "all", approximation2 = "rectangle", intersection2 = "empty", exclusion2 = "empty") {
    .Call(`_FPOPplus_TestTwoFPOPplus`, data, penalty, approximation1, intersection1, exclusion1, approximation2, intersection2, exclusion2)
}

