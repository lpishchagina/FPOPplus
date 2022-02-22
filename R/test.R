N <- 1000
Noise <- 1
Dim <- 2
Penalty <- 2*Dim*log(N)
time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)


NbCand <- list()

PELT <-FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'empty', exclusion = 'empty',  NbOfCands = TRUE)
NbCand[[1]]<-PELT$NumberOfCandidats

PELTplus <-FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere', intersection = 'last', exclusion = 'all',  NbOfCands = TRUE)
NbCand[[2]]<-PELTplus$NumberOfCandidats

LastFPOP <-FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'last', exclusion = 'all',  NbOfCands = TRUE)
NbCand[[3]]<-LastFPOP$NumberOfCandidats

PELTplusLastFPOP <-FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere_rectangle', intersection = 'last', exclusion = 'all',  NbOfCands = TRUE)
NbCand[[4]]<-PELTplusLastFPOP$NumberOfCandidats

FPOP <-FPOPplus(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'all', exclusion = 'all',  NbOfCands = TRUE)
NbCand[[5]]<-FPOP$NumberOfCandidats

PELTplusFPOP <-FPOPplus(data = time_series, penalty = Penalty, approximation = 'sphere_rectangle', intersection = 'all', exclusion = 'all',  NbOfCands = TRUE)
NbCand[[6]]<-PELTplusFPOP$NumberOfCandidats


#NbCand

#NbCand[[1]]-NbCand[[5]]#PELT - FPOP

#NbCand[[1]]-NbCand[[2]]#PELT - PELT+

#NbCand[[2]]-NbCand[[5]]#PELT+ -FPOP

#NbCand[[2]]-NbCand[[3]]#PELT+ -LastFPOP

#NbCand[[3]]-NbCand[[4]]#LastFPOP - PELT+LastFPOP

#NbCand[[5]]-NbCand[[6]]#FPOP - PELT+FPOP

#NbCand[[4]]-NbCand[[6]]#PELT+LastFPOP - PELT+FPOP


