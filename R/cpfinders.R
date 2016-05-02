#' Find putative change point in discrete-time cumulative records
#' @param Cum Input vector of cumulative response counts
#' @return R A vector with the same length as Cum, with putative change points for each
#' trial. These correspond to the preceding trial at which the deviation of the
#' observed count from the expected count is maximal
#' @details Not normally called directly, but via the cp_wrapper function instead.

cpd <- function(Cum)
{
  N <- 1:length(Cum) #Trial count vector
  Slopes <- Cum/N #Average count or measure per trial for trials 1 to N

  M <- diag(length(Cum))
  M[upper.tri(M)] <- 1 # Mask with ones on and above diagonal & zeros below

  Diag <- M * matrix(rep(Slopes, length(Slopes)), nrow = length(Slopes), byrow = TRUE)
  # Creates an array in which successive cols have successive slopes of the cumulative record.
  # The slope for a given col fills all the cells on and above the main diagonal

  Preds <- 1:length(Slopes) * Diag
  # Predicted (expected) cumulative values in a diagonal array

  Obs <- Cum * M
  # Diagonal array of observed cumulations

  Devs<- abs(Obs-Preds) # Diagonal array of deviations from expectations

  mx <- M * matrix(rep((apply(Devs,2,max)), length(Slopes)), nrow = length(Slopes), byrow = TRUE)

  # mx is a matrix listing the maxima of the deviations

  R <- left_join(data.frame(col=1:length(Cum)),
                 as.data.frame(which(Devs == mx & mx > 0, arr.ind=TRUE)),by="col") %>%
    rename_(R = "row", N = "col")
  # R at this stage is a data.frame with columns N, and R. R has the trial numbers of the
  # putative change points and N - the trial numbers. This and further data manipulations
  # depend on the dplyr package

  unlist(R %>% group_by(N) %>% summarise_(R = ~min(R)) %>% select_("R"), use.names = FALSE)
  #return only R as a vector
}

#' Find putative change point in continuous-time cumulative records
#' @param Cum Input vector of the cumulative interevent intervals
#' @return R A vector with the same length as Cum, with putative change points for each
#' event. The putative change point corresponding to the Nth
#' event is the preceding event at which the deviation of the observed event count
#' from the expected event count is maximal. The expected event count
#' at any earlier event, n, is Cum[n], the interval up to the nth event,
#' divided by the average interevent interval over the range from n = 0 to n = N
#' The deviation from expectation is n - this expectation. R is the value of n at
#' which this deviation is maximal
#' @details Not normally called directly, but via the cp_wrapper function instead

cpc <- function(Cum)
{
  N <- 1:length(Cum) # Event count vector
  Slopes <- N/Cum # Average slope up to given point in cumulative function

  M <- diag(length(Cum))
  M[upper.tri(M)] <- 1 # Mask with ones on and above diagonal & zeros below

  Diagonal <- M * matrix(rep(Slopes, length(Slopes)), nrow = length(Slopes), byrow = TRUE)
  # Creates an array in which successive cols have successive slopes of the cumulative record.
  # The slope for a given col fills all the cells on and above the main diagonal

  Preds <-  Cum * Diagonal
  # Creates diagonal array of the predicted numbers of events at each time in Cum

  Obs <- 1:length(Slopes) * M
  # Diagonal array with actual numbers of events

  Devs<- abs(Obs-Preds) # Diagonal array of deviations from expectations

  mx <- M * matrix(rep((apply(Devs,2,max)), length(Slopes)), nrow = length(Slopes), byrow = TRUE)

  # mx is a matrix listing the maxima of the deviations

  R <- left_join(data.frame(col=1:length(Cum)),
                 as.data.frame(which(Devs == mx & mx > 0, arr.ind=TRUE)),by="col") %>%
    rename_(R = "row", N = "col")
  # R at this stage is a data.frame with columns N, and R. R has the cumulative inter-event intervals
  # of the putative change points and N - the event numbers. This and further data manipulations
  # depend on the dplyr package.

  unlist(R %>% group_by(N) %>% summarise_(R = ~min(R)) %>% select_("R"), use.names = FALSE)
  #return only R as a vector
}

#' Truncate the cumulative record at significant change point
#' @param Cum Input vector of cumulative response counts
#' @param R A vector with trial numbers of putative change points
#' @param L A vector giving for each trial the pseudologit, approximately
#' the log of the odds that there has been a change point
#' @param Crit A real-valued decision criterion on logit
#' @return A list of Cumt, the truncated cumulative record,
#' Lt, the truncated L vector up to the row at which a change point was detected,
#' and r, the row at which L was truncated.
#' @details Not normally called directly, but via the cp_wrapper function instead.
#' Works for both discrete and continuous cumulative records. All input arguments obligatory.

trun <- function(Cum, R, L, Crit)
{
  # Cumt is truncated cumulative record,
  # the record as it would be if observation began at time CP+ or after trial CP
  # Lt is the L vector truncated at Alert, which is the row at which an CP was detected
  # r is row at which it was truncated.

  La <- abs(L)
  if(length(which(La>Crit))>0) {
    Alert <- min(which(La>Crit))# Finds first row where decision criterion is exceeded
  } else { # If there is no such row, then Alert will be NULL
    Alert <- NULL
  }

  if (is.null(Alert)) {
    list(Cumt = Cum, Lt = L, r = Alert)
  }
  else
  {
    r <- R[Alert] # The putative change point at the value of Cum that first yields
    # significant logit.This putative change point is always the number of an earlier
    # event or trial

    I <- c(Cum[1],diff(Cum)) # Interevent interval vector OR vector
    # of responses on successive trials

    Cumt <- cumsum(I[(r+1):(length(I))]) # Truncated cumulative record

    Lt<- L[1:Alert] # Truncated logit vector

    list(Cumt = Cumt, Lt = Lt, r = r)
  }

}
