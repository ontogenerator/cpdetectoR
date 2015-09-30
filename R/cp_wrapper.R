#' Find change points in a sequence of choices based on given test and decision criterion
#' @param input A vector of correct or incorrect decisions, coded as 1s and 0s
#' @param test Name of the test to be performed. Currently the following tests are implemented: "binomial", "chisquare".
#' @param Crit A real-valued decision criterion on logit. Values recommended by Gallistel et al. (2004)
#' are between 1.3 and 6, which correspond to p values of 0.05 and 0.000001, respectively.
#' The values are the logarithms of the odds against the null (no-change) hypothesis.
#' logit = log[(1-p)/p], where p is the desired significance level
#' @return A data frame with the number of Trials (Trial) for each change point detected,
#' the cumulative number of correct responses, or successes, (CumSs), and the ratio of correct/total choices (Slopes)
#' for the trials since the last change point and the current change point.
#' @examples
#' # generate dummy data with a change point at trial 120
#' input <- c(rbinom(120,1,0.5),rbinom(200,1,0.8))
#' cp_wrapper(input, "chisquare", 3)
#' cp_wrapper(input, "binomial", 3)
#'
#' # use included eyeblink data set:
#' eyeblink[,] # inspect data set
#' cp_wrapper(eyeblink, "chisquare", 4)
#' # using small criterion, e.g. 2, fails due to computational problems
#' # better use "binomial" instead
#' cp_wrapper(eyeblink, "binomial", 2)
#'
#' # use included plusmaze data set:
#' plusmaze[,] # inspect data set
#' (cp.1 <- cp_wrapper(plusmaze, "binomial", 1.3))
#' # decrease sensitivity and detect different change points
#' (cp.2 <- cp_wrapper(plusmaze, "binomial", 2))
#' # the chisquare test detects more change points even with higher criterion
#' (cp.3 <- cp_wrapper(plusmaze, "chisquare", 2))
#'
#' # plotting data with ggplot2
#' library(ggplot2)
#' # first generate data.frame with cumulative responses:
#' my.data <- data.frame(Trial = 1:length(plusmaze[,]),
#' CumRespMeasure = cumsum(plusmaze)[,])
#' # plot cumulative learning curve with change points
#' ggplot(my.data) + geom_line(aes(Trial,CumRespMeasure)) +
#'  geom_point(data=cp.1, aes(Trial, CumSs), shape=1, size = 6, color = "blue") +
#'  geom_point(data=cp.2, aes(Trial, CumSs), shape=2, size = 3, color = "green") +
#'  geom_point(data=cp.3, aes(Trial, CumSs), shape=3, size = 3, color = "red")
#'
#' # plot average response rate per trial, with dplyr::lead
#' ggplot() + geom_step(data=cp.1, aes(Trial,dplyr::lead(Slopes))) +
#'  ylab("Average Response Rate per Trial") # ignore Warning message
#'
#' @export

cp_wrapper <- function(input, test, Crit)
{
  input <- unlist(input, use.names = FALSE)
  Cum <- unlist(cumsum(input),use.names = FALSE)


  if ((test=="binomial")&(!identical(input %% 1,mat.or.vec(length(input),nc=1))))
    stop("When the binomial test is used with discrete-trial data,
         the data must be integer valued")

  if ((test=="chisquare")&(sum(input + !input)!=length(input)))
    stop("When chi square test is used, data must be binary, i.e. 0 or 1")

  outputlist <- c()
  outputlist <- list(Cumt = Cum) # Initializiing for while loop. The Cumt vector will
  # be truncated as change points are found

  switch(test,
         binomial = {
           CritLength <- 1 # When binomial test is used, there must be at least two data
         },
         chisquare = {
           CritLength <- 7 # A test of differences of frequency cannot be significant when
           # the total number of observations is less than 8
         },
         stop("Test can only be 'binomial' or 'chisquare' ")
  )

  CP <- matrix(c(0,0),ncol=2)
  outputlist$r <- 1 # Initializing for while loop

  while (!is.null(outputlist$r)&(length(outputlist$Cumt)>CritLength))
  {
    R <- cpd(outputlist$Cumt) # putative inflection points

    switch(test,
           binomial = { # if binomial test is to be used
             L <- rrd(outputlist$Cumt,R) # logit vector
           },
           chisquare={ # if chisquare test is to be used
             L <- chi2logit(outputlist$Cumt,R) # logit vector
           },
           stop("Test can only be 'binomial' or 'chisquare' ")
    )

    outputlist <- trun(outputlist$Cumt,R,L,Crit) #Cumt is the truncated cumulative
    # record; Lt is the logit vector up to the point of
    # truncation (not used); r is the change point; r is empty if there is
    # no significant change point

    if (!is.null(outputlist$r)){ # if there is a change point, update change-point array
      # In the discrete case, the row data go in the first column of CP
      cumrCP <- CP[length(CP[,1]),1]+outputlist$r # Add Cumt row for latest change point
      # to last change point to get Cum row of latest change point
      CP <- rbind(CP,c(cumrCP, Cum[cumrCP])) # Value of cumulative
      #record at the change point
    } # end of updating change-point array

  } # end of while loop for finding successive change points when the binomial test is used

  # Adding final point to output array
  CP <- as.data.frame(rbind(CP, c(length(Cum), Cum[length(Cum)])))
  names(CP) <- c("Trial", "CumSs")
  # last row of CP array gives coordinates of final point in cumulative record
  CP %>% mutate_(Slopes = ~(CumSs-lag(CumSs))/(Trial-lag(Trial)))
  #depends on the dplyr package
}
