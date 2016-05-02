#' Find change points in a sequence of choices based on given test and decision criterion
#' @param input A vector of correct or incorrect decisions, coded as 1s and 0s
#' @param isDiscrete A boolean parameter, true if the data to be analyzed are discrete,
#' false if continuous, e.g. succesive intervals
#' @param test Name of the test to be performed. The following tests are implemented:
#' "binomial" (random rate), "chisquare", "ttest", "KS" (Komolgorov-Smirnov). The binomial test is
#' used for finding changes in the rate parameter of a random rate process, either intermittently
#' (isDiscrete = TRUE) or continuously (isDiscrete = FALSE) sampled. The chi square test is used with
#' data where there has been a frequency change, as in, e.g. correct choices. The t test should be
#' used for normally distributed data. The Komolgorov-Smirnof test is used for data with nonnormal or
#' unknown distribution.
#' @param Crit A real-valued decision criterion on logit. Values recommended by Gallistel et al. (2004)
#' are between 1.3 and 6, which correspond to p values of 0.05 and 0.000001, respectively.
#' The values are the logarithms of the odds against the null (no-change) hypothesis.
#' logit = log[(1-p)/p], where p is the desired significance level
#' @return A data frame with the number of Trials (Trial) for each change point detected,
#' the cumulative number of correct responses, or successes, (CumSs),
#' and the ratio of correct/total choices (Slopes) for the trials since the last change point
#' and the current change point
#' @references Gallistel CR, Fairhurst S, Balsam P (2004) The learning curve:
#' Implications of a quantitative analysis. PNAS 101:13124-13131. doi: 10.1073/pnas.0404965101
#' @examples
#' # generate dummy data with a change point at trial 120
#' input <- c(rbinom(120,1,0.5),rbinom(200,1,0.8))
#' cp_wrapper(input, TRUE, "chisquare", 3)
#' cp_wrapper(input, TRUE, "binomial", 3)
#'
#' # use included eyeblink data set:
#' eyeblink[,] # inspect data set
#' cp_wrapper(eyeblink, TRUE, "chisquare", 4)
#' # using small criterion, e.g. 2, fails due to computational problems
#' # better use "binomial" instead
#' cp_wrapper(eyeblink, TRUE, "binomial", 2)
#'
#' # use included plusmaze data set:
#' plusmaze[,] # inspect data set
#' (cp.1 <- cp_wrapper(plusmaze, TRUE, "binomial", 1.3))
#' # decrease sensitivity and detect different change points
#' (cp.2 <- cp_wrapper(plusmaze, TRUE, "binomial", 2))
#' # the chisquare test detects more change points even with higher criterion
#' (cp.3 <- cp_wrapper(plusmaze, TRUE, "chisquare", 2))
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
#' # plot average response rate per trial
#' ggplot() + geom_step(data=cp.1, aes(Trial,Slopes)) +
#'  ylab("Average Response Rate per Trial")
#'
#' @export

cp_wrapper <- function(input, isDiscrete, test, Crit)
{
  if(!is.logical(isDiscrete)) stop("Second input (isDiscrete) must be boolean")

  input <- unlist(input, use.names = FALSE)
  Cum <- unlist(cumsum(input),use.names = FALSE)


  if ((test=="binomial")&(!identical(input %% 1,mat.or.vec(length(input),nc=1)))&(isDiscrete))
    stop("When the binomial test is used with discrete-trial data,
         the data must be integer valued")

  if (isDiscrete==0 & test=="chisquare")
    stop("Cannot use chi square test when data are successive intervals")

  if ((test=="chisquare")&(sum(input + !input)!=length(input)))
    stop("When chi square test is used, data must be binary, i.e. 0 or 1")

  switch(test,
         binomial = {
           CritLength <- 1 # When binomial test is used, there must be at least two data
         },
         chisquare = {
           CritLength <- 7 # A test of differences of frequency cannot be significant when
           # the total number of observations is less than 8
         },
         KS = {
           CritLength <- 7 # K-S test is not valid when there are fewer than 4 data in either
           # of the two samples
         },
         ttest = {
           CritLength <- 2 # when t test is used there must be at least 3 data
         },
         stop("Test can only be 'binomial', 'chisquare', 'ttest', or 'KS' ")
  )

  outputlist <- c() # resetting the output

  if (test=="chisquare"|test=="binomial")
  {
    CP <- matrix(c(0,0),ncol=2)
    outputlist <- list(Cumt = Cum) # Initializiing for while loop. The Cumt vector will
    # be truncated as change points are found
    outputlist$r <- 1 # Initializing for while loop

    while (!is.null(outputlist$r)&(length(outputlist$Cumt)>CritLength))
    {
      if (!isDiscrete){ # data are continuous
        R <- cpc(outputlist$Cumt) # putative inflection points
        L <- rrc(outputlist$Cumt, R) # logit vector for continuous case
      } else { # data are discrete
        R <- cpd(outputlist$Cumt) # putative inflection points

        switch(test,
               binomial = { # if binomial test is to be used
                 L <- rrd(outputlist$Cumt,R) # logit vector
               },
               chisquare={ # if chisquare test is to be used
                 L <- chi2logit(outputlist$Cumt,R) # logit vector
               },
               stop("For discrete data test can only be 'binomial' or 'chisquare' ")
        )
      }
      outputlist <- trun(outputlist$Cumt,R,L,Crit) #Cumt is the truncated cumulative
      # record; Lt is the logit vector up to the point of
      # truncation (not used); r is the change point; r is empty if there is
      # no significant change point

      if (!is.null(outputlist$r)){ # if there is a change point, update change-point array
        if(!isDiscrete) # In the continuous case, the row count goes in the
          # y-column of the output array (the event count); in all other cases, it goes
          # in the x column. In the continuous case, the x column
          # contains the successive event times
        {
          cumrCP <- CP[length(CP[,1]),2] + outputlist$r # Add Cumt row for latest change point
          # to last change point to get Cum row of latest change point
          CP <- rbind(CP,c(Cum[cumrCP], cumrCP)) # Value of cumulative record at the change point
        }else{# In the discrete case, the row data go in the first column of CP
          cumrCP <- CP[length(CP[,1]),1] + outputlist$r # Add Cumt row for latest change point
          # to last change point to get Cum row of latest change point
          CP <- rbind(CP,c(cumrCP, Cum[cumrCP])) # Value of cumulative record at the change point
        }

      } # end of updating change-point array

    } # end of while loop for finding successive change points when the binomial test is used
  } # end of section that computes change-point array when binomial or chi square test is used

  if (test=="ttest"|test=="KS")
  {
    outputlist <- list(newinput = input) # Initializing for while loop. These vectors will be
    # truncated as CPs are found out
    CP <- matrix(c(0,0),ncol=2) # Initializing for while loop
    outputlist$r <- 1 # Initializing for while loop

    while (!is.null(outputlist$r)&(length(outputlist$newinput)>CritLength))
    {
      outputlist$newCum <- cumsum(outputlist$newinput)
      if (!isDiscrete){ # data are continuous
        R <- cpc(outputlist$newCum) # putative inflection points
      } else { # data are discrete
        R <- cpd(outputlist$newCum) # putative inflection points
      } #computing R vector

      switch(test,
             KS = { # if K-S test is to be used
               outputlist$r <- KS(outputlist$newinput,R, Crit) # logit vector
             },
             ttest={ # if t test is to be used
               outputlist$r <- cpt(outputlist$newinput,R, Crit) # logit vector
             }
      ) # end of computing new changepoint
      if (!is.null(outputlist$r)){ # if there is a change point, update change-point array and
        # truncate newinput
        cumrCP <- CP[length(CP[,1]),1] + outputlist$r # Add Cumt row for latest change point
        # to last change point to get Cum row of latest change point
        CP <- rbind(CP,c(cumrCP, Cum[cumrCP])) # Value of cumulative record at the change point
        outputlist$newinput <-
          outputlist$newinput[(outputlist$r+1):length(outputlist$newinput)] # Truncated data vector

      } # end of updating change-point array & truncating
    } # end of while loop for computing CP array when K-S or t test are used
  } # end of section that computes CP array when K-S or t test are used

  # Adding final point to output array
  if(isDiscrete) #
  {
    CP <- as.data.frame(rbind(CP, c(length(Cum), Cum[length(Cum)]))) # last row of CP array
    # gives coordinates of final point in cumulative record
    names(CP) <- c("Trial", "CumSs")
    # last row of CP array gives coordinates of final point in cumulative record
    CP %>% mutate_(Slopes = ~(CumSs-lead(CumSs))/(Trial-lead(Trial)),
                   Slopes = ~ifelse(is.na(Slopes), lag(Slopes), Slopes))
    #depends on the dplyr package
  } else { # in continuous case, row count goes in y column
    CP <- as.data.frame(rbind(CP, c(Cum[length(Cum)], length(Cum)))) # last row of CP array
    # gives coordinates of final point in cumulative record
    names(CP) <- c("Time", "Events")
    # last row of CP array gives coordinates of final point in cumulative record
    CP %>% mutate_(Slopes = ~(Events-lead(Events))/(Time-lead(Time)),
                   Slopes = ~ifelse(is.na(Slopes), lag(Slopes), Slopes))
  }
}
