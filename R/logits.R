#' Compute the pseudologit vector for random rate discrete data, using binomial test
#' @param Cum Input vector of cumulative response counts
#' @param R A vector with trial numbers of putative change points
#' @return L A vector giving for each trial the pseudologit, approximately
#' the log of the odds that there has been a change point
#' @details Not normally called directly, but via the cp_wrapper function instead

rrd <- function(Cum, R)
{
  # L is the vector giving for each trial the pseudologit--log[P / (1 - P + p)]--
  # for the probability of observing n_a or fewer post-CP responses if
  # lambda_a / lambda_b = 1
  # where lambda_a and lambda_b are the rates of responding before and after
  # the putative change point; P is the probability of observing Cum - Cum(R) or fewer
  # responses after the putative change point; and p is the probability of
  # observing exactly Cum - Cum[R] responses after the putative change point.
  # Note that the numerator and denominator of the pseudologit are not complementary
  # probabilities, as they are in a true logit; they both include the probability of
  # observing exactly Cum - Cum[R] post-CP responses, where Cum is the cumulative number
  # of responsess at the moment of calculation.
  # The putative change point is the row (R = trial number) at which the difference
  # between the observed and expected number of responses is maximal.
  # Cum / N is the slope of the cumulative record up to N (i.e., the average response rate
  # up to the trial of calculation
  # (Cum / N) * R is the expected number of responses up to trial R; Cum[R] is the
  # observed number. A negative pseudologit means that n_a is less than expected.
  Trial <- 1:length(Cum) # The trial count vector
  Ta <- Trial - R # the trials since the putative change point
  p <- Ta / Trial # probability of any one response occuring during the trials since the
  # putative CP
  Na <- Cum - Cum[R] # number of responses since putative CP

  Peqorl <- stats::pbinom(Na, Cum, p)[-1]
  # Probability of observing Na or fewer total responses on the trials since putative CP

  Peqorm <- 1 - Peqorl + stats::dbinom(Na, Cum, p)[-1] #Probability of observing Na or
  # more total responses on the trials since putative CP. Note that this probability
  # overlaps the "complementary" probability; both include the probability of observing
  # exactly Na events.

  c(0, log10(Peqorl / Peqorm)) # Vector of the pseudologits

}

#' Compute the pseudologit vector for continuous random rate case
#' @param Cum A cumulative interevent interval vector
#' @param R A vector with trial numbers of putative change points
#' @return L A vector giving for each trial the pseudologit, approximately
#' the log of the odds that there has been a change point
#' @details  For use when finding changes in the rate parameter of a random rate process.
#' Not normally called directly, but via the cp_wrapper function instead

rrc <- function(Cum, R)
{
  # The output is the vector giving for each event the pseudologit log[P / (1 - P + p)]
  # for the probability of observing n_a or fewer post-CP responses if
  # lambda_a / lambda_b = 1
  # where lambda_a and lambda_b are the rates of responding before and after
  # the putative change point; P is the probability of observing N - R or fewer
  # events after the putative change point; and p is the probability of
  # observing exactly N - R events after the putative change point.
  # Note that the numerator and denominator of the pseudologit are not complementary
  # probabilities, as they are in a true logit; they both include the probability of
  # observing exactly N - R events, where N is the total number of events at the moment of
  # calculation. The putative change point is the row (R = event number) at which
  # the difference between (N / Cum[N]) * (Cum[R]) and Cum[R] is maximal.
  # N / Cum[N] is the slope of the cumulative record up to the Nth event.
  # Cum[R] is the time up to the Rth event.
  # A negative pseudologit means that n_a is less than expected.

  Ta <- Cum - Cum[R] # the interval elapsed since the putative CP
  p <- Ta / Cum # probability of any one event falling after the putative CP
  N <- 1:length(Cum) # event count vector
  Na <- N - R # vector giving number of events since putative CP



  Peqorl <- stats::pbinom(Na, N, p)[-1]
  # Probability of observing Na or fewer events in the interval Ta

  Peqorm <- 1 - Peqorl + stats::dbinom(Na, N, p)[-1] # Probability of observing Na or more
  # events in the interval Ta. Note that this probability overlaps the "complementary"
  # probability; both include the probability of observing exactly Na events.

  c(0, log10(Peqorl / Peqorm)) # Vector of the pseudologits
}


#' Compute the pseudologit vector for random discrete data using Chi squared and Fisher
#' exact tests
#' @param Cum Input vector of cumulative response counts
#' @param R A vector with trial numbers of putative change points
#' @return Lgt A vector giving for each trial the pseudologit, approximately
#' the log of the odds that there has been a change point
#' @details Not normally called directly, but via the cp_wrapper function instead.
#' For cases, in which the chi square test is not valid (expected obersavations all
#' smaller than 5), the Fisher exact test is used instead. Occasionally, for noisy data
#' even this test will fail. In that case use the rrd function instead. Notice that the p
#' values delivered by this test are estimated from a 1000 Monte Carlo simulations.

chi2logit <- function(Cum, R)
{
  # Level is the point beyond which the
  # Chi square test is not valid (because no expectation < 5)
  # and NV gives the rows for which the chi square test cannot
  # validly be performed (because at least one cell has an expectation
  # less than 5).
  N <- 1:length(Cum) # The trial count vector

  Level1 <- min(which(N > 7)) # Finds the row
  # below which even the Fisher exact test should not be applied, because it
  # cannot yield a p value lower than .1 with fewer than 7 observations, and,
  # moreover, for fewer than 4 observations, fishexct returns p values
  # greater than 1.

  Ta <- N - R # the trials since the putative change point
  Onespre <- Cum[R] # number of (corect) responses up to the putative CP
  Onespre[which(is.na(Onespre))] <- 0
  Onespost <- Cum - Cum[R] # number of (correct) responses since the putative CP
  Onespost[which(is.na(Onespost))] <- 0
  Zeroespre <- R - Onespre # number of incorrect responses (or no responses) up to the
  # putative CP
  Zeroespre[which(is.na(Zeroespre))] <- 0
  Zeroespost <- Ta - Onespost # number of incorrect responses (or no responses) since the
  # putative CP
  Zeroespost[which(is.na(Zeroespost))] <- 0

  Smallest <- pmin(Cum, N - Cum) * pmin(R, Ta) / N # The smallest expectation is the
  # smallest row total times the smallest column total, divided by N

  NV <- which(Smallest < 5) # Critical value for cell expectations is 5. NV is the
  # vector of rows whose contingency tables have a cell with less than the
  # critical value

  Level <- max(NV) # The highest level at which the  smallest expectation is less than 5.
  # Chi square is not valid when smallest expectation is less than 5.

  p_chisq <- rep(1, length(N)) # initiate pvalue vector
  # mat <- matrix(rep(0, 4), nrow = 2, ncol = 2)

  for (i in 1:length(N)) # loop over the trial count vector
  {
    mat <- cbind(c(Onespre[i], Onespost[i]),
                 c(Zeroespre[i], Zeroespost[i]))

    p_chisq[i] <-
      ifelse(N[i] > Level,
             stats::chisq.test(mat, simulate.p.value = TRUE, B = 1000)$p.value,
             ifelse(N[i] <= Level & N[i] > Level1,
                    # For the cases in which chi square is invalid, calculate p values
                    # from Fisher's exact test.
                    stats::fisher.test(mat, simulate.p.value = TRUE, B = 1000)$p.value,
                    1
                    )
             )
  }

  Lgt <- ifelse(p_chisq == 1, 0, log10(p_chisq / (1 - p_chisq)))
  Lgt
}


#' Uses t test to find first significant change point
#' @param Data A vector of trial by trial measures or successive intervals
#' @param R A vector of putative change points
#' @param Crit Decision criterion, the value the logit must exceed
#' for the function to return a significant change point
#' @return CP The first significant change point
#' @details This test is appropriate if one is looking for a change in the expectation of
#' a renewal event-generating process, where the interevent intervals are normally (rather
#' than exponentially) distributed. Not normally called directly, but via the cp_wrapper
#' function instead.

cpt <- function(Data, R, Crit)
{
  Data <- unlist(Data, use.names = FALSE)
  L <- c(0, 0) # initialization of the logit vector up to and including
  # the row where the significance criterion (Crit) was exceeded
  r <- 3 # Initializing for while loop, with index r

  while (abs(L[length(L)]) < Crit)
  { # loop that ends when critical L found or end of data reached
    if (!is.na(R[r])) # if R[r] is NA skip to next r
    {
      if ((sum(stats::sd(Data[1:R[r]]), stats::sd(Data[(R[r] + 1):r]),
               na.rm = TRUE) > 0) &
          (length(Data[1:R[r]]) > 1) & (length(Data[(R[r] + 1):r]) > 1))# test cannot be
        # run when there is no variance on either side of putative CP, for example, in the
        # sequence 0 0 0 3 3 3 or from a single observation e.g. 0 vs 3 2 5
      {
        pb <- stats::t.test(Data[1:R[r]], Data[(R[r] + 1):r], "greater")$p.value
        # Probability that the mean after the putative change point is greater than the
        # mean up to and including the putative change point
        pl <- stats::t.test(Data[1:R[r]], Data[(R[r] + 1):r], "less")$p.value
        # Probability that the mean after the putative change point is less than the
        # mean up to and including the putative change point
        L[r] <- log10(pb / pl) # Latest logit
      } else {
        L[r] <- 0
      } # end of if that computes individual logit (L) values
    } # end of if that checks if R[r] is NA


    r <- r + 1 # Incrementing latest row for next iteration
    if (r > length(Data)) break # end of data reached

  } # end of while loop

  if (abs(L[length(L)]) > Crit) CP <- R[r - 1] else CP <- NULL # if no significant change
  # point, CP empty
  CP
}


#' Uses Komolgorov-Smirnov test to find first significant change point
#' @param Data A vector of trial by trial measures or successive intervals
#' @param R A vector of putative change points
#' @param Crit Decision criterion, the value the logit must exceed
#' for the function to return a significant change point
#' @return r1 the change point row when the decision criterion is exceeded
#' @details Not normally called directly, but via the cp_wrapper function instead.

KS <- function(Data, R, Crit)
{
  # L is the pseudologit vector for rows where the
  # approximation formula for the Kolmogorov-Smirnov p is valid.
  # Its final value is the first value to exceed the decision criterion
  # t is the col vector of rows for which L is defined
  # r1 is the change point row when the decision criterion is exceeded
  # r2 is the row at which the decision criterion is exceeded
  # If there are no testable rows or if the decision criterion is never
  # exceeded, the variables are returned empty

  N <- 1:length(Data) # Number of rows in Data vector
  L <- rep(0, length(N)) # initialization of the logit vector

  Na <- N - R # Col vector giving for each row in Data the number of rows
  # after the putative change point. So R gives the number of
  # rows before the change point and Na the number after

  Test <- which(Na * R / (Na + R) >= 4) # For the approximation to the KS
  # probability to be valid, the product of the two n's divided by
  # the sum must be greater than or equal to 4. Test is the col
  # vector of rows that satisfy this constraint--the row numbers(!),
  # not the entries themselves
  r1 <- c()
  r2 <- c()
  t <- rep(0, length(Test))

  if (length(Test) > 0) {# check if any rows satisfy the constraint and proceed
    for (T in 1:length(Test)) {
      pval <- suppressWarnings(stats::ks.test(Data[1:R[Test[T]]],
                                       Data[(R[Test[T]] + 1):Test[T]],
                                       exact = FALSE)$p.value)
      L[T] <- log10((1 - pval) / pval)
      t[T] <- Test[T] # The row to which the latest value of L[T] "belongs"

      if (L[T] > Crit) # Value of logit exceeds decision criterion
      {
        r2 <- Test[T] # the row (in Data) at which the criterion is exceeded
        r1 <- R[Test[T]] # the change point when the criterion is exceeded
        L <- L[1:T] # Drop rows of L after row in which decision criterion reached
        break # Break out of loop when decision criterion exceeded
      }
    }
  }
  r1
}
