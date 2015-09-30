#' Compute the pseudologit vector for random rate discrete data, using binomial test
#' @param Cum Input vector of cumulative response counts
#' @param R A vector with trial numbers of putative change points
#' @return L A vector giving for each trial the pseudologit, approximately
#' the log of the odds that there has been a change point
#' @details Not normally called directly, but via the cp_wrapper function instead

rrd <- function(Cum, R)
{

  # L is the vector giving for each trial the pseudologit--log[P/(1-P+p)]--
  # for the probability of observing n_a or fewer post-CP responses if lambda_a/lambda_b = 1
  # where lambda_a and lambda_b are the rates of responding before and after
  # the putative change point; P is the probability of observing Cum-Cum(R) or fewer
  # responses after the putative change point; and p is the probability of
  # observing exactly Cum-Cum[R] responses after the putative change point.
  # Note that the numerator and denominator of the pseudologit are not complementary
  # probabilities, as they are in a true logit; they both include the probability of observing
  # exactly Cum-Cum[R] post-CP responses,
  # where Cum is the cumulative number of responsess at the moment of calculation.
  # The putative change point is the row (R=trial number) at which
  # the difference between the observed and expected number of responses is maximal
  # Cum/N is the slope of the cumulative record up to N (i.e., the average response rate
  # up to the trial of calculation
  # (Cum/N)*R is the expected number of responses up to trial R; Cum[R] is the observed number.
  # A negative pseudologit means that n_a is less than expected
  Trial <- 1:length(Cum) # The trial count vector
  Ta <- Trial-R # the trials since the putative change point
  p <- Ta/Trial # probability of any one response occuring during the trials since the putative CP
  Na <- Cum-Cum[R] # number of responses since putative CP

  Peqorl <- pbinom(Na,Cum,p)[-1]
  # Probability of observing Na or fewer total responses on the trials since putative CP

  Peqorm <- 1 - Peqorl + dbinom(Na,Cum,p)[-1] #Probability of observing Na or more
  # total responses on the trials since putative CP. Note that this probability overlaps
  # the "complementary" probability; both include the probability of observing exactly Na events.

  c(0,log10(Peqorl/Peqorm)) # Vector of the pseudologits


}


#' Compute the pseudologit vector for random discrete data using Chi squared and Fisher exact tests
#' @param Cum Input vector of cumulative response counts
#' @param R A vector with trial numbers of putative change points
#' @return Lgt A vector giving for each trial the pseudologit, approximately
#' the log of the odds that there has been a change point
#' @details Not normally called directly, but via the cp_wrapper function instead. For cases, in which
#' the chi square test is not valid (expected obersavations all smaller than 5), the Fisher exact test
#' is used instead. Occasionally, for noisy data even this test will fail. In that case use the
#' rrd function instead.

chi2logit <- function(Cum, R)
{
  # Level is the point beyond which the
  # Chi square test is not valid (because no expectation < 5)
  # and NV gives the rows for which the chi square test cannot
  # validly be performed (because at least one cell has an expectation
  # less than 5).
  N<-1:length(Cum) # The trial count vector


  Level1<- min(which(N>6)) # Finds the row
  # below which even the Fisher exact test should not be applied, because it
  # cannot yield a p value lower than .1 with fewer than 7 observations, and,
  # moreover, for fewer than 4 observations, fishexct returns p values
  # greater than 1.

  Ta <- N-R # the trials since the putative change point
  Onespre <- Cum[R] # number of (corect) responses up to the putative CP
  Onespre[which(is.na(Onespre))] <- 0
  Onespost <- Cum-Cum[R] # number of (correct) responses since the putative CP
  Onespost[which(is.na(Onespost))] <- 0
  Zeroespre <- R - Onespre # number of incorrect responses (or no responses) up to the putative CP
  Zeroespre[which(is.na(Zeroespre))] <- 0
  Zeroespost <- Ta - Onespost # number of incorrect responses (or no responses) since the putative CP
  Zeroespost[which(is.na(Zeroespost))] <- 0

  Smallest <- pmin(Cum, N-Cum)*pmin(R,Ta)/N # The smallest expectation is the
  # smallest row total times the smallest column total, divided by N

  NV <- which(Smallest<5) # Critical value for cell expectations is 5. NV is the
  # vector of rows whose contingency tables have a cell with less than the
  # critical value

  Level<- max(NV) # The highest level at which the  smallest expectation is less than 5.
  # Chi square is not valid when smallest expectation is less than 5.

  Output <- data.frame(N=N,Onespre = Onespre, Onespost = Onespost,
                       Zeroespre = Zeroespre, Zeroespost = Zeroespost) %>% rowwise() %>%
    # a data frame for the output is created, with each row representing a contingency table
    # for the subsequent tests. Further calculations depend on the dplyr package.
    mutate_(pchisq =
              ~ifelse(N > Level, chisq.test(cbind(c(Onespre, Onespost),
                                                  c(Zeroespre, Zeroespost)),
                                            simulate.p.value = TRUE, B = 1000)$p.value,1)) %>%
    # a column with the (simulated) p values from chi square test is added
    mutate_(pchisq =
             ~ifelse(N <= Level & N > Level1,fisher.test(
               cbind(c(Onespre, Onespost),c(Zeroespre, Zeroespost)),
               simulate.p.value = TRUE, B = 1000)$p.value,pchisq),
           Lgt = ~ifelse(pchisq == 1,0,log10(pchisq*(1-pchisq))))
  # For the cases in which chi square is invalid, calculate p values from Fisher's exact test
  # Zeros logit values for the initial string of observations within which
  # there are too few observations to do even Fisher's exact test
  unlist(Output %>% select_("Lgt"), use.names = FALSE)
}
