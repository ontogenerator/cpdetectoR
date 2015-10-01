Description
===========

This package is for finding change points in sequences of responses. Putative change points are first detected, then a statistical test at a predefined significance level (Criterion) is applied to decide if the change point is supported or not. Based on algorithm published by Gallistel et al. (2004).

Examples
========

There is a single wrapper function, called `cp_wrapper`, which can be used to analyze sequences of responses. Currently only binary responses (containing ones and zeroes) are supported. The second argument of the function gives the statistical test (either "binomial", or "chisquare"). The different tests usually give similar, but often different results. The third and last argument is the criterion, with values between 1.3 and 6.

The input can be either a vector:

``` r
library(cpdetectorr) # load package first
cp_wrapper(c(rbinom(50,1,0.3), rbinom(50,1,0.8)), "binomial",2)
```

or a dataframe:

``` r
d_responses <- data.frame(Responses = c(rbinom(50,1,0.3),
rbinom(50,1,0.8)))

cp_wrapper(d_responses, "chisquare",2)
```

For the same value of the criterion, the chi square test usually gives a higher number of change points (in this case, false positives). The value of the criterion should lie between 1.3 and 6, correpsonding to p values of 0.05 and 0.000001, respectively. These values are the logarithms of the odds against the null (no-change) hypothesis.

Let us look first at the included eyeblink data set:

``` r
eyeblink[,] # inspect data set
```

Gallistel et al. advise against using the chi square test on these data. Indeed, with a criterion of 2, the test fails:

``` r
cp_wrapper(eyeblink, "chisquare", 2)
```

However, using either a higher criterion or the "binomial" test gives a result:

``` r
cp_wrapper(eyeblink, "chisquare", 4)

cp_wrapper(eyeblink, "binomial", 2)
```

With ggplot we can visualize the results by first generating a `data.frame` with the cumulative responses:

``` r
library(ggplot2) # load the ggplot package
eyeblinkdata <- data.frame(Trial = 1:length(eyeblink[,]),
                           CumRespMeasure = cumsum(eyeblink)[,])
changepoints <- cp_wrapper(eyeblink, "binomial", 4)

ggplot(eyeblinkdata) + geom_line(aes(Trial, CumRespMeasure)) +
  geom_point(data = changepoints, aes(Trial, CumSs), size = 3)
```

Another type of plot one can look at is the average response rate per trial vs trial. First let us look at the plusmaze data set included with the package.

``` r
plusmaze[,] # inspect data set
(cp.1 <- cp_wrapper(plusmaze, "binomial", 1.3)) #find the change points
# plot average response rate per trial, with dplyr::lead
ggplot() + geom_step(data=cp.1, aes(Trial,dplyr::lead(Slopes))) +
  ylab("Average Response Rate per Trial") # ignore Warning message
```

References
==========

1.  Gallistel CR, Fairhurst S, Balsam P (2004) The learning curve: Implications of a quantitative analysis. PNAS 101:13124-13131. doi: 10.1073/pnas.0404965101