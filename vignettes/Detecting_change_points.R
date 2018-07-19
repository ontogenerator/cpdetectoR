## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(cpdetectoR) # load package first
set.seed(25)
cp_wrapper(c(rbinom(50, 1, 0.3), rbinom(50, 1, 0.8)), TRUE, "binomial", 2)

## ------------------------------------------------------------------------
d_responses <- data.frame(Responses = c(rbinom(50, 1, 0.3), rbinom(50, 1, 0.8)))
cp_wrapper(d_responses, TRUE, "chisquare", 2)

## ------------------------------------------------------------------------
eyeblink[,] # inspect data set

## ---- eval = FALSE-------------------------------------------------------
#  cp_wrapper(eyeblink, TRUE, "chisquare", 2)

## ------------------------------------------------------------------------
cp_wrapper(eyeblink, TRUE, "chisquare", 3)

cp_wrapper(eyeblink, TRUE, "binomial", 2)

## ------------------------------------------------------------------------
library(ggplot2) # load the ggplot package
eyeblinkdata <- data.frame(Trial = 1:length(eyeblink[,]),
                           CumRespMeasure = cumsum(eyeblink)[,])
changepoints <- cp_wrapper(eyeblink, TRUE, "binomial", 4) # save the output of the change point analysis
#generate a cumulative response vs trial plot:
ggplot(eyeblinkdata) + geom_line(aes(Trial, CumRespMeasure)) +
  geom_point(data = changepoints, aes(Trial, CumSs), size = 3)

## ------------------------------------------------------------------------
plusmaze[,] # inspect data set
(cp.1 <- cp_wrapper(plusmaze, TRUE, "binomial", 1.3)) #find the change points
# plot average response rate per trial
ggplot() + geom_step(data = cp.1, aes(Trial,Slopes)) +
  ylab("Average Response Rate per Trial")
# for comparison, the cumulative response vs trial plot, as in the example above:
plusmazedata <- data.frame(Trial = 1:length(plusmaze[,]),
                           CumRespMeasure = cumsum(plusmaze)[,])
ggplot(plusmazedata) + geom_line(aes(Trial, CumRespMeasure)) +
  geom_point(data = cp.1, aes(Trial, CumSs), size = 3)

## ------------------------------------------------------------------------
(cp.2 <- cp_wrapper(hopperentry, TRUE, "ttest", 4)) #find the change points
# cumulative response vs trial plot
hedata <- data.frame(Trial = 1:length(hopperentry[,]),
                           CumRespMeasure = cumsum(hopperentry)[,])
pl1 <- ggplotGrob(ggplot(hedata) + geom_line(aes(Trial, CumRespMeasure)) +
  geom_point(data = cp.2, aes(Trial, CumSs), size = 3))
# plot average response rate per trial
pl2 <- ggplotGrob(ggplot(cp.2) + geom_step(aes(Trial, Slopes)) +
  ylab("Average Response Rate per Trial"))
# stack the two plots vertically using the grid package
grid::grid.draw(rbind(pl1, pl2, size = "first"))

## ------------------------------------------------------------------------
(cp.3 <- cp_wrapper(matching, FALSE, "binomial", 2)) #find the change points
# cumulative response vs trial plot
matchingdata <- data.frame(Events = 1:length(matching[,]),
                           Time = cumsum(matching)[,])
pl3 <- ggplotGrob(ggplot(matchingdata) + geom_line(aes(Time, Events)) +
  geom_point(data = cp.3, aes(Time, Events), size = 3))
# plot average response rate per trial
pl4 <- ggplotGrob(ggplot(cp.3) + geom_step(aes(Time, Slopes)) +
  ylab("Events per Unit Time"))
grid::grid.draw(rbind(pl3, pl4, size = "first"))

