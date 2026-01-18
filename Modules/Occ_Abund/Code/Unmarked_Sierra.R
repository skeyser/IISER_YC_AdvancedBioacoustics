## -------------------------------------------------------------
##
## Script name: Fitting Occupancy Models using Unmarked
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2026-01-17
##
## Email: srk252@cornell.edu
##
## Github: https://github.com/skeyser
##
## -------------------------------------------------------------
##
## Notes:
##
##
## -------------------------------------------------------------

## Defaults
options(scipen = 10, digits = 10)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(here)
library(unmarked)
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Unmarked Model Prep
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load in the data we saved earlier
occ.dat <- readRDS(here("./Modules/Occ_Abund/Data/Occ_Data/UnmarkedData_SpeciesThresh_975minMaxPrex_NewVars.rds"))
dim(occ.dat$y)

## To reduce the amount of data let's subset to 500 random locations
X <- seq(1:nrow(occ.dat$siteCovs))
rs <- sort(sample(x = X, size = 500, replace = F))

## Break the data apart for demonstration purposes
y <- occ.dat$y 
siteCovs <- occ.dat$siteCovs
obsCovs <- occ.dat$obsCovs

## Subsampled sites
y <- y[,rs,]
siteCovs <- siteCovs[rs,]
obsCovs <- lapply(obsCovs, FUN = function(x) x[rs,])

# Take one species since we are doing a single species model
y <- y[1,,]

## Package these elements into an unmarked object
unmk <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
class(unmk)
summary(unmk)

## -------------------------------------------------------------
##
## Begin Section: Model Fitting
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Null model or Intercept-only model
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Start with the simplest possible model (i.e., a null model)
fm1 <- occu(~1 ~1, data = unmk)

## Take a minute to look at these values...what do they mean?
summary(fm1)

## Mean detection probability for this species is...
plogis(1.49) # 81.6% probability that we detect the species if it's present
plogis(-2.02) # The mean probability that a species occurs at a site is 11.7%

## We can get these value transformed returned from unmarked itself with 95% CIs
predict(fm1, type = "state")[1,]
predict(fm1, type = "det")[1,]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Detection covariates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fm2 <- occu(~ scale(eff.hrs) + scale(eff.jday) ~1, data = unmk)

## Before we were able to directly estimate backtransformed detection and
## occupancy because they were not a function of covariates...now we have
## to select some values since we have a detection model with covariates

## occupancy - note we didn't use covariates in this part of the model!!
predict(fm2, type = "state")[1,]

## detection
## make a new dataframe
nd <- data.frame(eff.hrs = 0, eff.jday = 0)
round(predict(fm2, type = "det", newdata = nd, appendData = TRUE), 2)

## Try this again but varying one of the predictors (also known as generate marginal or conditional response values)
nd <- data.frame(eff.hrs = 10:30, eff.jday = 0)
round(predict(fm2, type = "det", newdata = nd, appendData = TRUE), 2)

## We can get confidence intervals for parameter estimates via confint
confint(fm2, type = "det")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Add some occupancy covariates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fm3 <- occu(~ scale(eff.hrs) + scale(eff.jday) ~ scale(topo_elev), data = unmk)

## make a new dataframe
nd <- data.frame(eff.hrs = 0, eff.jday = 0, topo_elev = seq(min(siteCovs$topo_elev), max(siteCovs$topo_elev), 100))
round(predict(fm3, type = "state", newdata = nd, appendData = TRUE), 2)

## We can get confidence intervals for parameter estimates via confint
confint(fm3, type = "det")
confint(fm3, type = "state")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Model Selection and Fit
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Unmarked fitList()
fms <- fitList('psi(.)p(.)' = fm1, 
               'psi(.)p(Hours+Date)' = fm2,
               'psi(elevation)p(Hours+Date)' = fm3)
modSel(fms)

## What do you see? Which model is more appropriate?

## To assess model fit we can use a Chi-square discrepancy test since our data is binary
chisq <- function(fm){
  umf <- fm@data
  y <- umf@y
  y[y>1] <- 1
  fv <- fitted(fm)
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

set.seed(182)
pb <- parboot(fm3, statistic = chisq, nsim = 100, parallel = FALSE)
pb

## -------------------------------------------------------------
##
## Begin Section: Visualizations
##
## -------------------------------------------------------------

## Quick plotting via unmarked
plotEffects(fm3, "state", "topo_elev")
plotEffects(fm3, "det", "eff.hrs")
plotEffects(fm3, "det", "eff.jday")


plot.df <- plotEffectsData(fm3, "state", "topo_elev")

ele.plot <- ggplot(data = plot.df, aes(x = covariateValue, y = Predicted)) + 
  geom_line(size = 1, color = "red") +
  geom_line(aes(x = covariateValue, y = upper), linetype = "dashed", color = "black") +
  geom_line(aes(x = covariateValue, y = lower), linetype = "dashed", color = "black") +
  geom_ribbon(aes(ymax = upper, ymin = lower), fill = "gray", alpha = 0.5) + 
  theme_bw() + 
  ylab(expression(psi)) + 
  xlab("Elevation (masl)")
  
ele.plot

## Marginal effect of elevation on occupancy assuming mean detection
nd <- data.frame(eff.hrs = mean(obsCovs$eff.hrs, na.rm = T), eff.jday = mean(obsCovs$eff.jday, na.rm = T), topo_elev = seq(min(siteCovs$topo_elev), max(siteCovs$topo_elev), length.out = 200))
plot.df <- predict(fm3, type = "state", newdata = nd, appendData = TRUE)

ele.plot <- ggplot(data = plot.df, aes(x = topo_elev, y = Predicted)) + 
  geom_line(size = 1, color = "red") +
  geom_line(aes(x = topo_elev, y = upper), linetype = "dashed", color = "black") +
  geom_line(aes(x = topo_elev, y = lower), linetype = "dashed", color = "black") +
  geom_ribbon(aes(ymax = upper, ymin = lower), fill = "gray", alpha = 0.5) + 
  theme_bw() + 
  ylab(expression(psi)) + 
  xlab("Elevation (masl)")

## -------------------------------------------------------------
##
## Begin Section: Single-season, multi-species version
##
## -------------------------------------------------------------
## Load in the data we saved earlier
occ.dat <- readRDS(here("./Modules/Occ_Abund/Data/Occ_Data/UnmarkedData_SpeciesThresh_975minMaxPrex_NewVars.rds"))
dimnames(occ.dat$y)
spoi <- c("Acorn Woodpecker", 
          "Hermit Warbler", 
          "Golden-crowned Kinglet", 
          "Lazuli Bunting", 
          "Olive-sided Flycatcher")
sp.ind <- which(dimnames(occ.dat$y)[[1]] %in% spoi)

## Break the data apart for demonstration purposes
dimnames(occ.dat$y)
y <- occ.dat$y[sp.ind,,] # Take species of interest

## y is expected to be formatted M x J x S (Sites x Reps x Species)
## We can use aperm to reorder
y.perm <- aperm(y, c(2, 3, 1))

## Site and obs covariates are the same for all species
siteCovs <- occ.dat$siteCovs
obsCovs <- occ.dat$obsCovs

## Package these elements into an unmarked object
unmk.comm <- unmarkedFrameOccuComm(y = y.perm, siteCovs = siteCovs, obsCovs = obsCovs, speciesCovs = NULL)
class(unmk.comm)

## Let's test the same best fitting model from above
fm4 <- occuComm(~ scale(eff.hrs) + scale(eff.jday) ~ scale(topo_elev), data = unmk.comm)

## Inspect the model
fm4

## We can still do AIC model selection for the multi-species models
## Add a polynomial term for the elevation in a slightly more complex model form
fm5 <- occuComm(~ scale(eff.hrs) + scale(eff.jday) ~ scale(topo_elev) + scale(I(topo_elev^2)), data = unmk.comm)

fm4@AIC > fm5@AIC
(fm4@AIC - fm5@AIC) > 4

## We can extract the AICs and calculate AIC weights with something like this
AICweight <- function(models){
  if(is.null(names(models))){
    mod_names <- paste("Model ", seq(1,length(models)))
  } else {
    mod_names <- names(models)
  }
  aics <- unlist(lapply(models, function(x) x@AIC))
  minAIC <- min(aics)
  delta.aic <- aics - minAIC
  weight <- exp((-1 * delta.aic)/2) / sum(exp((-1 * delta.aic)/2))
  weight <- data.frame(Model = mod_names, AIC = aics, DeltaAIC = delta.aic, AICweight = round(weight, 3))
  rownames(weight) <- NULL
  return(weight)
}

## Overall model selection - would be interesting to see if this is the case for individual species models too!
AICweight(models = list("Elev" = fm4, "Elev Sq" = fm5))

# Look at species-specific random intercepts and slopes
rt <- randomTerms(fm4)
rt

rtmean <- randomTerms(fm4, addMean = TRUE)
rtmean


## Compare Naive to Estimated
## We can see some of the issues with diversity estimates
## without accounting for imperfect detection...
## Some sites have lower observed richness than expected!
## This problem might become more extreme with more species
## and rare species
## We can get richness values from this model
r <- richness(fm4, posterior = T) ## posterior = T gets us uncertainty
est <- apply(r@samples, 1, mean)
low <- apply(r@samples, 1, quantile, 0.025)
high <- apply(r@samples, 1, quantile, 0.975)


naive.rich <- apply(apply(y.perm, c(1,3), function(x) if(all(is.na(x))) return(NA) else return(max(x, na.rm = T))), 1, sum, na.rm = T)
plot(naive.rich, est, xlab="Observed Richness", ylab="Estimated Richness", main="Species Richness", pch = 19, col = rgb(0,0,0,alpha = 0.2))
segments(naive.rich, low, naive.rich, high)
abline(a=0, b=1)
