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
## Notes: Workshop material
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
occ.dat$y[1,,]

## Reorder the array so it's easier to use
occ.dat$y <- aperm(occ.dat$y, c(2,3,1))

## To reduce the amount of data let's subset to 500 random locations with a minimum of 1km between points
siteCovs <- occ.dat$siteCovs
siteCovs$Cell <- as.factor(gsub("_U\\d{1}", "", siteCovs$Cell_Unit))
coords <- siteCovs |> 
  sf::st_as_sf(coords = c("X", "Y"), crs = 4326) 

## Subsample
spsub <- spatialEco::subsample.distance(coords, size = 500, d = 1000)
siteCovs <- spsub
subARU <- as.numeric(rownames(siteCovs))
## Break the data apart for demonstration purposes
y <- occ.dat$y 
#siteCovs <- occ.dat$siteCovs
obsCovs <- occ.dat$obsCovs

## Subsampled sites
y <- y[subARU,,]
#siteCovs <- siteCovs[rs,]
obsCovs <- lapply(obsCovs, FUN = function(x) x[subARU,])

all(which(is.na(y[,,1])) == which(is.na(obsCovs$eff.hrs)))

# Take one species since we are doing a single species model
## Let's focus on Hermit warbler
sp.ind <- which(dimnames(y)[[3]] == "Mountain Quail")
y <- y[,,sp.ind]

## Remove the spatial features of the siteCovs
siteCovs <- sf::st_drop_geometry(siteCovs)
all(rownames(y) == siteCovs$Cell_Unit)

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
plogis(0.593) # ~64% probability the species occupies a site
plogis(1.44) # ~81% probability we detect the species if it occupies a site

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
eff.hrs.max <- max(obsCovs$eff.hrs, na.rm = T)
eff.jday.max <- max(obsCovs$eff.jday, na.rm = T)
nd <- data.frame(eff.hrs = eff.hrs.max, eff.jday = eff.jday.max)
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
## Subsection: Even more realistic model
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fm4 <- occu(~ scale(eff.hrs) + 
              scale(eff.jday) + 
              scale(I(eff.jday^2)) + 
              scale(cc_cfo_mn) + 
              scale(topo_tpi) + 
              (1|Cell) 
            ~ scale(utmn) + 
              scale(topo_elev) + 
              scale(I(topo_elev^2)) + 
              scale(ch_cfo_mn), 
            data = unmk)

fm4

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Model Selection and Fit
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Unmarked fitList()
fms <- fitList('psi(.)p(.)' = fm1, 
               'psi(.)p(Hours+Date)' = fm2,
               'psi(elevation)p(Hours+Date)' = fm3,
               'psi(ele+ele^2+cc)p(Hours+Date+Date^2+cc)' = fm4)
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

## Following MacKenzie and Bailey GoF (2004) the p-value
## is the number of simulations in which our model Chi-square is
## greater than the simulated data (> 0.05 is adequate model fit)
set.seed(182)
pb <- parboot(fm4, statistic = chisq, nsim = 999, parallel = FALSE)

## We have some strong unmodelled detection heterogeneity
## What could we do to improve the fit here? Maybe go back to the acoustics for some information?
mbgof <- mb.gof.test(fm4, nsim = 10, print.table = T)

## Adjust by overdispersion
c_hat <- 2.16
aic_table <- modSel(fms)
qaic_table <- aic_table
qaic_table@Full$AIC <- qaic_table@Full$AIC / c_hat  # Convert AIC to QAIC
qaic_table@Full$delta <- qaic_table@Full$AIC - min(qaic_table@Full$AIC)
print(qaic_table)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: False-positive models
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Code below is for a false-positive flavor of occupancy model
## False-positive models are likely critical for occupancy modelling
## with bioacoustics data derived from ML classifiers like (BirdNET)
## and detection thresholds. Given the rigidity of assumptions
## in standard occupancy models that there are zero FPs this
## extension is particularly relevant for all participants.
## However, these models are particularly tricly in lieu of
## extensive validation efforts and suffer from parameter
## identifiability issues. That is to say that it's impossible
## for us to distinguish between scenarios in which p11 and p10
## such that the model is multi-modal (Link and Royle, 2006). 
## We can use some assumptions
## and numerical "tricks" to get the model to avoid sign switching
## but this is dependent on the species, acoustic quality, classifier performance, 
## amount of validated data, and model structure.
## The models below correspond to a "Type 2" situation where we don't have 
## explicitly validated data with the acoustic detections. "Type 3" can 
## improve parameter estimation and multi-modality by integrating
## "known" detections (e.g., expert listeners detect the species of interest).
## This data can then be integrated into the model as a multi-state model
## to help resolve the likelihood of an observation belonging to FP class vs
## TP class sensu Miller et al., (2013) and Clements et al, (2017).

## Some relevant sources
## Link and Royle (2006)
## Miller et al., (2013)
## Clements et al., (2017)
## Rhinehart et al., (2023)

## Adding species richness as an indicator for potentially "noisy" environments
## Can you think of other indicators to provide the model to help delineate
## FPs? What might they be?
rich <- occ.dat$y
rich <- apply(rich, c(1,2), sum)
rich <- rich[subARU,]
## We assume all data can be impacted by FP
type <- c(0,5,0)
obsCovsFP <- list(eff.hrs = obsCovs$eff.hrs,
                eff.jday = obsCovs$eff.jday,
                rich = rich,
                METH = y)
unmk.fp <- unmarkedFrameOccuFP(y = y, siteCovs = siteCovs, obsCovs = obsCovsFP, type = type)
class(unmk.fp)
summary(unmk.fp)

## We add starting values to try to avoid incorrect modality
## We state that p11 >> p10 (i.e., FP should be quite small)
largerp11 <- qlogis(c(0.5, 0.7, 0.1))
largerp11 <- qlogis(c(0.5, 0.5, 0.5, 0.5, 0.5, #occupancy (int + covs) 
                      0.7, 0.5, 0.5, 0.5, 0.5, #detection (int + covs)
                      0.1, 0.5 #fp (int + covs)
                      ))

fm5 <- occuFP(detformula = ~1, 
              stateformula = ~1,
              FPformula = ~1,
              data = unmk.fp,
              starts = largerp11)

fm6 <- occuFP(detformula = ~scale(eff.hrs) + scale(eff.jday) + scale(I(eff.jday^2)) + scale(cc_cfo_mn), 
              stateformula = ~scale(utmn) + scale(topo_elev) + scale(I(topo_elev^2)) + scale(ch_cfo_mn),
              FPformula = ~scale(rich),
              data = unmk.fp,
              starts = largerp11)

fm6

## False-positive rate confidence
unmarked::confint(fm6, type = "fp")
plogis(-7.2)
plogis(-3.5)

## Expected number of potential FPs at the average site
## ~10 observations
sum(!is.na(y)) * plogis(-5.34)

## AIC (False-positive model is not inherently better than best fitting non-FP model)
fm6@AIC > fm4@AIC
(fm4@AIC - fm6@AIC) > 2

fm4; plogis(-0.254)
fm6; plogis(-0.34)

## -------------------------------------------------------------
##
## Begin Section: Visualizations
##
## -------------------------------------------------------------

## Quick plotting via unmarked
plotEffects(fm4, "state", "topo_elev")
plotEffects(fm4, "det", "eff.hrs")
plotEffects(fm4, "det", "eff.jday")
plotEffects(fm6, "fp", "rich")
plotEffects(fm6, "state", "ch_cfo_mn")


plot.df <- plotEffectsData(fm4, "state", "topo_elev")

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

ele.plot

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
