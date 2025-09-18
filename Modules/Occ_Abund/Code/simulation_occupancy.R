set.seed(123)


#Manually 
N <- 100 # Number of sites
psi <- 0.7 # True occupancy probability
occupancy_status <- rbinom(N, 1, psi) # 1 for present, 0 for absent

J <- 3 # Number of surveys per site
p <- 0.4 # True detection probability
detection_history <- matrix(NA, nrow = N, ncol = J)

for (i in 1:N) {
  if (occupancy_status[i] == 1) {
    detection_history[i, ] <- rbinom(J, 1, p) # Simulate detections for present sites
  } else {
    detection_history[i, ] <- 0 # No detections for absent sites
  }
}

# Example with a site-level covariate for occupancy
site_covariate <- runif(N, 0, 1)
beta_psi <- c(-1, 2) # Intercept and slope for occupancy
logit_psi <- beta_psi[1] + beta_psi[2] * site_covariate
psi_variable <- plogis(logit_psi) # Transform to probability
occupancy_status_cov <- rbinom(N, 1, psi_variable)


#Using unmarked

library(unmarked) 

  # See ??unmarked::simulate
  # Simulation of an occupancy dataset from scratch 

  # First create an unmarkedFrame with the correct design
  
  M <- 300 # number of sites
  J <- 5   # number of occasions
  
  # The values in the y-matrix don't matter as they will be simulated
  # We can supply them as all NAs
  y <- matrix(NA, M, J)
  
  # Site covariate
  x <- rnorm(M)
  
  # Create unmarkedFrame
  umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(x = x))
  
  # Must specify model = occu since unmarkedFrameOccu is also used for occuRN
  # the formula species the specific model structure we want to simulate
  # If we don't specify coefs, unmarked will generate a template you can copy and use
  simulate(umf, model = occu, formula = ~1~x)
  
  # Now set coefs
  # Here we imply a mean occupancy and mean detection of 0.5
  # (corresponding to values of 0 on the inverse link scale) and a positive effect of x
  s <- simulate(umf, model = occu, formula = ~1~x, 
                coefs = list(state = c(0,0.3), det = 0))
  
  head(s[[1]])
  
  occu(~1~x, s[[1]])
  
  # For some models we can also include a random effect
  # add a factor covariate
  umf@siteCovs$x2 <- factor(sample(letters[1:10], M, replace=TRUE))
  
  # The final value in coefs now represents the random effect SD for x2
  s <- simulate(umf, model = occu, formula = ~1~x+(1|x2), 
                coefs = list(state = c(0,0.3, 1), det = 0))
  
  head(s[[1]])
  
  occu(~1~x+(1|x2), s[[1]])
  
  # Here's a more complicated example simulating a gdistsamp dataset
  # using a negative binomial distribution
  M <- 100
  J <- 3
  T <- 2
  y <- matrix(NA, M, J*T)
  umf2 <- unmarkedFrameGDS(y=y, 
                           siteCovs=data.frame(x=rnorm(M)),
                           dist.breaks = c(0, 10, 20, 30), unitsIn='m',
                           numPrimary = T, survey="point")
  
  cf <- list(lambda=c(1, 0.3), phi=0, det=c(log(20), 0), alpha=log(1))
  
  # Note we now also supply another argument mixture="NB" to ... 
  s2 <- simulate(umf2, coefs=cf, lambdaformula=~x, phiformula=~1, pformula=~x,
                 mixture="NB")
  head(s2[[1]])
  
  gdistsamp(~x, ~1, ~x, s2[[1]], mixture="NB")
  