## -------------------------------------------------------------
##
## Script name: Occupancy Model Data Prep
##
## Script purpose: Prepare data for Bayesian MSOM model fit.
##
## Author: Spencer R Keyser
##
## Date Created: 2024-11-27
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
options(scipen = 6, digits = 4)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(stringr)
library(ggplot2)
library(here)
library(abind)
library(lubridate)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Species Data
##
## -------------------------------------------------------------

## Load in the bird data for 2021
load(here("./Modules/Occ_Abund/Data/Thresh_By_Species/2021_99Conf_OccSppList.RData"))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Handling problem species
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## *************************************************************
##
## Section Notes: Analytical decisions
## 1. Remove raptors due to violation of occupancy assumptions (large home ranges)
## 2. Remove Clark's nutcracker for simialr reason to raptors
## 3. Combine acoustically similar species
##  a. All sapsuckers need to be collapsed
##  b. Plumbeous and Cassin's vireo 
## 
##
## *************************************************************

## ADDED removing species and name change for PACFLY
sp.det.list <- sp.det.list[which(!names(sp.det.list) %in% c("Red-tailed Hawk",
                                                            "Osprey",
                                                            "Red-shouldered Hawk",
                                                            "Clark's Nutcracker",
                                                            "American Kestrel"
))]

names(sp.det.list)[which(names(sp.det.list) == "Pacific-slope Flycatcher")] <- "Western Flycatcher"

combine_species <- function(species_list, species_to_combine, new_name, date_cols) {
  # Start with first species
  combined_df <- species_list[[species_to_combine[1]]]
  
  # Sum values across species
  for(sp in species_to_combine[-1]) {
    combined_df[, date_cols] <- combined_df[, date_cols] + species_list[[sp]][, date_cols]
  }
  
  # Binarize: convert all values > 0 to 1
  combined_df <- combined_df %>%
    mutate(across(all_of(date_cols), ~as.integer(. > 0)))
  
  # Remove original species and add combined
  species_list[species_to_combine] <- NULL
  species_list[[new_name]] <- combined_df
  
  return(species_list)
}

## Combine Sapsuckers
sp.det.list <- combine_species(species_list = sp.det.list, 
                               species_to_combine = c("Red-breasted Sapsucker", "Red-naped Sapsucker", "Williamson's Sapsucker"), 
                               new_name = "Sphyrapicus spp.",
                               date_cols = colnames(sp.det.list[[1]])[str_detect(colnames(sp.det.list[[1]]), "Cell_Unit", negate = T)])

## Combine vireos
sp.det.list <- combine_species(species_list = sp.det.list, 
                               species_to_combine = c("Cassin's Vireo", "Plumbeous Vireo"), 
                               new_name = "Vireo spp.",
                               date_cols = colnames(sp.det.list[[1]])[str_detect(colnames(sp.det.list[[1]]), "Cell_Unit", negate = T)])

names(sp.det.list)
str(sp.det.list)

## Cell Unit mapping file
cu.map <- data.frame(ID = 1:length(sp.det.list[[1]]$Cell_Unit), Cell_Unit = sp.det.list[[1]]$Cell_Unit)
cu.map <- cu.map |>
  mutate(Cell_Unit = ifelse(
    stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
    gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
    Cell_Unit
  ))

## Reduce the number of dates for sampling
## Check the duration of the sampling period
## 60 sampling periods
length(colnames(sp.det.list[[1]])) - 1

## Array for the data
## D1 (i) = Site, D2 (j) = Sampling Date, D3 (k) = species
samp.cols <- colnames(sp.det.list[[1]])[str_detect(colnames(sp.det.list[[1]]), "\\d")]
samp.cols <- as.Date(samp.cols, format = "%Y-%m-%d")

## Create different sampling periods
second_samp <- function(DAT, interval, id_col, eff = F, e.var = "Days"){
  tmp.cols <- colnames(DAT)[str_detect(colnames(DAT), "\\d")]
  tmp.splits <- split(tmp.cols, ceiling(seq_along(tmp.cols) / interval))
  tmp.new <- as.data.frame(cbind(DAT[, id_col], matrix(data = NA, nrow = nrow(DAT), ncol = length(tmp.splits), dimnames = list(NULL, paste0("J", seq(1:length(tmp.splits)))))))
  for(i in 1:length(tmp.splits)){
    if(!eff){
      Jsum <- rowSums(DAT[, tmp.splits[[i]]], na.rm = T)
      Jsum <- ifelse(Jsum > 0, 1, 0)
    } 
    if(eff & e.var == "Days"){
      Jsum <- ifelse(DAT[, tmp.splits[[i]]] > 0, 1, 0)
      Jsum <- rowSums(Jsum, na.rm = T)
    }
    if(eff & e.var == "Hrs"){
      Jsum <- rowSums(DAT[, tmp.splits[[i]]], na.rm = T)
    }
    if(eff & e.var == "FirstJDay"){
      mjd <- as.Date(gsub("[[:punct:]]", "-", tmp.splits[[i]]), format = "%Y-%m-%d")
      mjd <- median(lubridate::yday(mjd))
      Jsum <- mjd
    }
    
    tmp.new[,i+1] <- Jsum
  }
  return(tmp.new)
}

## Apply the function
sp.det.list.r <- lapply(sp.det.list, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = F))

samp.cols <- colnames(sp.det.list.r[[1]])[str_detect(colnames(sp.det.list.r[[1]]), "J")]
nsite <- nrow(sp.det.list.r[[1]])
nrep <- ncol(sp.det.list.r[[1]]) - 1
nspec <- length(names(sp.det.list.r))
sp.det <- lapply(sp.det.list.r, function(x) x |> select(all_of(samp.cols)))
y <- array(unlist(lapply(sp.det, as.matrix)), dim = c(nsite, nrep, nspec))
dimnames(y) <- list(cu.map$Cell_Unit, NULL, names(sp.det.list.r))

## No NAs present in the data
table(nsurveys <- apply(y[,,1], 1, function(x) sum(!is.na(x))))

## Species with 0 occurrences
## 9 species
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == -Inf] <- NA

## Naive Occupancy Estimates
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE) / nrow(tmp))

drop.sp <- which(obs.occ == 0)
y <- y[,,-drop.sp]
sp.names <- dimnames(y)
sp.df <- data.frame(Index = 1:length(sp.names[[3]]),
                    Species = sp.names[[3]])
#write.csv(sp.df, file = here::here("Code/Occupancy_Modeling/SpeciesIndex_Filtered_VarThresh.csv"))

## Redefine nspec
nspec <- dim(y)[3]
# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum)) # Compute and print sorted species counts
hist(C, breaks = 30, main = "Naive Species Richness")

## We notice some sites with zero detections for all species
## we should remove those
zero.obs <- which(C == 0)
aru.drop <- names(zero.obs)
cu.map <- cu.map[!cu.map$Cell_Unit %in% aru.drop,]
y <- y[-zero.obs,,]
nrow(y)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Data checking
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Double check the different cutoffs for each species
print(y)
sp.ind <- which(dimnames(y)[[3]] == "Acorn Woodpecker")

test <- y[,,sp.ind]
test <- apply(test, 1, function(x) max(x, na.rm = T))
test <- data.frame(Det = test, Cell_Unit = cu.map$Cell_Unit)
nrow(test[test$Det > 0,])/nrow(test)

## Load in the spatial coordinates
source(here("./Code/Acoustic_Data_Prep/Sierra_functions.R"))
aru_sf <- aru_sf_query(years = c(2021))
aru_sf <- aru_sf |> 
  filter(Cell_Unit %in% test$Cell_Unit) |>
  left_join(test)

mapview::mapview(aru_sf, zcol = "Det")

## -------------------------------------------------------------
##
## Begin Section: Detection Covariates
##
## -------------------------------------------------------------

## Get the number of hours per survey for the detection covariate
eff.dat <- read.csv(here("./Modules/Occ_Abund/Data/Thresh_By_Species/2021_99Conf_OccEffortFileSubset.csv"))
eff.dat$Cell_Unit <- ifelse(stringr::str_detect(eff.dat$Cell_Unit, "^C\\d{4}"), eff.dat$Cell_Unit, gsub("C", "C0", eff.dat$Cell_Unit))
eff.dat <- eff.dat[,-1]
colnames(eff.dat) <- gsub("[[:punct:]]", "_", gsub("X", "", colnames(eff.dat)))
eff.dat <- eff.dat[eff.dat$Cell_Unit %in% dimnames(y)[[1]],]

## Summarize the effort data at the same temporal interval
eff.days <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Days")
colnames(eff.days)[colnames(eff.days) == "V1"] <- "Cell_Unit"
eff.days <- eff.days[,-1]
colnames(eff.days) <- NULL

## Summarize the total number of hours surveyed per sampling unit
eff.hrs <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Hrs")
colnames(eff.hrs)[colnames(eff.hrs) == "V1"] <- "Cell_Unit"
eff.hrs <- eff.hrs[,-1]
colnames(eff.hrs) <- NULL


## Summarize the median Jdate
eff.jday <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "FirstJDay")
colnames(eff.jday)[colnames(eff.jday) == "V1"] <- "Cell_Unit"
eff.jday <- eff.jday[,-1]
colnames(eff.jday) <- NULL

## Check that the biological data matches the effort data
all.equal(eff.dat$Cell_Unit, dimnames(y)[[1]])

## Set y to NA for days without sampling
## Iterate through the matrix and set corresponding values in the 3D array to NA
for (i in 1:nrow(eff.days)) {
  for (j in 1:ncol(eff.days)) {
    if (as.numeric(eff.days[i, j]) == 0) {
      # Set all layers (third dimension) of the array at position (i, j) to NA
      y[i, j, ] <- NA
    }
  }
}

# Print the modified 3D array
print(y)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: "Ragged array" data input for skipping NAs
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ragged <- FALSE
spOcc <- FALSE
unmarked <- TRUE
## *************************************************************
##
## Section Notes:
## We want a way to skip the pesky NAs that comprise ~30% of the
## data and will inevitably slow the model down. The solution
## is to link the data only the the sites with sampling effort
## and skip the NAs. Because the ARUs don't sample at a given site
## this creates a similar scheme for all species and can be handled
## for one species and attributed to all remaining species.
##
## *************************************************************
if(ragged){
  ## One slice of the species matrix
  y_sub <- y[,,1]
  
  ## Find the indices with NAs
  has_data <- which(
    !is.na(y_sub),
    arr.ind = T
  )
  
  ## Make the covariates into long format
  obs_cov_long <- matrix(
    NA,
    nrow(has_data),
    ncol = 3
  )
  
  ## Place the correct data into the new long format
  eff.hrs <- as.matrix(eff.hrs)
  eff.days <- as.matrix(eff.days)
  eff.jday <- as.matrix(eff.jday)
  
  ## Make these in long format
  for(i in 1:nrow(has_data)){
    obs_cov_long[i,1] <- eff.hrs[
      has_data[i,1], # site
      #1,
      has_data[i,2] # rep
    ]
    
    obs_cov_long[i,2] <- eff.days[
      has_data[i,1], # site
      #1,
      has_data[i,2] # rep
    ]
    
    obs_cov_long[i,3] <- eff.jday[
      has_data[i,1], # site
      #1,
      has_data[i,2] # rep
    ]
    
  }
  
  head(obs_cov_long)
  
  ## Scale the detection covariates
  ## Hours ARU was active
  eff.hrs <- obs_cov_long[,1]
  mean.eff.hrs <- mean(eff.hrs)
  sd.eff.hrs <- sd(eff.hrs)
  eff.hrs.scale <- (eff.hrs - mean.eff.hrs) / sd.eff.hrs
  
  ## Days ARU is active
  eff.days <- obs_cov_long[,2]
  mean.eff.days <- mean(eff.days)
  sd.eff.days <- sd(eff.days)
  eff.days.scale <- (eff.days - mean.eff.days) / sd.eff.days
  
  ## DOY
  eff.jday <- obs_cov_long[,3]
  mean.eff.jday <- mean(eff.jday)
  sd.eff.jday <- sd(eff.jday)
  eff.jday.scale <- (eff.jday - mean.eff.jday) / sd.eff.jday
  
  ## Make the response variable in long format
  y_long <- matrix(data = NA, nrow = length(y_sub[!is.na(y_sub)]), ncol = dim(y)[3])
  for(i in 1:dim(y)[3]){
    y.tmp <- y[,,i]
    y_long[,i] <- y.tmp[!is.na(y.tmp)]
  }
  
} #if ragged

## -------------------------------------------------------------
##
## Begin Section: Occupancy Covariates
##
## -------------------------------------------------------------

## Right now lets pull in the ARU meta data
aru_meta <- read.csv(here("Modules/Occ_Abund/Data/ARU_120m_New.csv"))
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)
names(aru_meta)

## Select the variables we are interested in for detection and occupancy
aru_meta <- aru_meta |> select(Cell_Unit, 
                               utme, 
                               utmn,
                               X, 
                               Y,
                               topo_elev, 
                               topo_tpi,
                               tmx_bcm_mn, 
                               ppt_bcm_mn,
                               #fire1yr_cbi_mn, 
                               #fire2_5yr_cbi_mn,
                               fire1_5yr_cbi_mn = fire5yr_cbi_mn, #correct column for 1-5 year fire
                               fire6_10yr_cbi_mn, 
                               fire11_35yr_cbi_mn,
                               #fire1yr_high_prop, 
                               #fire1yr_lowmod_prop,
                               fire1_5yr_high_prop = fire5yr_high_prop, #correct column for 1-5 year fire
                               fire1_5yr_lowmod_prop = fire5yr_lowmod_prop,
                               #fire2_5yr_high_prop, 
                               #fire2_5yr_lowmod_prop,
                               fire6_10yr_high_prop, 
                               fire6_10yr_lowmod_prop,
                               fire11_35yr_high_prop, 
                               fire11_35yr_lowmod_prop,
                               standage_f3_mn, 
                               cpycovr_f3_mn,
                               ch_cfo_mn,
                               ch_cfo_sd,
                               cc_cfo_mn, 
                               cc_cfo_sd
) |> 
  filter(Cell_Unit %in% cu.map$Cell_Unit) |> 
  arrange(match(Cell_Unit, cu.map$Cell_Unit))

aru_cor <- cor(aru_meta |> select(topo_elev,
                                  ppt_bcm_mn,
                                  fire1_5yr_cbi_mn,
                                  fire6_10yr_cbi_mn,
                                  fire11_35yr_cbi_mn,
                                  ch_cfo_mn,
                                  cc_cfo_mn,
                                  Y,
                                  X))

## CC and CFO are highly correlated
corrplot::corrplot(aru_cor,
                   method = "number")

## Create orthogonal predictor by lm
plot(aru_meta$cc_cfo_mn, aru_meta$ch_cfo_mn)
ch_res <- lm(ch_cfo_mn ~ cc_cfo_mn, aru_meta)$residual

aru_meta <- aru_meta |> mutate(ch_res = ch_res)

coords <- aru_meta |> 
  select(Cell_Unit, Y, X, ch_res) |> 
  sf::st_as_sf(coords = c("X", "Y"), crs = 4326) |> 
  sf::st_transform(crs = 3310) |> 
  sf::st_coordinates() |> 
  as.data.frame()

## Scale the preds of interest
{
  ## Spatial data
  utme <- as.vector(scale(aru_meta$utme))
  utmn <- as.vector(scale(aru_meta$utmn))
  Lat <- as.vector(scale(aru_meta$Y))
  Long <- as.vector(scale(aru_meta$X))
  
  ## Elevation and climate data
  ele <- as.vector(scale(aru_meta$topo_elev))
  ppt <- as.vector(scale(aru_meta$ppt_bcm_mn))
  tmx <- as.vector(scale(aru_meta$tmx_bcm_mn))
  
  ## Original CBI mean data
  #cbi1 <- as.vector(scale(aru_meta$fire1yr_cbi_mn))
  cbi1_5 <- as.vector(scale(aru_meta$fire1_5yr_cbi_mn))
  #cbi2_5 <- as.vector(scale(aru_meta$fire2_5yr_cbi_mn))
  cbi6_10 <- as.vector(scale(aru_meta$fire6_10yr_cbi_mn))
  cbi11_35 <- as.vector(scale(aru_meta$fire11_35yr_cbi_mn))
  
  ## Adding in the proportional fire data
  #hsf_prop1 <- as.vector(scale(aru_meta$fire1yr_high_prop))
  hsf_prop1_5 <- as.vector(scale(aru_meta$fire1_5yr_high_prop))
  #hsf_prop2_5 <- as.vector(scale(aru_meta$fire2_5yr_high_prop))
  hsf_prop6_10 <- as.vector(scale(aru_meta$fire6_10yr_high_prop))
  hsf_prop11_35 <- as.vector(scale(aru_meta$fire11_35yr_high_prop))
  #lmsf_prop1 <- as.vector(scale(aru_meta$fire1yr_lowmod_prop))
  lmsf_prop1_5 <- as.vector(scale(aru_meta$fire1_5yr_lowmod_prop))
  #lmsf_prop2_5 <- as.vector(scale(aru_meta$fire2_5yr_lowmod_prop))
  lmsf_prop6_10 <- as.vector(scale(aru_meta$fire6_10yr_lowmod_prop))
  lmsf_prop11_35 <- as.vector(scale(aru_meta$fire11_35yr_lowmod_prop))
  
  ## Forest characteristics
  stage <- as.vector(scale(aru_meta$standage_f3_mn)) 
  cc_f3 <- as.vector(scale(aru_meta$cpycovr_f3_mn))
  cc_cfo <- as.vector(scale(aru_meta$cc_cfo_mn))
  cc_cfo_sd <- as.vector(scale(aru_meta$cc_cfo_sd))
  ch_cfo <- as.vector(scale(aru_meta$ch_cfo_mn))
  ch_res <- as.vector(scale(aru_meta$ch_res))
  ch_cfo_sd <- as.vector(scale(aru_meta$ch_cfo_sd))
  }

## -------------------------------------------------------------
##
## Begin Section: Prep Data for JAGS
##
## -------------------------------------------------------------
if(ragged){
  dimnames(y) <- NULL
  
  ## Win data new
  win.data.rag <- list(y = y_long,
                       nsite = dim(y)[1],
                       N = nrow(y_long),
                       nspec = ncol(y_long),
                       site_id = has_data[,1],
                       eff.hrs = eff.hrs.scale,
                       eff.jday = eff.jday.scale,
                       #bmass = bm.scale,
                       #beak.pc1 = pc1,
                       #beak.pc2 = pc2,
                       utmn = utmn,
                       lat = Lat,
                       ele = ele,
                       ppt = ppt,
                       tmx = tmx,
                       #cbi1 = cbi1,
                       cbi1_5 = cbi1_5,
                       #cbi2_5 = cbi2_5,
                       cbi6_10 = cbi6_10,
                       cbi11_35 = cbi11_35,
                       stage = stage,
                       cc_f3 = cc_f3,
                       cc_cfo = cc_cfo,
                       ch_cfo = ch_cfo,
                       ch_res = ch_res
  )
  str(win.data.rag)
}
#                  nrep = dim(y)[2],
#                  nspec = dim(y)[3],
#                  eff.days = eff.days,
#                  eff.hrs = eff.hrs,
#                  utmn = utmn,
#                  ele = ele,
#                  ppt = ppt,
#                  tmx = tmx,
#                  cbi1 = cbi1,
#                  cbi2_5 = cbi2_5,
#                  cbi6_10 = cbi6_10,
#                  cbi11_35 = cbi11_35,
#                  stage = stage,
#                  cc = cc
# )

if(spOcc){
  
  eff.hrs.m <- as.matrix(eff.hrs)
  eff.hrs.m[eff.hrs.m == 0] <- NA
  eff.jday.m <- as.matrix(eff.jday)
  eff.jday.m[is.na(eff.hrs.m)] <- NA
  
  y.new <- y
  dimnames(y.new) <- list(
    sites = aru_meta$Cell_Unit,
    reps = paste0("Rep_", seq(1:5)),
    species = dimnames(y.new)[[3]]
  )
  y.new <- aperm(y.new, c(3,1,2))
  str(y.new)
  y.new[1:10, 1:10, 1]
  
  spOcc.data <- list(
    y = y.new,
    det.covs = list(eff.hrs = eff.hrs.m,
                    eff.jday = eff.jday.m),
    occ.covs = data.frame(ele = ele,
                          ppt = ppt,
                          tmx = tmx,
                          cbi1 = cbi1,
                          cbi1_5 = cbi1_5,
                          cbi2_5 = cbi2_5,
                          cbi6_10 = cbi6_10,
                          cbi11_35 = cbi11_35,
                          stage = stage,
                          cc_f3 = cc_f3,
                          cc_cfo = cc_cfo,
                          ch_cfo = ch_cfo,
                          ch_res = ch_res),
    coords = coords)
}

if(unmarked){
  
  ## Make incomplete sampling events NA to avoid estimating them
  eff.hrs.m <- as.matrix(eff.hrs)
  eff.hrs.m[eff.hrs.m == 0] <- NA
  eff.jday.m <- as.matrix(eff.jday)
  eff.jday.m[is.na(eff.hrs.m)] <- NA
  
  ## Det/Non-det data
  y.new <- y
  dimnames(y.new) <- list(
    sites = aru_meta$Cell_Unit,
    reps = paste0("Rep_", seq(1:5)),
    species = dimnames(y.new)[[3]]
  )
  y.new <- aperm(y.new, c(3,1,2))
  str(y.new)
  y.new[1:10, 1:10, 1]
  
  unmarked.data <- list(
    y = y.new,
    obsCovs = list(eff.hrs = eff.hrs.m,
                    eff.jday = eff.jday.m),
    siteCovs = aru_meta)
}

## -------------------------------------------------------------
##
## End Section:
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Cleaning file and saving .RDATA
##
## -------------------------------------------------------------
if(ragged){
  to_keep <- c("y", "y_long", "win.data.rag")
  to_remove <- setdiff(ls(), to_keep)
  
  rm(list = to_remove)
  rm(to_remove)
  
  ## Save the RDATA
  save.image(file = here("./Data/JAGS_Data/MSOM_Ragged_2021_SpeciesThresh_975minMaxPrex_NewVars.RData"))
}

if(spOcc){
  saveRDS(spOcc.data, here("./Data/SpOccupancy_Data/SpOccData_SpeciesThresh_975minMaxPrex_NewVars.rds"))
}

if(unmarked){
  saveRDS(unmarked.data, here("./Modules/Occ_Abund/Data/Occ_Data/UnmarkedData_SpeciesThresh_975minMaxPrex_NewVars.rds"))
}
