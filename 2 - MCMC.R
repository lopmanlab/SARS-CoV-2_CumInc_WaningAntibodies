################################################################################
# 
#   PROJECT: Estimating the cumulative incidence of SARS-CoV-2 infection and 
#            the infection fatality ratio in light of waning antibodies
#
#      DATE: October 2020
#  CODED BY: Kayoko Shioda
#    E-MAIL: kayoko.shioda@emory.edu; kayoko.shioda@aya.yale.edu
# 
################################################################################

#------------------------------------------------------------------------------#
# Set up
#------------------------------------------------------------------------------#

# Which site are you analyzing?
select.site <- "CT" # "NYC" or "CT"

# Which data do you want to use to inform the delay between symptom onset and death?
time_onset_death <- "GA" # "GA" (GA data) or "ODriscoll" (Wuhan data)

# Do you want to save MCMC outputs?
save.results <- "yes" # "yes" or "no"

# If yes, specify a folder path where you want to save your RData
res.save.folder <- "./Results//"

#------------------------------------------------------------------------------#
# Create a function to estimate time between seroconversion and seroreversion
#------------------------------------------------------------------------------#

FQ_dist <- function (t_lag_all, shape1, shape2, scale1, scale2) { 
 pr <- {}
 for (i in 1:length(t_lag_all)) {
  pr[i] <- pweibull(t_lag_all[i], shape1, scale1) - integrate(f = function (y) {
   dweibull(y, shape1, scale1) * (pweibull(t_lag_all[i]-y, shape2, scale2))}, 
   lower=0, upper=t_lag_all[i])$value 
 }
 return(pr)
}

#------------------------------------------------------------------------------#
# Load and reformat the COVID-19 mortality data and CDC serosurvey data
#------------------------------------------------------------------------------#

# Select the type of death counts
if (select.site=="NYC") {
 # "Combined", "Probable", or "Confirmed"
 # "Combined" was used in the main analysis in the paper for NYC 
 select.death.outcome <- "Combined" 
} else if (select.site=="CT") {
 # "Daily.deaths" or "Deaths.7dAvg"
 # "Deaths.7dAvg" was used in the main analysis in the paper for CT
 select.death.outcome <- "Deaths.7dAvg" 
}

# Load the COVID-19 mortality data and CDC serosurvey data, and reformat them
if (select.site=="NYC") {
 source("0 - data formatting NYC.R")
} else if (select.site=="CT") {
 source("0 - data formatting CT.R")
}

# Estimate the timing of symptom onset for fatal cases
if (select.site=="NYC") {
 source("1 - redistributed mortality NYC.R")
} else if (select.site=="CT") {
 source("1 - redistributed mortality CT.R")
}
if (time_onset_death == "GA") {
 timing_onset <- timing_onset.ga
} else if (time_onset_death == "ODriscoll") {
 timing_onset <- timing_onset.ODriscoll
}

# Load the reported number of COVID-19 associated deaths
rep_counts <- d1[,select.death.outcome] 

# Specify the population size in each site
if (select.site=="NYC") {
 pop <- 8.336817*1000000 # From https://www.census.gov/quickfacts/fact/table/newyorkcitynewyork,bronxcountybronxboroughnewyork,kingscountybrooklynboroughnewyork,newyorkcountymanhattanboroughnewyork,queenscountyqueensboroughnewyork,richmondcountystatenislandboroughnewyork/PST045219
} else if (select.site=="CT") {
 pop <- 3.7*1000000
}

#------------------------------------------------------------------------------#
# Create a log-likelihood function for MCMC
#------------------------------------------------------------------------------#

weibull_hazard_LL <- function(set.rho, set.mean, set.sd) {
 
 #-----*-----*-----*-----*-----*-----*-----*-----*-----#
 # 1. True number of infections
 #-----*-----*-----*-----*-----*-----*-----*-----*-----#
 
 # Calculate the number of true infections by adjusting the reported number of 
 # deaths for time lag between symptom onset and death and IFR (set.rho)
 true_inf <- rep(NA, num_date) 
 for (t1 in 1:(num_date)) {
  true_inf[t1] <- min(round(timing_onset[t1] / set.rho), pop) # True infections should not exceed the whole population
 }

 #-----*-----*-----*-----*-----*-----*-----*-----*-----#
 # 2. Timeline of seroconversion & seroreversion
 #-----*-----*-----*-----*-----*-----*-----*-----*-----#
 
 # Specify the shape and scale parameters for a Weibull distribution for time 
 # from symptom onset to seroconversion
 shape_conv <- 2.092141 # From Iyer, et al. medRxiv 2020
 scale_conv <- 13.00029 # From Iyer, et al. medRxiv 2020
 
 # Calculate the shape and scale parameters for a Weibull distribution for time 
 # from seroconversion to seroreversion
 mean_serorev <- set.mean # From MCMC sample
 sd_serorev <- set.sd     # From MCMC sample
 sum.r <- (mean_serorev^2 + sd_serorev^2)/mean_serorev^2
 fn <- function (x, sum) {lgamma(1+2/x) - 2*lgamma(1+1/x) - log(sum)}
 shape_rev <- uniroot(fn, sum=sum.r, c(0.01,100), extendInt = "yes", maxiter = 10000)$root
 scale_rev <- mean_serorev/(gamma(1+1/shape_rev))
 
 #-----*-----*-----*-----*-----*-----*-----*-----*-----#
 # 3. Seropositive individuals on Day t
 #-----*-----*-----*-----*-----*-----*-----*-----*-----#
 
 # Calculate the probability that individuals are seropositive on each day
 FQ_cdf <- FQ_dist(1:num_date, 
                   shape1=shape_conv, scale1=scale_conv, 
                   shape2=shape_rev, scale2=scale_rev) 
 
 # Calculate the number of seropositive individuals on each day
 seropos_cases <- rep(NA, num_date)
 seropos_cases[1] <- 0 # No individuals have seroconverted on Day 1
 for (t2 in 2:num_date) { # On Day 2, 3, ..., last day,
  seropos_cases_day_i <- c() 
  for (k in 1:(t2-1)) {
   seropos_cases_day_i[k] <- true_inf[k] * FQ_cdf[t2-k]
  }
  seropos_cases[t2] <- sum(seropos_cases_day_i)
  seropos_cases[t2] <- min(seropos_cases[t2], pop) # Seropositive cases should not exceed the whole population
 }
 
 #-----*-----*-----*-----*-----*-----*-----*-----*-----#
 # 4. Log likelihood
 #-----*-----*-----*-----*-----*-----*-----*-----*-----#
 
 log_likl <- c()
 for (round in 1:nrow(cdc_seroprev)) { # For each round of the CDC serosurvey...
  # Calculate the log-likelihood, comparing the simulated seroprevalence and 
  # observed seroprevalence reported by CDC
  log_likl[round] <- dbinom(x = cdc_seroprev$n_pos[round], # Number of seropositives (CDC data)
                            prob = seropos_cases[(d1$date)==cdc_seroprev$mid[round]]/pop, # Simulated seroprevalence
                            size = cdc_seroprev$n_sample[round], # Number of samples collected (CDC data)
                            log=T)
 }
 
 log_LL <- sum(log_likl) # Sum of log-likelihood for all rounds of the CDC serosurvey
 
 return(list(log_LL, seropos_cases))
}

#------------------------------------------------------------------------------#
# Use MCMC to optimize parameters
#------------------------------------------------------------------------------#

# Give starting values for each parameter (arbitrary choice)
rho_0  <- 0.02   # IFR
mean_0 <- 80     # Average time from seroconversion to seroreversion (days)
sd_0   <- 50     # Standard deviation of the time from seroconversion to seroreversion (days)
old_rho  <- rho_0
old_mean <- mean_0
old_sd   <- sd_0

# Number of iterations to run the MCMC sampler
burn_in_steps <- 20000 
run_steps     <- 100000 
total_steps   <- burn_in_steps + run_steps

# Get a starting log_likelihood value
old_log_LL <- weibull_hazard_LL(old_rho, old_mean, old_sd)[[1]] # [[1]] is log_LL. [[2]] is seropos_cases

# Make an empty vector to store sampled parameters in
sampled_pars  <- matrix(0, total_steps, 2) # 3 columns for rho (IFR), mean
sampled_logLL <- matrix(0, total_steps, 2) # 3 columns for rho (IFR), mean
dt.seropos1   <- matrix(0, total_steps, num_date) # Estimated number of seropositives on each day
dt.seropos2   <- matrix(0, total_steps, num_date) # Estimated number of seropositives on each day
dt.seropos3   <- matrix(0, total_steps, num_date) # Estimated number of seropositives on each day
accept_rho  <- rep(0, total_steps) # 1 if a new candidate value is accepted. 0 if rejected
accept_mean <- rep(0, total_steps) # 1 if a new candidate value is accepted. 0 if rejected
accept_sd   <- rep(0, total_steps) # 1 if a new candidate value is accepted. 0 if rejected

# Implement MCMC
for (i in 1:total_steps) {
 
 ###############################
 ########## rho (IFR) ##########
 ###############################
 
 # 1. Propose a new value for rho (IFR)
 current_rho  <- -999
 while (current_rho < sum(timing_onset)/pop | current_rho > 1) { 
  current_rho <- old_rho + 0.001*rnorm(1,0,1)
 }
 
 # 2. Proposal ratio (0 because the proposal distribution was symmetric)
 log_proposal_ratio  <- 0.0
 
 # 3. Get the log-likelihood for the newly sampled parameter
 current_log_LL <- weibull_hazard_LL(current_rho, old_mean, old_sd)[[1]]
 ts_seropos     <- weibull_hazard_LL(current_rho, old_mean, old_sd)[[2]] # This stores the estimated number of seropositives on each day
 
 # 4. Log-likelihood ratio for the last and current sampled parameters
 log_likelihood_ratio <- current_log_LL - old_log_LL
 
 # 5. Get the acceptance probability
 acceptance_prob <- exp(log_proposal_ratio + log_likelihood_ratio)
 
 # 6. Draw a random nunber from (0, 1] and if it's less than the acceptance
 # probability, accept the current parameter. Otherwise, stay on the last one.
 if (runif(1) < acceptance_prob) { 
  old_rho    <- current_rho
  old_log_LL <- current_log_LL
  accept_rho[i] <- 1 # 1 if accepted. 0 if rejected.
 } 
 
 # 7. Save the current parameter to the trace
 sampled_pars[i,1]  <- old_rho
 sampled_logLL[i,1] <- old_log_LL
 dt.seropos1[i,]    <- ts_seropos
 
 ##########################
 ########## Mean ##########
 ##########################
  
 # 1. Propose a new value for the mean of the duration of seroreversion
 current_mean <- -999
 # while (current_mean < 15 | current_mean > 300) {
 current_mean <- old_mean + 1*rnorm(1,0,1)
 # }
 
 # 2. Proposal ratio (0 because the proposal distribution was symmetric)
 log_proposal_ratio  <- 0.0
 
 # # Prior for mean
 # log_prior_new <- dexp(current_mean, rate=0.0001, log=TRUE) # non-informative prior
 # log_prior_old <- dexp(old_mean, rate=0.0001, log=TRUE)
 # log_prior_ratio <- log_prior_new - log_prior_old
 
 # 3. Get the log-likelihood for the newly sampled parameter
 current_log_LL <- weibull_hazard_LL(old_rho, current_mean, old_sd)[[1]]
 ts_seropos     <- weibull_hazard_LL(old_rho, current_mean, old_sd)[[2]] # This stores the estimated number of seropositives on each day
 
 # 4. Log-likelihood ratio for the last and current sampled parameters
 log_likelihood_ratio <- current_log_LL - old_log_LL
 
 # 5. Get the acceptance probability
 acceptance_prob <- exp(log_proposal_ratio + log_likelihood_ratio)
 
 # 6. Draw a random nunber from (0, 1] and if it's less than the acceptance
 # probability, accept the current parameter. Otherwise, stay on the last one.
 if (runif(1) < acceptance_prob) { 
  old_mean   <- current_mean
  old_log_LL <- current_log_LL
  accept_mean[i] <- 1 # 1 if accepted. 0 if rejected.
 } 
 
 # 7. Save the current parameter to the trace
 sampled_pars[i,2]  <- old_mean
 sampled_logLL[i,2] <- old_log_LL
 dt.seropos2[i,]    <- ts_seropos
 
 print(round(i/total_steps, 2))
 
}


# Save MCMC outputs
dt.seropos <- rbind(dt.seropos1, dt.seropos2) 
if (save.results == "yes") {
 save(sampled_pars,  file=paste0(res.save.folder,"sampled_pars.RData"))
 save(sampled_logLL, file=paste0(res.save.folder,"sampled_logLL.RData"))
 save(dt.seropos,    file=paste0(res.save.folder,"dt.seropos.RData"))
}