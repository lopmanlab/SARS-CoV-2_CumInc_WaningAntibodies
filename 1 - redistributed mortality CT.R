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
# Load the GA mortality data
#------------------------------------------------------------------------------#

ga <- read.csv('./Data/GADeaths.csv')

#------------------------------------------------------------------------------#
# Reformat and clean the data
#------------------------------------------------------------------------------#

# Date variables
ga$onset <- as.Date(ga$OnsetDate, format="%Y/%m/%d")
ga$death <- as.Date(ga$DateDeath, format="%Y/%m/%d")

# Some dates were very old, so remove them
ga1 <- ga[which(ga$onset > as.Date("2020-01-01") & ga$death > as.Date("2020-01-01")),]

# Calculate the time from symptom onset to death
ga1$time_onset_death <- as.numeric(difftime(ga1$death, ga1$onset, units = c("days")))

# Remove negative values (i.e., deaths happened before symptom onset)
# Also, for Weibull and gamma, t should not be zero, so remove zeros as well
ga2 <- ga1[which(ga1$time_onset_death>0), ] 

#------------------------------------------------------------------------------#
# Fit gamma, Weibull, lognormal distributions to time from onset to death
#------------------------------------------------------------------------------#

library(MASS)
gamma.par <- fitdistr(ga2$time_onset_death, densfun = 'gamma',   start=list('shape'=0.2,'rate'=0.1))[[1]]
weib.par  <- fitdistr(ga2$time_onset_death, densfun = 'weibull', start=list('shape'=0.2,'scale'=0.1))[[1]]
lnorm.par <- fitdistr(ga2$time_onset_death, densfun = 'lognormal')[[1]]

# Check the fit
g <- dgamma(  seq(0, max(ga2$time_onset_death), by=0.01), shape=gamma.par[1], rate=gamma.par[2])
w <- dweibull(seq(0, max(ga2$time_onset_death), by=0.01), shape=weib.par[1],  scale=weib.par[2])
ln <- dlnorm( seq(0, max(ga2$time_onset_death), by=0.01), meanlog=lnorm.par[1],  sdlog=lnorm.par[2])
par(mfrow=c(1,1))
hist(ga2$time_onset_death, breaks=50, col="grey") # Observed data
par(new=T) 
plot(y=g, x=seq(0, max(ga2$time_onset_death), by=0.01), # Fitted gamma distribution
     xlim=c(0,max(ga2$time_onset_death)), axes=F, ann=F, type='l', lwd=3, col='red')
par(new=T)
plot(y=w, x=seq(0, max(ga2$time_onset_death), by=0.01), # Fitted Weibull distribution
     xlim=c(0,max(ga2$time_onset_death)), axes=F, ann=F, type='l', lwd=3, col='blue')
par(new=T)
plot(y=ln, x=seq(0, max(ga2$time_onset_death), by=0.01), # Fitted lognormal distribution
     xlim=c(0,max(ga2$time_onset_death)), axes=F, ann=F, type='l', lwd=3, col='orange')

# Goodness of fit test
library(actuar)
library(fitdistrplus)
gamma.par1 <- fitdist(ga2$time_onset_death, 'gamma')
weib.par1  <- fitdist(ga2$time_onset_death, 'weibull')
lnorm.par1 <- fitdist(ga2$time_onset_death, 'lnorm')
gofstat(list(gamma.par1, weib.par1, lnorm.par1), fitnames = c("gamma", "weibull", "lognormal"))
# Log normal is selected as best by all tests

#------------------------------------------------------------------------------#
# Load the CT mortality data and reformat it
#------------------------------------------------------------------------------#

# Select the site
select.site <- "CT"

# Select the type of death counts: 
select.death.outcome <- "Deaths.7dAvg" 

# Load the mortality data and CDC serosurvey data, and reformat them
source("0 - data formatting CT.R")

# Add 200 days at the beggining of the mortality data to account for the time 
# lag between symptom onset and death
d <- 200
num_date <- nrow(d1) + d
d0 <- as.data.frame(matrix(NA, d, ncol(d1)))
names(d0) <- names(d1)
d0$date <- seq(d1$date[1]-d, d1$date[1]-1, "day") # Add d days 
d1 <- rbind(d0, d1)

# Load the reported number of deaths
rep_counts <- d1[,select.death.outcome] 

#------------------------------------------------------------------------------#
# Estimate the timing of symptom onset for fatal cases using the GA data
#------------------------------------------------------------------------------#

# CDF for the log-normal distribution fitted to the GA mortality data
Ft <- plnorm(1:num_date, meanlog=lnorm.par[1],  sdlog=lnorm.par[2]) # CDF

# Empty matrix to store results
p2 <- matrix(0, nrow=num_date, ncol=num_date) 

# Estimate the timing of symptom onset for these cases died on day a
for (a in (d+2):(num_date)) { # "d" is just extra days added at the beginning of the mortality data so that we can go back to extra days before the pandemic (e.g., January 2020)
 for (k in 2:a) { # k is from *2* to a, becasuse F[k-1] appears below. k-1 cannot be zero.
  p2[a-k+1, a] <- rep_counts[a] * (Ft[k] - Ft[k-1]) # Prob that cases who died on day a were infected on day a-k
 }
}

# Save the result
timing_onset.ga <- rowSums(p2)

#------------------------------------------------------------------------------#
# Estimate the timing of symptom onset for fatal cases using the Wuhan, China 
# data (O'Driscoll, et al. Nature 2020 paper)
#------------------------------------------------------------------------------#

set.mean <- 20 # From O'Driscoll, et al. Nature 2020 
set.sd   <- 10 # From O'Driscoll, et al. Nature 2020
rate_odriscoll  <- set.mean / set.sd^2
shape_odriscoll <- rate_odriscoll * set.mean

# CDF for the gamma distribution from O'Driscoll, et al. 
Ft.ODriscoll <- pgamma(1:num_date, shape=shape_odriscoll, rate=rate_odriscoll)

# Empty matrix to store results
p3 <- matrix(0, nrow=num_date, ncol=num_date) 

# Estimate the timing of symptom onset for these cases died on day a
for (a in (d+2):(num_date)) { # "d" is just extra days added at the beginning of the mortality data so that we can go back to extra days before the pandemic (e.g., January 2020)
 for (k in 2:a) { # k is from *2* to a, becasuse F[k-1] appears below. k-1 cannot be zero.
  p3[a-k+1, a] <- rep_counts[a] * (Ft.ODriscoll[k] - Ft.ODriscoll[k-1]) # Prob that cases who died on day a were infected on day a-k
 }
}

timing_onset.ODriscoll <- rowSums(p3)

