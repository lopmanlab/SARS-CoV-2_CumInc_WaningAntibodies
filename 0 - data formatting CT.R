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
# Import and reformat the mortality data
#------------------------------------------------------------------------------#

# Downloaded from here https://data.ct.gov/Health-and-Human-Services/COVID-19-Tests-Cases-Hospitalizations-and-Deaths-S/rf3k-f8fg
# at 9:43PM (ET) on Friday, October 6, 2020
d <- read.csv('./Data/CT_COVID-19_Tests__Cases__Hospitalizations__and_Deaths__Statewide_.csv')

# Create a copy of the original dataset
d1 <- d

#------------------------------------------------------------------------------#
# Date
#------------------------------------------------------------------------------#

# Reformat a date variable
d1$date <- as.Date(d1$Date, format="%m/%d/%Y")

# Order by date
d1 <- d1[order(d1$date),]

# Weekends are missing. Add them so that the date has no missings.
date <- seq(d1$date[1], d1$date[nrow(d1)], by="days")
seqnum <- 1:length(date)
alldays <- data.frame(seqnum, date)
d1 <- merge(d1, alldays, by="date", all.y = T)

# On these missing dates (weekends), correct Total.deaths (currently NA)
d1$Total.deaths1 <- d1$Total.deaths
for (t in 2:nrow(d1)) {
 if (is.na(d1$Total.deaths1[t])) {
  d1$Total.deaths1[t] <- d1$Total.deaths1[t-1]
 }
}

#------------------------------------------------------------------------------#
# Reformat the death counts
#------------------------------------------------------------------------------#

# Calculate the daily, new deaths from the cumulative counts (Total.deaths1)
d1$Daily.deaths[1] <- d1$Total.deaths1[1]
for (t in 2:nrow(d1)) {
 d1$Daily.deaths[t] <- d1$Total.deaths1[t] - d1$Total.deaths1[t-1]
}

# Calculate the 7-day rolling average death counts
d1$Deaths.7dAvg[1] <- d1$Daily.deaths[1]
for (t1 in 2:6) {
 d1$Deaths.7dAvg[t1] <- mean(d1$Daily.deaths[1:t1])
}
for (t in 7:nrow(d1)) {
 d1$Deaths.7dAvg[t] <- mean(d1$Daily.deaths[(t-6):t])
}

#------------------------------------------------------------------------------#
# Import and reformat the CDC seroprevalence data
#------------------------------------------------------------------------------#

# Reference: https://covid.cdc.gov/covid-data-tracker/?CDC_AA_refVal=https%3A%2F%2Fwww.cdc.gov%2Fcoronavirus%2F2019-ncov%2Fcases-updates%2Fcommercial-labs-interactive-serology-dashboard.html#serology-surveillance
cdc_seroprev.0 <- read.csv('./Data/CDC_seroprevalence.csv')

# Select the CT data
cdc_seroprev <- cdc_seroprev.0[which(cdc_seroprev.0$location=="CT"),]

# Reformat the date variables
cdc_seroprev$start <- as.Date(cdc_seroprev$start, format="%Y/%m/%d")
cdc_seroprev$end   <- as.Date(cdc_seroprev$end,   format="%Y/%m/%d")

# Calculate the middle day of each survey round
cdc_seroprev$mid <- cdc_seroprev$start + floor((cdc_seroprev$end-cdc_seroprev$start)/2)

# Calculate the weighted number of positive samples in each round
cdc_seroprev$n_pos <- round(cdc_seroprev$n_sample * cdc_seroprev$median/100)

