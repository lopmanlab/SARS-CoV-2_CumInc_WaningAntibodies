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

# Downloaded from here https://www1.nyc.gov/site/doh/covid/covid-19-data-deaths.page
# at 9:00AM (ET) on Friday, October 2, 2020
d <- read.csv('./Data/data-sS9Wp.csv')

# Create a copy of the original dataset
d1 <- d

# Reformat a date variable
d1$date <- as.Date(d1$date_of_interest, format="%m/%d/%Y")

# Calculate the total deaths (confirmed + probable)
d1$Combined <- d1$DEATH_COUNT + d1$DEATH_COUNT_PROBABLE

#------------------------------------------------------------------------------#
# Import and reformat the CDC seroprevalence data
#------------------------------------------------------------------------------#

# Reference: https://covid.cdc.gov/covid-data-tracker/?CDC_AA_refVal=https%3A%2F%2Fwww.cdc.gov%2Fcoronavirus%2F2019-ncov%2Fcases-updates%2Fcommercial-labs-interactive-serology-dashboard.html#serology-surveillance
cdc_seroprev.0 <- read.csv('./Data/CDC_seroprevalence.csv')

# Select the NYC data
cdc_seroprev <- cdc_seroprev.0[which(cdc_seroprev.0$location=="NYC"),]

# Reformat the date variables
cdc_seroprev$start <- as.Date(cdc_seroprev$start, format="%Y/%m/%d")
cdc_seroprev$end   <- as.Date(cdc_seroprev$end,   format="%Y/%m/%d")

# Calculate the middle day of each survey round
cdc_seroprev$mid <- cdc_seroprev$start + floor((cdc_seroprev$end-cdc_seroprev$start)/2)

# Calculate the weighted number of positive samples in each round
cdc_seroprev$n_pos <- round(cdc_seroprev$n_sample * cdc_seroprev$median/100)
