library(neonUtilities)
library(tidyverse)
library(lubridate)
library(ggplot2)

# NEON site names
sites <- c("SERC", "ORNL", "OSBS", "SCBI")

# download data from neon
tick.data <- loadByProduct("DP1.10093.001",                  # Ticks sampled using drag cloths
                           site = sites,                     # specify sites
                           check.size = FALSE,               # so I dont have to say yes to downloading data
                           token = Sys.getenv("NEON_TOKEN")) # using personal neon api token, set in .Renviron
# note: data will download without a token, but rate of download will be limited

save(tick.data, file = "data/NEON_tick_NEFI2020.RData")      # save for later usage


load("data/NEON_tick_NEFI2020.RData")       # load tick.data
dat.field <- tick.data[["tck_fielddata"]]   # get field data
dat.tax <- tick.data[["tck_taxonomyProcessed"]] # get taxonomy data

dat.tax$siteID <- gsub("[[:punct:]]\\d*", "", dat.tax$plotID) # need to create siteID column
spp <- "Amblyomma americanum" # most abundant species
dat.aa <- dat.tax %>%       # use records that have Amblyomma americanum ID
  filter(scientificName == spp)

# only adult ticks have a sex ID
# add age column (convert male and female to "Adult")
dat.aa$age <- gsub("Male", "Adult", dat.aa$sexOrAge)
dat.aa$age <- gsub("Female", "Adult", dat.aa$age)

# ===================
#     Nymph Data 
# ===================

# need to retrieve zeros from field data
nymph.0 <- dat.field %>% 
  mutate(day = format(collectDate, "%Y-%m-%d")) %>% # create day column (drop hours:min:sec)
  group_by(siteID) %>%              # group by site first
  group_by(day, .add = TRUE) %>%    # group by sampling day within site (some plots sampled on same day)
  summarise(totalNymph = sum(nymphCount)) %>%  # total nymphs counted each day at each site (across plots)
  filter(totalNymph == 0)        # grab zero rows
 
# Nymph at the site level
aa.nymph <- dat.aa %>% 
  filter(age == "Nymph") %>%        # grab nymphs
  mutate(day = format(collectDate, "%Y-%m-%d")) %>%   # drop hours:min:sec
  group_by(siteID) %>%              # need to group by site first
  group_by(day, .add = TRUE) %>%    # group by sampling day within site (some plots sampled on same day)
  summarise(totalNymph = sum(individualCount))      # sum over sampling day 

# all AA nymph data sorted by date and into train / validate
nymph.data <- bind_rows(aa.nymph, nymph.0) %>% 
  arrange(day) %>%  
  mutate(trainValidate = ifelse(day < ymd("2018-05-01"), "Train", "Validate")) %>% 
  unite(col = siteDay, siteID, day, sep = "_", remove = FALSE) # create siteDay column

# make axis labesl and plot
x.seq <- nymph.data %>% 
  filter(siteID == "OSBS") %>% 
  filter(trainValidate == "Train") %>% 
  select(day)

at <- round(seq(1, nrow(x.seq), length.out = 10))

ggplot(nymph.data %>% filter(trainValidate == "Train"), aes(day, totalNymph)) +
  geom_point(size = 3) +
  facet_grid(rows = vars(siteID),scales = "free") +
  scale_x_discrete(breaks = x.seq$day[at]) +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "A. americanum nymphs",
       y = "Individuals Observed",
       x = "Date") 


# ===================
#     Adult Data 
# ===================

# need to retrieve zeros from field data
adult.0 <- dat.field %>% 
  mutate(day = format(collectDate, "%Y-%m-%d")) %>% # create day column (drop hours:min:sec)
  group_by(siteID) %>%              # group by site first
  group_by(day, .add = TRUE) %>%    # group by sampling day within site (some plots sampled on same day)
  summarise(totalAdult = sum(adultCount)) %>%  # total nymphs counted each day at each site (across plots)
  filter(totalAdult == 0)        # grab zero rows

# Nymph at the site level
aa.adult<- dat.aa %>% 
  filter(age == "Adult") %>%        # grab nymphs
  mutate(day = format(collectDate, "%Y-%m-%d")) %>%   # drop hours:min:sec
  group_by(siteID) %>%              # need to group by site first
  group_by(day, .add = TRUE) %>%    # group by sampling day within site (some plots sampled on same day)
  summarise(totalAdult = sum(individualCount))      # sum over sampling day 

# all AA nymph data sorted by date and into train / validate
adult.data <- bind_rows(aa.adult, adult.0) %>% 
  arrange(day) %>%  
  mutate(trainValidate = ifelse(day < ymd("2018-05-01"), "Train", "Validate"))%>% 
  unite(col = siteDay, siteID, day, sep = "_", remove = FALSE)

# make axis labesl and plot
x.seq <- adult.data %>% 
  filter(siteID == "OSBS") %>% 
  filter(trainValidate == "Train") %>% 
  select(day)

at <- round(seq(1, nrow(x.seq), length.out = 10))

ggplot(adult.data %>% filter(trainValidate == "Train"), aes(day, totalAdult)) +
  geom_point(size = 3) +
  facet_grid(rows = vars(siteID),scales = "free") +
  scale_x_discrete(breaks = x.seq$day[at]) +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "A. americanum adults",
       y = "Individuals Observed",
       x = "Date") 






# ===================
#     Weather Data 
# ===================


# ran the download and formating on BUs computing cluster, 
# download was taking forever on local machine
# ==============================================================================
met.ids <- c("DP1.00006.001", # precipitation
             "DP1.00098.001", # relative humidity
             "DP1.00003.001") # air temperature

met.data <- list()
for(i in seq_along(met.ids)){
  met.data[[i]] <- loadByProduct(met.ids[i],
                                 site = sites,
                                 check.size = FALSE,
                                 token = Sys.getenv("NEON_TOKEN"))
  
}
save(met.data, file = "NEON_Met_NEFI2020.RData")


# loadByProduct() and save() ran as a batch job, loading output
load("NEON_Met_NEFI2020.RData")

# extract raw data and remove large list
precip <- met.data[[1]]
relativeHumidity <- met.data[[2]]
airTemperature <- met.data[[3]]
rm(met.data)


## precipitation
vars <- precip$variables_00006 %>%  # look at variables for 30 minute data
  filter(table == "THRPRE_30min")
through.30 <- precip$THRPRE_30min  # extract 30 min data

# same general workflow as ticks above
daily.precip <- through.30 %>% 
  mutate(day = format(through.30$endDateTime, "%Y-%m-%d")) %>% 
  group_by(siteID) %>% 
  group_by(day, add = TRUE) %>% 
  filter(TFPrecipExpUncert != -9999) %>%   # remove days with -9999 as percent uncertainty
  summarise(dailyTotalPrecip = sum(TFPrecipBulk), # total precip on each day
            dailyPrecipExpUncert = mean(TFPrecipExpUncert, na.rm = TRUE)) # mean uncertainty on each day

x.seq <- daily.precip %>% 
  filter(siteID == "SERC") %>% 
  select(day)

at <- round(seq(1, nrow(x.seq), length.out = 10))

ggplot(daily.precip, aes(day, dailyTotalPrecip)) +
  geom_point(size = 1.5) +
  facet_grid(rows = vars(siteID), scales = "free") +
  scale_x_discrete(breaks = x.seq$day[at]) +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Daily Precipitation",
       y = "mm",
       x = "Date") 



## relative humidity, same workflow as precip
vars <- relativeHumidity$variables_00098 %>% 
  filter(table == "RH_1min")
sensors <- relativeHumidity$sensor_positions_00098
rh.1min <- relativeHumidity$RH_1min # using 1 minute data

daily.rh <- rh.1min %>% 
  mutate(day = format(rh.1min$endDateTime, "%Y-%m-%d")) %>% 
  group_by(siteID) %>% 
  group_by(day, add = TRUE) %>% 
  filter(RHExpUncert != -9999) %>%
  summarise(dailyMaxRH = max(RHMaximum, na.rm = TRUE),  # maximum rh each day
            dailyMinRH = min(RHMinimum, na.rm = TRUE),  # minumun rh each day
            dailyAvgExpUncertRH = mean(RHExpUncert, na.rm = TRUE))

x.seq <- daily.rh %>% 
  filter(siteID == "SCBI") %>% 
  select(day)

at <- round(seq(1, nrow(x.seq), length.out = 10))

ggplot(daily.rh, aes(day, dailyMinRH)) +
  geom_point(size = 1.5) +
  facet_grid(rows = vars(siteID), scales = "free") +
  scale_x_discrete(breaks = x.seq$day[at]) +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Daily Min RH",
       y = "%",
       x = "Date") 

## temperature
vars <- airTemperature$variables_00003
sensors <- airTemperature$sensor_positions_00003
temp.1min <- airTemperature$TAAT_1min
dailyAirTemp <- temp.1min %>% # 1 minute temperature data
  mutate(day = format(temp.1min$endDateTime, "%Y-%m-%d")) %>% 
  group_by(siteID) %>% 
  group_by(day, add = TRUE) %>% 
  filter(tempTripleExpUncert != -9999) %>% 
  summarise(dailyMaxAirTemp = max(tempTripleMaximum, na.rm = TRUE),
            dailyMinAirTemp = min(tempTripleMinimum, na.rm = TRUE),
            dailyAvgExpUncertTemp = mean(tempTripleExpUncert)) %>% 
  filter(!is.infinite(dailyMaxAirTemp))

x.seq <- dailyAirTemp %>% 
  filter(siteID == "SCBI") %>% 
  select(day)

at <- round(seq(1, nrow(x.seq), length.out = 10))

ggplot(dailyAirTemp, aes(day, dailyMaxAirTemp)) +
  geom_point(size = 1.5) +
  facet_grid(rows = vars(siteID), scales = "free") +
  scale_x_discrete(breaks = x.seq$day[at]) +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Daily Max Air Temp",
       y = "Deg. C",
       x = "Date") 


NEFI.met <- list()
NEFI.met$dailyPrecip <- daily.precip
NEFI.met$dailyRH <- daily.rh
NEFI.met$dailyAirTemp <- dailyAirTemp

save(NEFI.met, file = "DailyMetWithUncert.RData")
# after saving NEFI.met, transfer to local machine
# ==============================================================================


# load raw met data after transfering it to local machine
load("data/DailyMetWithUncertNEFI2020.RData")
  
str(NEFI.met)

dailyPrecip <- pluck(NEFI.met, "dailyPrecip") %>% 
  unite(col = siteDay, siteID, day, sep = "_", remove = FALSE) # make a siteDay column for matching datasets
dailyRH <- pluck(NEFI.met, "dailyRH") %>% 
  unite(col = siteDay, siteID, day, sep = "_", remove = FALSE)
dailyTemp <- pluck(NEFI.met, "dailyAirTemp") %>% 
  unite(col = siteDay, siteID, day, sep = "_", remove = FALSE)

# join data by siteDay
nymphData <- left_join(nymph.data, dailyPrecip, by = "siteDay")
nymphData <- left_join(nymphData, dailyRH, by = "siteDay")
nymphData <- left_join(nymphData, dailyTemp, by = "siteDay")
nymphData <- nymphData %>% 
  select(-contains(".")) %>% # remove duplicate columns
  separate(siteDay, c("Site", "Day"), sep = "_") # separate siteDay back to two columns

# join data by siteDay, same workflow as nymphs
adultData <- left_join(adult.data, dailyPrecip, by = "siteDay")
adultData <- left_join(adultData, dailyRH, by = "siteDay")
adultData <- left_join(adultData, dailyTemp, by = "siteDay")
adultData <- adultData %>% 
  select(-contains(".")) %>% 
  separate(siteDay, c("Site", "Day"), sep = "_")

# write data to csv
write.csv(nymphData, file = "nymphDataWithMet.csv")
write.csv(adultData, file = "adultDataWithMet.csv")





