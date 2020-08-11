# Lawlor and Arellano 2020
# Scientific Reports

# Script 1. Prepare and Clean code

# 1. Call libraries
#-------------------------------------
library(tidyverse)
library(PNWColors)
library(here)
library(broom)
library(akima)
library(ggrepel)
library(patchwork)
library(egg)
#----


# Upload data
#-------------------------------------
growth <- read.csv(here::here("data_raw","GrowthData_Master.csv")) #6586 rows
mort <- read.csv(here::here("data_raw","OlyMortality2018.csv"))
growthchem <- read.csv(here::here("data_raw","GrowthChemistry_Master_DIC_full.csv"))

#----


# find day when each cup reaches 70% competence
#-------------------------------------
days.before.comp <-
  growth %>%
  group_by(Cup,Date) %>%
  arrange(rev(Cup),Date) %>%
  summarise(N = n(), # find number larvae sampled in that cup / day
            Eyes = sum(Eyes=="y"), # find number of larvae in sample with eyes
            Eyes = round(Eyes/N * 100,digits=2)) %>% # find % eyed larvae
  ungroup() %>% 
  group_by(Cup) %>%
  slice(if(any(Eyes>=70)) 1:which.max(Eyes >= 70) else row_number()) %>% # cut off at first day > 70%
  summarise(days.before.comp = n()) %>% # find number of sample days before > 70% competence
  as.data.frame()
#----



# find day when cups are >95% mortality
#-------------------------------------
days.before.mort <-
  mort %>% 
  mutate(calcmort = round(Dead/(Alive+Dead)*100,digits=2)) %>%
  mutate(calcmort = ifelse(is.na(calcmort), 95, calcmort)) %>%# change NA values to 95% mortality (NA was used when cups were far too dead to count)
  group_by(Cup) %>%
  arrange(Cup, Date) %>%
  select(-Mortality,-Eyed,-Counter,-Percent.eyed) %>%
  slice(if(any(calcmort>=95)) 1:which.max(calcmort >= 95) else row_number()) %>% # find first day > 95% mortality
  summarise(days.before.mort=n()) %>%
  select(Cup,days.before.mort) %>%
  as.data.frame()
# we are missing one day of data for cup D5, 
# so it's saying D5 lasted 6 days instead of 7
# For now, we'll just manually correct:
days.before.mort[days.before.mort$Cup %in% c("D5"),2] <- 7
#----



# merge "days until comp" and "days until mort" dfs, add column for whichever is lower
#-------------------------------------
dates <- growth %>% arrange(Date) %>% distinct(Date) # vector of all sampling dates
daylist <- c(0,3,5,7,10,12,14,17) # vector of days of experiment
sampledays <- right_join(days.before.comp,days.before.mort,by="Cup") %>%
  rowwise() %>%
  mutate(sampledays = min(days.before.comp,days.before.mort)) %>% # add day that we cut off the cup (becasue either high mort or high competence)
  rowwise() %>%
  mutate(enddate = sapply(dates, "[", sampledays)) %>%
  mutate(sampleday = daylist[sampledays]) %>%
  as.data.frame()
# remove merged dfs
rm(days.before.comp,days.before.mort,dates,daylist)
#----


# trim growth and mortality datasets after cutoff date
#-------------------------------------
growth_trimmed <- inner_join(growth,sampledays,by="Cup") %>%
  arrange(Cup,Date) %>%
  group_by(Cup) %>%
  filter(Date <= enddate) %>%
  select(-enddate,-sampledays) %>%
  as.data.frame() %>%
  mutate(day = (Date-050418)/100)

mort_trimmed <- inner_join(mort,sampledays,by="Cup") %>%
  arrange(Cup,Date) %>%
  group_by(Cup) %>%
  filter(Date <= enddate) %>%
  select(-enddate,-Mortality) %>%
  as.data.frame()

rm(growth, mort)
#----



# order dfs in order of cup target values
#-------------------------------------
cuporder <- c("I4","J5","I1","J2","I5","J1","I2","I3","J3","J4",
              "G1","H1","G3","G2","H3","H2","G5","G4","H5","H4",
              "F2","F5","F3","E5","E4","E3","E1","E2","F4","F1",
              "C1","C2","D3","C4","D4","C5","D5","D2","C3","D1",
              'B2','B4','B1','B5',"A3","A1","A5","A4","B3","A2")

sampledays <- sampledays %>%
  mutate(Cup_f = factor(Cup,levels=cuporder))

growth_trimmed <- growth_trimmed %>%
  mutate(Cup_f = factor(Cup, levels = cuporder))

mort_trimmed <- mort_trimmed %>%
  mutate(Cup_f = factor(Cup,levels=cuporder))
#----




# make growth summary (sizes per day)
#-------------------------------------
growthsummary <- growth_trimmed %>%
  group_by(Cup,day,CupNumber,Date,Cup_f) %>%
  
  # find average size per sample day
  summarise(N = n(),
            Eyes = sum(Eyes=="y"),
            Eyes = round(Eyes/N * 100,digits=2),
            um2 = mean(um),
            sd = sd(um),
            se = sd/sqrt(N)
  ) %>%
  
  # rename reokace old size column with new
  rename(um = um2) %>%
  
  # change cup_f to ordered (cup_f = cup as factor)
  ungroup() %>%
  mutate(Cup_f = factor(Cup_f, levels=cuporder))
#----




# find time to success (>25% competence)
#-----------------------------------------
## make a df of time until success (>25% eyes for each cup)
successtime <- c()
for(i in 1:50){
  df <- growthsummary[growthsummary$Cup %in% unique(growthsummary$Cup)[i],]
  df <- df[order(df$Date),] 
  
  if(any(df$Eyes >= 25)){
    successtime <<- cbind(successtime,df$day[min(which(df$Eyes >= 25))])
  }
  else{
    successtime <<- cbind(successtime,NA)
  }
}
successtime <- as.vector(successtime) # vector of days to success
# make DF of first day of success for each cup
successdf <- data.frame(successtime,unique(growthsummary$Cup_f))
colnames(successdf) <- c("successtime","Cup_f")
rm(df,successtime,i)
#----



# make growth regressions for each cup
#--------------------------------------
regressions <- 
  growthsummary %>%
  
  # trim to relevant columns
  select(Cup,day,um) %>%
  group_by(Cup) %>%
  
  # collapse
  nest() %>%
  
  mutate(
    
    # fit regression model for each cup
    fit = map(data, ~ lm(um ~ 0 + day,  # 0 forces through intercept
                         offset = rep(156.05,length(day)), # offset accounts for starting size (156.05um)
                         data = .x)),
    
    # expand regressions
    tidied = map(fit, broom::tidy),
    
    # pull r2
    summary = map(fit,summary),
    r_sq = map_dbl(summary, "r.squared"),
  )  %>%
  
  ungroup() %>%
  
  # get summary stats
  unnest(tidied) %>%
  
  # remove unneeded rows
  select(-data,-fit,-summary,-statistic,-term) %>%
  
  # round to 2 decimals
  mutate_at(vars(-Cup), .funs= list(~round(.,2))) %>%
  
  # add starting size
  mutate(intercept = 156.05) %>%
  
  # rename
  rename(Gr=estimate,p = p.value)%>%
  
  # make factor to order
  mutate(Cup_f = factor(Cup,levels=cuporder)) %>%
  as.data.frame()
#----



# Clean chemistry data
#--------------------------------------
growthchem <- inner_join(growthchem,sampledays,by="Cup") %>%
  arrange(Cup,Date) %>%
  group_by(Cup) %>%
  filter(Date <= enddate) %>%
  select(-enddate,-sampledays,-sampleday,-days.before.comp,-days.before.mort) %>%
  as.data.frame() %>%
  mutate(day = (Date-050418)/100)

# make df of mean chem trts
growthchem_means <- growthchem %>% 
  group_by(Cup) %>% 
  summarise_at(vars(Temp,Sal,pH,pCO2,Ar,DIC), list(mean,sd)) %>%
  as.data.frame() %>%
  mutate_at(2:13, round, 2) %>%
  mutate(Cup_f = factor(Cup,levels=cuporder))
colnames(growthchem_means) <- c("Cup","Temp","Sal","pH","pCO2","Ar","DIC","Tempsd","Salsd","pHsd","pCO2sd","Arsd","DICsd","Cup_f")
#----





save(growth_trimmed,file = here::here("data_processed","growth_trimmed.Rdata"))
save(mort_trimmed,file = here::here("data_processed","mort_trimmed.Rdata"))
save(regressions, file = here("data_processed","regressions.Rdata"))
save(growthsummary, file=here("data_processed","growthsummary.Rdata"))
save(sampledays, file = here("data_processed","sampledays.Rdata"))
save(successdf, file = here("data_processed","successdf.Rdata"))
save(growthchem, file = here("data_processed","growthchem.Rdata"))
save(growthchem_means, file = here("data_processed","growthchem_means.Rdata"))













