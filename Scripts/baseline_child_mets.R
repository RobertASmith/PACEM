# Script to read in data on physical activity behaviours
# convert minutes to met-mins
# and then use this to create a synthetic population
# by age, sex and IMD.
rm(list = ls())

library(foreign)

# Manual Parameter Inputs #
unfeasible <- 10000
pemets <- 5             # METs from PE assumed to be 6

# input parameters list..
df_mets_lookup <- read.csv("Data/Cleaned/activity_mets.csv")

# variable list, include some generic variables and all those in the activity list...
varlist <- c("SerialA","Age35g","Sex","nssec8","qimd","wemwbs",
             "Urban14b","eqv5","wt_int","wt_child","SampType",
             "Ag015g4","sch15","totalPA15","Normal","BMI",
             "sport08","InfAct08x","wlktot08",
             df_mets_lookup$Activity)
# read in data from hse 2015
df_hse <- foreign::read.dta(
  "Data/Cleaned/hse2015ai.dta",
  convert.factors = FALSE)  %>%
  select(varlist) %>%   # select the variables listed above.
  filter(Age35g == 5, # only those aged 11 and 12
         sch15 >= 0 &
         totalPA15 >= 0 &
         totalPA15 <= unfeasible &
         Normal == 3) %>%
  mutate(
    Sex = replace(
      x = Sex,
      list = which(Sex == 2),
      values = 0
    ),
    WlkScWT = replace(
      x = WlkScWT,
      list = which(WlkScWT < 0),
      values = 0
    )
  )

#====  Estimate of METs per individual  ====#
# 1.  Using compendium of METs from BUTTE et al. 2018, create a dataframe linking activities and METs.
# 2.  Create a matrix of individuals (rows) and activities (columns)
# 3.  Use matrix multiplication (%*%) to calculate the total physical activity (met-mins/wk) for each individual.


#1. Set up a 34*2 Matrix of activities and associated METs from excel sheet.
# make METs from P.E. equal to amount specified in input parameters above.
METlevels <- df_mets_lookup %>% 
  #as.data.frame() %>% 
  #filter(Activity %in% temp) %>%
  mutate(METs = replace(x = METs,
                        list =which(Activity %in% 
                                      c("schMonMVPA","schTueMVPA","schWedMVPA",                                                                       
                                        "schThurMVPA","schFriMVPA","schSatMVPA",
                                        "schSunMVPA")),
                        values = 0)) # removed due to double counting



#2. Set up an empty matrix with a column for each activity type and a row for each individual.
Activities <- matrix(ncol = length(METlevels$Activity), nrow = nrow(df_hse)) 

# Insert the data from hse2015 for each activity - the number in each...
for (i in 1:(length(METlevels$Activity))){
  Activities[,i]<- df_hse[ , METlevels$Activity[i]] #eval(parse(text = paste0(text="df_hse[,", METlevels$Activity[i])))
}

#3. Multiply the minutes in each activity by the mets associated with it...
# returns a vector of METs
df_hse$METmins <- as.numeric(Activities %*% 
  METlevels$METs)

# end up with an unweighted distribution of PA for males and for females...


### LSOA DATA - 1 is most deprived.

imd_2015 <- read.csv(file = "Data/Cleaned/IMD_2015.csv",
                     header = TRUE) %>% 
  as.data.frame() %>% 
  rename(lsoa_code ="LSOA.code..2011.",
         lsoa_name ="LSOA.name..2011.",
         lad_code = "Local.Authority.District.code..2013.",
         lad_name = "Local.Authority.District.name..2013.",
         imd_score = "Index.of.Multiple.Deprivation..IMD..Score",
         imd_decile = "Index.of.Multiple.Deprivation..IMD..Decile..where.1.is.most.deprived.10..of.LSOAs.",
         n_depchild = "Dependent.Children.aged.0.15..mid.2012..excluding.prisoners.",
         n_pop = "Total.population..mid.2012..excluding.prisoners.") %>% 
  select(lsoa_code,lsoa_name,lad_code,lad_name,imd_score,imd_decile,n_depchild,n_pop) %>% 
  mutate(imd_quintile = round((imd_decile+0.5)/2,digits=0))

#==== 
# PLOT DEPRIVATION IN HSE VS ONS.
#====
# deprivation higher for children,
# makes sense because having children makes you poor :)
ggplot(data = imd_2015)+ 
  geom_histogram(aes(x = imd_decile, 
                     y = ..density.., 
                     weight = n_depchild),
                 col = "red",
                 fill = "red",
                 alpha = 0.2) +
  geom_histogram(aes(x = imd_decile, 
                     y = ..density.., 
                     weight = n_pop),
                 col = "blue",
                 fill = "blue",
                 alpha=0.2) +
  theme_minimal()

# Step 4. Use the ONS data to weight the HSE sample.
# In order to get back to a representative sample we need to 
# overweight those in more deprived areas, 
# in particular males. 
# We assume 50% male and female in every area but this could be changed.
df_props <- imd_2015 %>% 
  group_by(imd_quintile) %>% 
  summarise(pop = sum(n_depchild)) %>%
  mutate(ons_male = pop/sum(pop)/2,
         ons_fmle = pop/sum(pop)/2)

# creates proportions of HSE data for 10 groups.
df_props$hse_fmle <- (table(df_hse$qimd,df_hse$Sex)/sum(table(df_hse$qimd,df_hse$Sex)))[,1]
df_props$hse_male <- (table(df_hse$qimd,df_hse$Sex)/sum(table(df_hse$qimd,df_hse$Sex)))[,2]

# creates weights to 'rebalance'.
weights <- df_props %>% mutate(wt_male = ons_male/hse_male,
                               wt_fmle = ons_fmle/hse_fmle) %>%
  select(wt_male,wt_fmle)

# use a loop to fill in weights.
for(i in 1:nrow(df_hse)){
  
  df_hse$SampType[i] <- as.numeric(weights[df_hse$qimd[i], df_hse$Sex[i]+1])
}

# store information in clean data
df_output <- df_hse %>% 
  select(Sex, qimd, SampType, BMI, METmins) %>% 
  rename(female = Sex, child_wt = SampType) %>% 
  as.data.frame

# write the data to the cleaned file ready to be read in by the model
write.csv(x = df_output,
          row.names = F, 
          file = "Data/Cleaned/hse15_mets.csv")
