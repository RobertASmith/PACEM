# Main script.
rm(list = ls())

# source external packages
library(dplyr)
library(ggplot2)
library(miceadds)

# source all functions for the model
miceadds::source.all("R")

# PARAMETER INPUTS
run_sens   <- T # whether to run the sensitivity analyses
run_main   <- F # whether to run the main analysis
n_psa      <- 5000    # number of PSA for main
n_psa_sens <- 2000    # number of PSA for sensitivity analysis
n_i        <- 10000   # number of individuals in the model

age_start   <- 11     # start age in model
age_adult   <- 18     # age at which classed as adult
age_end     <- 80     # end age in the model
met_min_max <- 10000  # upper truncation met mins
met_min_min <- 0      # lower truncation met mins

# READ IN DATA ----
## PA Data ----
cost_intervention <- 57

# intervention effect
v_pa_decay_baseline      <- rnorm(n = n_psa, mean = 0.061, calculate_std_from_ci(lwr = 0.039, upr = 0.083))
#v_pa_decay_intervention  <- v_pa_decay_baseline - 0.061 # not used
v_intervention_additional <- rnorm(n = n_i, mean =  186, calculate_std_from_ci(lwr = -36, upr = 408))

# set random number seeds for the model runs.
v_HSseed <- sample(x = 10000, size = n_psa, replace = F)

# Dataframe containing, in each row, the percentile of PA corresponding
# to each age & sex (in column).
df_PA_percentiles <- read.csv("Data/Cleaned/PA_Percentiles", 
                              row.names = 1)

# Baseline characteristics: Sex (1= female), IMD, and Initial METs
# randomly sample from weights:
m_fc <- dplyr::sample_n(tbl = read.csv(file = "Data/Cleaned/hse15_mets.csv"),
                  replace = TRUE,
                  size = n_i,
                  weight = child_wt) %>% 
  rename(init_mets = METmins) %>%
  mutate(age = age_start,
         init_mets = ifelse(test = init_mets<=0, yes =  10, no = init_mets),
         init_mets = ifelse(test = init_mets>10000, yes =  10000, no = init_mets)) %>% 
  select(female, age, qimd, init_mets)


# Visual - Create the plot for the initial PA level
plot_init_dist <- ggplot()+
  geom_density(data = m_fc, 
               mapping = aes(x = init_mets, col = as.factor(female)))+
  geom_vline(data = m_fc,
             aes(xintercept = 2520, col = as.factor(female)),
             linetype = 2) +
  facet_wrap(~female,
             labeller = as_labeller(x = c("1" = "Female", "0" = "Male")))+
  theme_minimal()+
  scale_color_manual(values = c("1" = "#AA336A", "0" = "#00008B"), 
                     guide = "none")+
  scale_y_continuous(name = "Density",labels = NULL)+
  scale_x_continuous(name = "Metabolic Equivalent of Task minutes per week",
                     limits = c(0,9500))

# calculate guideline proportions...
guideline_props <- m_fc %>% 
  group_by(female) %>% 
  summarise(prop_guidelines = mean(init_mets >= 2520))

# save the file in 'figures'
ggsave(filename = "Figures/Initial_Dist.png",
       plot = plot_init_dist, 
       width = 12, height = 8,
       scale = 0.7)


## HS data ----
# read in data on relative risks.
df_rr2 <- readRDS("Data/Cleaned/rr.rds")

# create a different set of parameters for each PSA run
m_c <- sapply(
  X = unique(df_rr2$cond),
  FUN = function(x, n, data, value) {
    mean <- getRRmean(data = data,
                      condition = x,
                      value = value)
    sd   <- getRRsd(data = data,
                    condition = x,
                    value = value)
    rnorm(n = n, mean = mean, sd = sd)
  },
  data = df_rr2,
  value =  "alpha",
  n = n_psa
)
m_ln <- sapply(
  X = unique(df_rr2$cond),
  FUN = function(x, n, data, value) {
    mean <- getRRmean(data = data,
                      condition = x,
                      value = value)
    sd   <- getRRsd(data = data,
                    condition = x,
                    value = value)
    rnorm(n = n, mean = mean, sd = sd)
  },
  data = df_rr2,
  value =  "beta",
  n = n_psa
)

# read in data on the incidence of each condition.
a_inc <- readRDS(file = "Data/Cleaned/GBD_2019_incidence.rds")

# read in the data on the condition mortality for each condition and OCM
df_cond_mort <-  read.csv(file = paste0("Data/Cleaned/Cond_mort"), 
                          row.names = 1) %>% 
                  replace(is.na(.), 0)

# create an array of NAs for condition 
a_cond_mort <- array(data = NA,
                     dim = c(dim(df_cond_mort), 1000))
dimnames(a_cond_mort) <- dimnames(df_cond_mort)

# input the baseline matrix into every slice
a_cond_mort[,,] <- as.matrix(df_cond_mort)
# add in random variation of 5% to every mortality rate ...
# random variation occurs by slice & condition, not by time!!
a_cond_mort[,1:5,] <- apply(a_cond_mort[,1:5,],
                            MARGIN = c(2, 3),
                            function(x)
                              x * runif(n = 1,
                                        min = 0.95,
                                        max = 1.05))

# create a set of durations, sampled at random, some uncertainty here around 10.
df_durations <- data.frame(
  "BC" = triangle::rtriangle(n = n_psa, a = 8, b = 12, c = 10),
  "CC" = triangle::rtriangle(n = n_psa, a = 8, b = 12, c = 10),
  "T2D" = 100,
  "IS" = 100,
  "IHD" = 100
)

## Econ data ----

# discount rate
dr_default <- 0.035

# read in the matrix of healthy health state quantiles derived from Ara & Brazier.
m_no_HS_quantiles <- readRDS("Data/Cleaned/AB_2011_Utilities.rds")$m_no_HS_quantiles

# Assign utility decrement multipliers from Gc et al. 2019.
df_u_multipliers <- data.frame(
  IHD = rbeta(n = n_psa, shape1 = 357, shape2 = 191),
  IS  = rbeta(n = n_psa, shape1 = 355, shape2 = 323),
  T2D = rbeta(n = n_psa, shape1 = 5032, shape2 = 2548),
  BC  = rbeta(n = n_psa, shape1 = 791, shape2 = 256),
  CC  = rbeta(n = n_psa, shape1 = 150, shape2 = 73)
) / 0.8 # for general population utilities...


# Cost inputs from Gc, V.S., et al., (2019). 
# Cost-effectiveness of physical activity interventions in adolescents: 
# model development and illustration using two exemplar interventions. 
# BMJ open, 9(8), p.e027566.
df_first_year_costs <-
  data.frame(
    "BC" =  rgamma(n = n_psa, shape = 100, scale = 122),
    "CC" =  rgamma(n = n_psa, shape = 100, scale = 170),
    "T2D" = rgamma(n = n_psa, shape = 100, scale = 13),
    "IS" =  rgamma(n = n_psa, shape = 100, scale = 101),
    "IHD" = rgamma(n = n_psa, shape = 100, scale = 56)
  ) * 1.1058 # uplift from 2013 to 2020 from https://www.pssru.ac.uk/pub/uc/uc2021/sourcesofinformation.pdf
df_subsequent_year_costs <- data.frame(
  "BC" = 0,
  "CC" = 0,
  "T2D" = rgamma(n = n_psa, shape = 100, scale = 13),
  "IS" = rgamma(n = n_psa, shape = 100, scale = 27),
  "IHD" = rgamma(n = n_psa, shape = 100, scale = 2)
)  * 1.1058 # uplift from 2013 to 2020 from https://www.pssru.ac.uk/pub/uc/uc2021/sourcesofinformation.pdf


# SINGLE RUN - CODE CHECKER ----
## Single Iteration run as a check to identify problems before running full PSA
x <- 2
results_mult <- runDeterm(
  output = "epi",
  n_i = n_i,
  v_init_mets = m_fc$init_mets,
  pa_decay_baseline = v_pa_decay_baseline[x],
  pa_decay_intervention = v_pa_decay_baseline[x],
  pa_additional_intervention = v_intervention_additional,
  df_PA_percentiles = df_PA_percentiles,
  age_start = age_start,
  age_adult = age_adult,
  age_end = age_end,
  v_female = m_fc$female,
  effect_decay_type = "exp",
  effect_decay_exp_rate = 0.5,
  effect_decay_linear_duration = 10,
  effect_decay_proportion = 2/10,
  v_constant = unlist(m_c[x, ]),
  v_ln = unlist(m_ln[x, ]),
  df_inc = a_inc[, , sample(1:1000, size = 1, replace = T)] ,
  v_durations = unlist(df_durations[x, ]),
  df_cond_mort = df_cond_mort,
  utility_method = "multiplicative",
  v_u_hp = m_no_HS_quantiles[, 500], #df_U$no_HS,
  v_first_year_costs = unlist(df_first_year_costs[x, ]),
  v_subsequent_year_costs = unlist(df_subsequent_year_costs[x, ]),
  v_HS_mult = unlist(df_u_multipliers[x,]),
  dr_costs = dr_default,
  dr_qalys = dr_default,
  HSseed = v_HSseed[x]
)

# take results and add the age and strategy
df_results <- as.data.frame(do.call("rbind", results_mult))
df_results$age <- 1:age_end
df_results$strategy <- rep(c("base", "int"), each = age_end)

# plot using ggplot
ggplot(data = df_results %>% 
         filter(age > 12) %>% 
         tidyr::pivot_longer(cols = -c(strategy, age)))+
  geom_line(aes(y = value, x = age, col = strategy))+
  facet_wrap(~name, scales = "free")

# DEFINE SENSITIVITY ANALYSIS ----

# create scenarios in a single data frame
# these define the scenarios run as sensitivity analyses
df_sensScenarios <-
  data.frame(
    name = c(
      "base",
      "habit_lower",
      "habit_upper",
      "decay_lin",
      "utils_minimum",
      "dr_1_5",
      "high_effect",
      "low_effect"
    ),
    label = c(
      "Baseline",
      "0% develop habit",
      "100% develop habit",
      "Linear decay",
      "Minimum HRQoL method",
      "1.5% Discount Rate",
      "279 METs",
      "93 METs"
    ),
    habit_prop = c(0.1, 0, 1, rep(0.1, 5)),
    decay_method = c(rep("exp", 3), "lin", rep("exp", 4)),
    utils_method = c(rep("multiplicative", 4), "minimum", rep("multiplicative", 3)),
    d_r = c(rep(dr_default, 5), 0.015, rep(dr_default, 2)),
    effect = c(rep(1, 6), 1.5, 0.5),
    n_psa = c(n_psa, rep(n_psa_sens, 7)),
    run = c(run_main, rep(run_sens, 7))
  )

# write a data frame of the scenarios for reference
write.csv(x = df_sensScenarios, 
          file = "Results/df_sensScenarios.csv")
 
# RUN MODEL ----
# Loop through the scenarios, running the model for each.
# store the results for each as a seperate .RDS file in 'Results'
for(scenario_index in 1:nrow(df_sensScenarios)){
  # initialise arrays  
  a_PSA_base <- a_PSA_int <- array(data = NA,
                                   dim = c(80,
                                           18,
                                           df_sensScenarios[scenario_index, "n_psa"]))
  dimnames(a_PSA_base) <- dimnames(a_PSA_int) <- list(1:80,
                                                      colnames(results_mult$base),
                                                      1:df_sensScenarios[scenario_index, "n_psa"])
 # get name of scenario
 scenario_name <- df_sensScenarios[scenario_index, "name"]
 # determine whether to run scenario, if skip then let me know
 if (!df_sensScenarios[scenario_index, "run"]) {
   message(paste0("\n Skipped over: ", scenario_name))
   next
 }
 # tell me running scenario
 message(paste0("\n Running: ", scenario_name))
 # initialise progress bar for scenario.
 pb <- txtProgressBar(min = 1, max = n_psa, style = 3)
 # run through each PSA iteration ... sampling from distributions and running the model
for(x in 1:df_sensScenarios[scenario_index, "n_psa"]){
 # change progress bar value
 setTxtProgressBar(pb, x)

 # run deterministic model with parameter inputs.
 temp <- runDeterm(
   output = "epi",
   n_i = n_i,
   v_init_mets = m_fc$init_mets,
   pa_decay_baseline = v_pa_decay_baseline[x],
   pa_decay_intervention = v_pa_decay_baseline[x], #v_pa_decay_intervention[x],
   pa_additional_intervention = v_intervention_additional * df_sensScenarios[scenario_index, "effect"],
   df_PA_percentiles = df_PA_percentiles,
   age_start = age_start,
   age_adult = age_adult,
   age_end = age_end,
   v_female = m_fc$female,
   effect_decay_type = df_sensScenarios[scenario_index, "decay_method"],
   effect_decay_exp_rate = 0.5,
   effect_decay_linear_duration = 10,
   effect_decay_proportion = 1 - df_sensScenarios[scenario_index, "habit_prop"],
   v_constant = unlist(m_c[x, ]),
   v_ln = unlist(m_ln[x, ]),
   df_inc = a_inc[, , sample(1:1000, size = 1, replace = T)],
   v_durations = unlist(df_durations[x, ]),
   df_cond_mort = a_cond_mort[,,sample(1:1000, size = 1, replace = T)],
   utility_method = df_sensScenarios[scenario_index, "utils_method"],
   v_u_hp = m_no_HS_quantiles[, sample(1:1000, size = 1, replace = T)],
   v_first_year_costs = unlist(df_first_year_costs[x, ]),
   v_subsequent_year_costs = unlist(df_subsequent_year_costs[x, ]),
   v_HS_mult = unlist(df_u_multipliers[x,]),
   dr_costs = df_sensScenarios[scenario_index, "d_r"],
   dr_qalys = df_sensScenarios[scenario_index, "d_r"],
   HSseed = v_HSseed[x]
 )
 # store matrices in arrays.
 # slow - speed up
 a_PSA_base[,,x] <- temp$base
 a_PSA_int[,,x] <- temp$int
}

# create a list of results.
l_results <- list(
 "int" = a_PSA_int,
 "base" = a_PSA_base,
 "params" = c(
   "n_psa" = df_sensScenarios[scenario_index, "n_psa"],
   "n_i" = n_i,
   "age_start" = age_start,
   "age_adult" = age_adult,
   "age_end" = age_end,
   "dr_costs" = df_sensScenarios[scenario_index, "d_r"],
   "dr_qalys" = df_sensScenarios[scenario_index, "d_r"],
   "cost_intervention" = cost_intervention,
   "effect multiplier" = df_sensScenarios[scenario_index, "effect"]
 )
)

# save as object for later.
saveRDS(object = l_results,
        file = paste0("Results/", scenario_name, ".rds"))
}
      

  

