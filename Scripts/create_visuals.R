rm(list = ls())

library(ggplot2)
library(dplyr)

# INCIDENCE & PREVALENCE ----

# simulation results
a_PSA_results <- readRDS("Results/base.rds")$base

# GBD incidence
df_inc <- read.csv(file = "Data/Cleaned/IHME-GBD_2019_DATA-c4b1175d-1.csv") %>% 
  filter(sex == "Both" & measure == "Incidence") %>% 
  merge(y = read.csv("Data/Cleaned/age_lookup.csv")) %>% 
  filter(!is.na(age_mid)) %>% 
  mutate(cause = dplyr::recode(cause,
                               "Breast cancer" = "inc_BC",
                               "Ischemic heart disease" = "inc_IHD",
                               "Ischemic stroke" = "inc_IS", 
                               "Diabetes mellitus type 2" = "inc_T2D", 
                               "Colon and rectum cancer" = "inc_CC"),
         val = val/100000,
         upper = upper/100000,
         lower = lower/100000,
         source = "GBD") %>% 
  mutate(name = cause, 
         age = age_mid,
         mean = val) %>% 
  select(age, name, upper, lower, mean, source)

# GBD prevalence updated
df_prev <- read.csv(file = "Data/Cleaned/IHME-GBD_2019_DATA-c4b1175d-1.csv") %>% 
  filter(sex == "Both" & measure == "Prevalence") %>% 
  merge(y = read.csv("Data/Cleaned/age_lookup.csv")) %>% 
  filter(!is.na(age_mid)) %>% 
  mutate(cause = dplyr::recode(cause,
         "Breast cancer" = "prev_BC",
         "Ischemic heart disease" = "prev_IHD",
         "Ischemic stroke" = "prev_IS", 
         "Diabetes mellitus type 2" = "prev_T2D", 
         "Colon and rectum cancer" = "prev_CC"),
         val = val/100000,
         upper = upper/100000,
         lower = lower/100000,
         source = "GBD") %>% 
  mutate(name = cause, 
         age = age_mid,
         mean = val) %>% 
  select(age, name, upper, lower, mean, source)

## Clean simulation data ----

# extract mean and intervals around...
mean_data <- apply(
  X = a_PSA_results,
  MARGIN =  c(1, 2),
  FUN = mean,
  na.rm = T
)

# convert to long format....
long_mean <- tidyr::pivot_longer(data = mean_data %>% 
                                   as.data.frame %>% 
                                   mutate(age = 1:80), 
                                  cols = -age)

lower_bound <- apply(
  X = a_PSA_results,
  MARGIN =  c(1, 2),
  FUN = quantile,
  na.rm = T,
  probs = 0.025
)

# convert to long format
long_lower <- tidyr::pivot_longer(data = lower_bound %>% 
                                    as.data.frame %>% 
                                    mutate(age = 1:80), 
                                  cols = -age)

upper_bound <- apply(
  X = a_PSA_results,
  MARGIN =  c(1, 2),
  FUN = quantile,
  na.rm = T,
  probs = 0.975
)

long_upper <- tidyr::pivot_longer(data = upper_bound %>% as.data.frame %>% mutate(age = 1:80), 
                    cols = -age)

# merge three datasets
m_sim_data_merged <- merge(x = long_upper, y = long_lower, by = c("age", "name"))
m_sim_data_merged <- merge(x = m_sim_data_merged, y = long_mean, by = c("age", "name"))
colnames(m_sim_data_merged) <- c("age", "name", "upper", "lower", "mean")

# set the conditions vectors & insert into labels for the plot.
v_conds <- c("BC", "CC", "T2D", "IS", "IHD")
v_labels <-  c(paste0("Prevalence ", v_conds), paste0("Incidence ", v_conds))
names(v_labels)  <-  c(paste0("prev_", v_conds), paste0("inc_", v_conds))

# limit data to metrics of interest interest and prevalence...
df_sim <- m_sim_data_merged %>% 
  filter(stringr::str_detect(name, "inc|prev") & age > 11)
df_sim$source <- "simulation"

## Combine datasets ----
df_plot <- rbind(df_sim, df_inc)
df_plot2 <- rbind(df_plot, df_prev)

## Plot Incidence & Prevalence ----
# containing incidence and prevalence vs data used to inform model.
p1 <- ggplot(data = df_plot2)+
  geom_line(mapping =  aes(y = mean, x = age, col = source), alpha = 0.7)+
  geom_ribbon(mapping =  aes(x = age, ymin = lower, ymax = upper, fill = source), 
              alpha = 0.2)+
  scale_y_continuous("Incidence & Prevalence Rate", expand = c(0,0))+
  scale_x_continuous("Age in years", limits = c(11,80), expand = c(0,0))+
  scale_fill_manual(
    "Source of Data",
    values = c("GBD" = "red", 
               "simulation" = "blue"),
    aesthetics = c("fill", "col"),
    labels = c("GBD" = "GBD Study", 
               "simulation" = "Model Results")
  )+
  facet_wrap(~name,
             labeller = labeller(name = v_labels),
             scales = "free", 
             nrow = 2)+
  theme_classic()+
  theme(legend.position = "top")

p1

# output figure ready for publication
ggsave(plot = p1, 
       filename = "Figures/baseline_epi_comparison.png",
       width = 15, 
       height = 8)


# SURVIVAL ------

## load data
tib_lifetable <- readxl::read_excel("Data/Cleaned/ons_lifetables.xls", 
                                    sheet = "2014-2016")

# calculate survival given survival at age 11
v_S_adjusted <- base::rowSums(tib_lifetable[tib_lifetable$age %in% 1:80, c("male_lx", "fmle_lx")]) / 
                                 sum(as.numeric(tib_lifetable[tib_lifetable$age == 11, c("male_lx", "fmle_lx")]))

df_survival_sim <- t(apply(
  X = a_PSA_results[,"v_S",],
  MARGIN =  1,
  FUN = quantile,
  na.rm = T,
  probs = c(0.025,0.975)
)) %>% 
  as.data.frame()

# add in the mean to the dataframe
df_survival_sim$mean   <- rowMeans(a_PSA_results[,"v_S",], na.rm = T)
df_survival_sim$age    <- 1:80
df_survival_sim$source <- "simulation"

colnames(df_survival_sim)[1:2] <- c("lower", "upper")

# create the equivalent for the Life-tables data
df_lifetables <- data.frame(mean = v_S_adjusted,
           lower = NA, 
           upper = NA,
           age = 1:length(v_S_adjusted),
           source = "lifetables")

# create the plotting dataframe
df_plot <- rbind(df_lifetables, df_survival_sim)

# create the ggplot using the data created.
p2 <- ggplot(data = df_plot %>% filter(age > 10),
       mapping = aes(x = age))+
  geom_line(mapping = aes(y = mean, 
                          col = source), 
            alpha = 0.8)+
  geom_ribbon(mapping =  aes(ymin = lower, 
                             ymax = upper, 
                             fill = source), 
              alpha = 0.2)+
  scale_y_continuous(name = "Survival from age 11", 
                     limits = c(0,1),
                     expand = c(0,0))+
  scale_x_continuous(name = "Age in years",
                     limits = c(11, 80))+
  scale_fill_manual(
    "Source of Data",
    values = c("lifetables" = "red", "simulation" = "blue"),
    aesthetics = c("fill", "col"),
    labels = c("lifetables" = "ONS Lifetables", "simulation" = "Model Results")
  )+
  theme_classic()+
  theme(legend.position = "top")

p2
  
ggsave(plot = p2, 
       filename = "Figures/survival_comparison.png",
       width = 15, 
       height = 8, 
       scale = 0.5)

# UTILITIES ----

source("R/econ_functions.R")
## load AB data ----
# baseline utiltiies.
df_AB <- read.csv(file = "Data/Cleaned/baseline_utilities.csv")
# vector containing the proportion of the decrement of utility for those with 
# one or more health conditions which is associated with included conditions. 
# 1 would be that the conditions included account for 100% of the burden at age
# 0 would be that the conditions account for 0% of the burden at age
# baseline assumption is that almost none of burden at age 11 is due to these
# conditions which generally occur later in life, but this increases to 50%
# over the lifecourse.
# https://ourworldindata.org/grapher/share-of-total-disease-burden-by-cause?country=~GBR
v_hs_burden <- seq(from = 0, to = 0.5, length.out = 80)
# use this to calculate a dataframe of utility values with three columns
# general population utility
# healthy population utility - utility in those with no reported health conditions
# no_HS - utility for those without the health conditions included in this analysis.
df_util_AB <- data.frame(
  age = 1:80,
  mean =
    calcPopUtils(
      v_age = df_AB$Age,
      v_U_gp = df_AB$GP_U,
      v_U_H = df_AB$H_U,
      HS_burden = v_hs_burden[1:80],
      age_range = 1:80
    )$gp,
  
  lower = calcPopUtils(
    v_age = df_AB$Age,
    v_U_gp = substr(
      x = df_AB$GP_R,
      start = 2,
      stop = 6
    ),
    v_U_H = df_AB$H_U,
    HS_burden = v_hs_burden[1:80],
    age_range = 1:80
  )$gp,
  upper =
    
    calcPopUtils(
      v_age = df_AB$Age,
      v_U_gp = substr(
        x = df_AB$GP_R,
        start = 8,
        stop = 12
      ),
      v_U_H = df_AB$H_U,
      HS_burden = v_hs_burden[1:80],
      age_range = 1:80
    )$gp,
  source = "AB"
)

## extract simulation data ----
# extract quantiles from the utility estimates in the simulation.
df_utility_sim <- t(apply(
  X = a_PSA_results[,"util",],
  MARGIN =  1,
  FUN = quantile,
  na.rm = T,
  probs = c(0.025,0.975)
)) %>% 
  as.data.frame()

# add in the mean to the dataframe
df_utility_sim$mean   <- rowMeans(a_PSA_results[,"util",], na.rm = T)
df_utility_sim$age    <- 1:80
df_utility_sim$source <- "simulation"
colnames(df_utility_sim)[1:2] <- c("lower", "upper")


# combine the two inputs
df_plot <- rbind(df_utility_sim, df_util_AB)

## plot utilities ----
p3 <- ggplot(data = df_plot %>% filter(age > 10 & age < 80),
             mapping = aes(x = age))+
  geom_line(mapping = aes(y = mean, col = source), alpha = 0.8)+
  geom_ribbon(mapping =  aes(ymin = lower, ymax = upper, fill = source), alpha = 0.2)+
  scale_y_continuous("Mean Utility", limits = c(0, 1), expand = c(0, 0))+
  scale_x_continuous("Age in years", limits = c(11, 80))+
  scale_fill_manual(
    "Source of Data",
    values = c("AB" = "red", "simulation" = "blue"),
    aesthetics = c("fill", "colour"),
    labels = c("AB" = "Ara & Brazier (2011)", "simulation" = "Model Results")
  )+
  theme_classic()+
  theme(legend.position = "top")

p3

# save plot of utilities
ggsave(plot = p3, 
       filename = "Figures/utilities_comparison.png",
       width = 15, 
       height = 8,
       scale = 0.5)


# ECON FIGS ----

rm(list = ls())

lambda <- 20000
lambda_range <- c("min" = 0, "max" = 100000)

library(darkpeak)
source("R/econ_functions.R")

l_PSA_results <- readRDS("Results/base.rds") 


# set of PSA results.
df_psa_res <- calcPSAoutcomes(l_PSA_results)



# combine into matrices for base plots
m_total_costs = cbind("Baseline" = df_psa_res$v_C_base,
                      "Intervention" = df_psa_res$v_C_int + l_PSA_results$params["cost_intervention"])
m_total_qalys = cbind("Baseline" = df_psa_res$v_U_base,
                      "Intervention" = df_psa_res$v_U_int)




## Create ICER Table ----
tab_ICER <- darkpeak::createICERtable(total_costs = m_total_costs,
                                      total_qalys = m_total_qalys,
                                      ref_index = 1)
tab_ICER

## Create CEPlane ----
plot_CEPlane <- darkpeak::makeCEPlane(total_costs = m_total_costs,
                      total_qalys = m_total_qalys,
                      treatment = "Intervention", 
                      comparitor = "Baseline", 
                      thresh = lambda,
                      show_ellipse = T)+
  scale_y_continuous(name = "Incremental Costs", 
                     labels = scales::dollar_format(prefix = "£", 
                                                    decimal.mark = ","))

plot_CEPlane

ggsave(plot = plot_CEPlane, 
       filename = "Figures/ceplane.png", 
       scale = 1.2)

## Create CEAC ----
plot_CEAC <- makeCEAC(
  total_costs = m_total_costs,
  total_qalys = m_total_qalys,
  treatment = c("Intervention", "Baseline"),
  lambda_min = lambda_range["min"],
  lambda_max = lambda_range["max"] / 2) +
  ggplot2::scale_color_manual(name = "Treatment",
                              values = c("Baseline" = "black",
                                         "Intervention" = "grey")) +
  ggplot2::scale_linetype_manual(name = "Treatment",
                                 values = c("Baseline" = 1,
                                            "Intervention" = 2))+
  ggplot2::scale_x_continuous(name = "Willingness to Pay",
                              labels = scales::dollar_format(prefix = "£",
                                                            decimal.mark = ",",
                                                            scale = 0.001,
                                                            suffix = "k"),
                              expand = c(0,0.1)
                              )+
  theme(plot.margin = margin(10, 10, 10, 10))

plot_CEAC

ggsave(plot = plot_CEAC, 
       filename = "Figures/ceac.png",
       width = 12, height = 8, scale = .5)


## Check stability ----

darkpeak::checkStability(
  total_costs = m_total_costs,
  total_qalys = m_total_qalys,
  strategies = c("Intervention","Baseline"),
  BS_samples = 1000,
  withinShiny = F,
  lambda = lambda)

# run analysis to assess stability of PSA results
# get inc costs & qalys
df_stability <- 
  data.frame(inc_costs = m_total_costs[, "Intervention"] - m_total_costs[ ,"Baseline"],
             inc_qalys = m_total_qalys[, "Intervention"] - m_total_qalys[ ,"Baseline"]) %>% 
  mutate(cummean_costs = cummean(inc_costs),# calculate cum mean vectors for each
         cummean_qalys = cummean(inc_qalys),
         v_cum_icer    = cummean_costs / cummean_qalys,
         mean_icer     = mean(inc_costs) / mean(inc_qalys),
         iterations    = 1:n())
  # get the ICER over iterations

# make plot to show how stable...
p_stability <- ggplot2::ggplot(data = df_stability %>% filter(iterations > 500),
                mapping = aes(x = iterations, y = v_cum_icer)) +
  theme_minimal() +
  geom_line() +
  geom_hline(yintercept = unique(mean_icer), 
             linetype = "dashed", 
             col = "red") +
  geom_ribbon(data = df_stability,
              aes(ymin = 0.95 * unique(mean_icer), 
                  ymax = 1.05 * unique(mean_icer)),
              alpha = 0.2, 
              fill = "red")+
  scale_y_continuous(name = "Incremental Cost-effectiveness Ratio",
                     labels = scales::dollar_format(prefix = "£",
                                           decimal.mark = ","))+
  scale_x_continuous(name = "Iterations", 
                     limits = c(0, nrow(df_stability)), 
                     expand = c(0,0))+
  theme(plot.margin = margin(10, 20, 10, 10))



ggsave(plot = p_stability, 
       filename = "Figures/stability_icer.png",
       width = 12, height = 8, scale = .5)
  


## EJPA ----
rm(list = ls())
source("R/econ_functions.R")
library(dplyr)
library(ggplot2)

lambda <- 20000
# Run model with three different levels of effectiveness to start with
# required three additional scenarios.
# Then run the below for each level of effectiveness and store the results
# Finally, combine the three EJPA vectors (prob CE) into a single 'long'
# dataframe and plot with three distinct colours.

# vector of feasible costs
getEJPAmat <- function(min_cost = 0,
                       max_cost = 1000,
                       n_i = l_PSA_results$params["n_i"],
                       lambda = 20000,
                       v_U_int = df_psa_res$v_U_int,
                       v_U_base = df_psa_res$v_U_base,
                       v_C_int = df_psa_res$v_C_int,
                       v_C_base = df_psa_res$v_C_base) {
  
  # a dataframe containing the intervention cost only
  df_plot <- data.frame(c_int = seq(min_cost, max_cost, length.out = 10000))
  
  # create matrix of NMB at each price in the vector.
  df_plot$pCE_50 <- sapply(
    X = df_plot$c_int,
    FUN = function(v_U_int,
                   v_U_base,
                   v_C_int,
                   v_C_base,
                   lambda,
                   c_int) {
      nmb <- (v_U_int - v_U_base) * lambda - (v_C_int + c_int - v_C_base)
      p_CE <- mean(nmb > 0)
      return(p_CE)
    },
    v_U_int = v_U_int,
    v_U_base = v_U_base,
    v_C_int = v_C_int,
    v_C_base = v_C_base,
    lambda = lambda
  )
  
  return(df_plot)
}

# READ IN DATA FOR THE DEFAULT SCENARIO...
n_i <- readRDS("Results/base.rds")$params["n_i"]
df_psa_res_CENTRAL <- readRDS("Results/base.rds") %>% calcPSAoutcomes
df_psa_res_UPPER <- readRDS("Results/high_effect.rds") %>% calcPSAoutcomes
df_psa_res_LOWER <- readRDS("Results/low_effect.rds") %>% calcPSAoutcomes

# maximum cost for analysis
max_cost <- 75

# effectiveness 1
df_plot_central <- getEJPAmat(
  min_cost = 0,
  max_cost = max_cost,
  n_i = n_i,
  lambda = lambda,
  v_U_int = df_psa_res_CENTRAL$v_U_int,
  v_U_base = df_psa_res_CENTRAL$v_U_base,
  v_C_int = df_psa_res_CENTRAL$v_C_int,
  v_C_base = df_psa_res_CENTRAL$v_C_base
) %>%
  mutate(scenario = "central")

# effectiveness 2 estimates
df_plot_upper <- getEJPAmat(
  min_cost = 0,
  max_cost = max_cost,
  n_i = n_i,
  lambda = lambda,
  v_U_int = df_psa_res_UPPER$v_U_int,
  v_U_base = df_psa_res_UPPER$v_U_base,
  v_C_int = df_psa_res_UPPER$v_C_int,
  v_C_base = df_psa_res_UPPER$v_C_base
) %>%
  mutate(scenario = "upper")

# effectiveness 3 estimates
df_plot_lower <- getEJPAmat(
  min_cost = 0,
  max_cost = max_cost,
  n_i = n_i,
  lambda = lambda,
  v_U_int = df_psa_res_LOWER$v_U_int,
  v_U_base = df_psa_res_LOWER$v_U_base,
  v_C_int = df_psa_res_LOWER$v_C_int,
  v_C_base = df_psa_res_LOWER$v_C_base
) %>%
  mutate(scenario = "lower")

# create long dataframe containing all lines...
df_plot <- rbind(df_plot_central, 
                 df_plot_lower, 
                 df_plot_upper)

# then calculate the price at which there is a 50% prob each is CE
df_line_cross <- df_plot %>% 
  group_by(scenario) %>%
  summarise("50perc" = c_int[which.min(abs(pCE_50 - 0.5))])


# create plot.
# col = effectiveness scenario
pEJPA <- ggplot()+
  theme_minimal()+
  geom_line(data    = df_plot, 
             mapping = aes(x = c_int, y = pCE_50, col = scenario),
             alpha = 1,
             size = 1.25)+
  geom_vline(data = df_line_cross,
             aes(xintercept  = `50perc`, col = scenario),
             linetype = 3)+
  theme_classic()+
  geom_hline(yintercept = 0.5, linetype = 2)+
  scale_y_continuous(name = "Probability Cost-effective", 
                     expand = c(0,0),
                     limits = c(0,1))+
  scale_x_continuous(name = "Cost of Intervention", 
                     expand = c(0,0), 
                     breaks = seq(0,75,15),
                     limits = c(0, max_cost),
                     labels = scales::dollar_format(prefix = "£",
                                                    decimal.mark = ","))+
  scale_colour_manual(name = "Intervention Effectiveness",
                     values = c("central" = "#FF7E33", 
                                "upper" = "#321CAA",
                                "lower" = "#15d438"),
                     labels = c("central" = "186 METs (Default)", 
                                "upper" = "279 METs",
                                "lower" = "93 METs"))+
  theme(legend.position = "top",
        plot.margin = margin(10,10,10,10))

pEJPA

ggsave(plot = pEJPA, 
       filename = "Figures/Final/evpa.png",
       width = 12, 
       height = 8, 
       scale = 0.7)

# SENSITIVITY ANALYSIS ----

rm(list = ls())

# economic functions loaded...
source("R/econ_functions.R")
# lambda 
lambda <- 20000

# read in scenario information for what has been run
df_sensScenarios <- read.csv("Results/df_sensScenarios.csv", 
                             row.names = 1)

# read in data from each scenario one at a time...
# here we read in all of them, if we want only 1 then:
# use df_sensScenarios$run in the row...
df_results_scenario <- lapply(
  X = df_sensScenarios[, "name"],
  FUN = function(x) {
    temp <- readRDS(paste0("Results/", x, ".rds")) 
    
    temp <- temp %>%
      calcPSAoutcomes %>%
      as.data.frame %>%
      mutate(scenario = paste(x),
             #n_psa = temp$params["n_psa"],
             cost_intervention = temp$params["cost_intervention"]
             )
    return(temp)
  }
) %>% 
  do.call(what = "rbind") %>% 
  # Calculate Incremental Costs and QALYs in each scenario...
  mutate(
    v_C_int = v_C_int + cost_intervention,
    inc_costs = v_C_int - v_C_base,
    inc_qalys = v_U_int - v_U_base
  )

df_results_scenario$scenario <- factor(x = df_results_scenario$scenario,
                                       levels = c("base", "decay_lin", 
                                                  "utils_minimum", "dr_1_5", 
                                                  "habit_lower", "habit_upper", 
                                                  "low_effect", "high_effect"))

# calc means
df_means <- df_results_scenario %>% 
  group_by(scenario) %>% 
  summarise(inc_costs = mean(inc_costs),
            inc_qalys = mean(inc_qalys))
  

# create Plane
v_plane_cols <- RColorBrewer::brewer.pal(n = length(df_sensScenarios[, "name"]),
                                         name = "Set1")
names(v_plane_cols) <- df_sensScenarios[, "name"]
v_plane_labs <- df_sensScenarios$label
#v_plane_labs <- factor(x = df_sensScenarios$label,
#                       ordered = T,
#                       labels = c("Baseline", "Linear decay", 
#                                  "Minimum HRQoL method", "1.5% Discount Rate", 
#                                  "0% develop habit", "100% develop habit", 
#                                  "279 METs", "93 METs"),
#                       levels = c("Baseline", "Linear decay", 
#                                  "Minimum HRQoL method", "1.5% Discount Rate", 
#                                  "0% develop habit", "100% develop habit", 
#                                  "279 METs", "93 METs"))
names(v_plane_labs) <- df_sensScenarios$name

# ggplot function to build main plot.
plot_cePlane <- 
  ggplot2::ggplot(mapping = ggplot2::aes(x = inc_qalys, y = inc_costs, col = scenario))  +
  # all points
  ggplot2::geom_point(data = df_results_scenario,
                      alpha = 0.2, 
                      shape = 20, 
                      size = 0.7) +  
  # means
  ggplot2::geom_point(data = df_means,
                      alpha = 1, 
                      shape = 17, 
                      size = 3, 
                      col = "black") +
  
  ggplot2::stat_ellipse(data = df_results_scenario,
                        type = "norm",
                        level = 0.5,
                        segments = 50,
                        size = 1) +
  
  ggplot2::geom_abline(slope = lambda, 
                       linetype = "dashed", intercept = 0)+
  
  # coordinate systems
  ggplot2::geom_vline(xintercept = 0, col = "grey") +
  ggplot2::geom_hline(yintercept = 0, col = "grey") +
  
  ggplot2::theme_minimal()  +

  ggplot2::scale_y_continuous(name = "Incremental Costs",
                     labels = scales::dollar_format(prefix = "£",
                                                    decimal.mark = ","),
                     limits = c(-300, 200) ) +
  ggplot2::scale_x_continuous(name = "Incremental Quality Adjusted Life Years",
                              limits = c(-0.05, 0.05))+
  ggplot2::scale_color_manual(name = "Scenario",
                              guide = "none",
                              labels = v_plane_labs,
                              values = v_plane_cols,
                              aesthetics = c("fill", "col"))+
  ggplot2::facet_wrap(~ scenario, 
                      labeller = ggplot2::as_labeller(v_plane_labs),
                      ncol = 2,
                      dir = "h")

  plot_cePlane

# save the cost-effectiveness plane created using ggplot2.
ggplot2::ggsave(plot = plot_cePlane, 
                filename = "Figures/sensitivity_CEPlane.png",
                width = 15,
                height = 20, 
                scale = 0.5)



# RESULTS TABLE
# create a results table that contains the costs and QALYs for each sensitivity 
# analysis
# store mean costs and qalys from the model output
c_f = function(x, d=2, ci=T){
  x = c(mean(x),quantile(x, probs = c(0.025,0.975)))
  x = formatC(x, digits=d, big.mark = ",",format = "f", drop0trailing = F)
  if(ci){
    x = paste0(x[1]," (",x[2],"; ",x[3],")")
  } else {
    x = paste0(x[1])
  }
  x[grepl("0 \\(", substr(x,1,3))] = "- (-, -)"
  x[grepl("0[.]0 \\(", substr(x,1,5))] = "- (-, -)"
  x[grepl("0[.]00 \\(", substr(x,1,6))] = "- (-, -)"
  x[grepl("0[.]000 \\(", substr(x,1,7))] = "- (-, -)"
  return(x)
}


tib_results_tab <- 
  df_results_scenario %>%
  group_by(scenario) %>%
  #summarise(format(round(mean(v_U_int, na.rm = T), 3), nsmall = 3))
  summarise(
    ICER = c_f(mean(inc_costs) / mean(inc_qalys), d = 2, ci = F),
    NMB  = c_f(inc_qalys * lambda - inc_costs, d = 0),
    C_noInt = c_f(v_C_base),
    C_Scen = c_f(v_C_int),
    U_noInt = c_f(v_U_base, d = 3),
    U_Scen = c_f(v_U_int, d = 3),
    inc_costs = c_f(inc_costs, d = 3),
    inc_qalys = c_f(inc_qalys, d = 3)
    )
# pivot longer and then rearrange to a data-frame of results in each scenario...
# COSTS
tib1 <- tib_results_tab %>% 
  tidyr::pivot_longer(cols = c("C_noInt", "C_Scen"), names_prefix = "C_", values_to = "Total Cost") %>%
  mutate(scenario_name = paste0(scenario, "_", name)) 
# UTILITIES
tib2 <- tib_results_tab %>% 
  tidyr::pivot_longer(cols = c("U_noInt", "U_Scen"), names_prefix = "U_", values_to = "Total Utility") %>%
  mutate(scenario_name = paste0(scenario, "_", name))

# TIDY
df_ICER_table <- merge(tib1[, c("scenario", "scenario_name", "Total Cost", "inc_costs", "name")], 
      tib2[, c("scenario", "scenario_name", "Total Utility", "inc_qalys", "ICER", "NMB")]) %>% 
  mutate(inc_costs = ifelse(test = name == "noInt", NA, inc_costs), 
         inc_qalys = ifelse(test = name == "noInt", NA, inc_qalys),   
         ICER = ifelse(test = name == "noInt", NA, ICER)) %>% 
  select(scenario, `Total Cost`, `Total Utility`, `inc_costs`, `inc_qalys`, ICER, NMB) %>% 
  rename(`Incremental Costs` = inc_costs,
         `Incremental Utility` = inc_qalys)

# save to csv
write.csv(x = df_ICER_table,
          file =  "Results/sensitivity_ICER_table.csv")



# TORNADO DIAGRAM

df_temp <- df_results_scenario %>%
  group_by(scenario) %>%
  #summarise(format(round(mean(v_U_int, na.rm = T), 3), nsmall = 3))
  summarise(
    NMB = mean(inc_qalys) * lambda - mean(inc_costs),
    ICER = mean(inc_costs) / mean(inc_qalys),
    C_noInt = mean(v_C_base),
    C_Scen = mean(v_C_int),
    U_noInt = mean(v_U_base),
    U_Scen = mean(v_U_int),
    inc_costs = mean(inc_costs),
    inc_qalys = mean(inc_qalys)
  ) %>% 
  mutate(icer_base = ICER[1],
         icer_sens = ICER,
         icer_lb = ifelse(test = icer_sens <= icer_base, yes =  icer_sens, no =  NA),
         icer_ub = ifelse(test = icer_sens >  icer_base, yes =  icer_sens, no =  NA),
         nmb_base = NMB[1],
         nmb_sens = NMB,
         nmb_lb = ifelse(test = nmb_sens <= nmb_base, yes =  nmb_sens, no =  NA),
         nmb_ub = ifelse(test = nmb_sens >  nmb_base, yes =  nmb_sens, no =  NA)) 

df_plot <- df_temp %>% 
  mutate(diff = abs(nmb_sens - nmb_base)) %>%
  arrange(diff) %>%
  mutate(rank = 1:n())
  
merge(df_plot, df_sensScenarios[,c("name", "label")])

plot_Tornado <- ggplot(df_plot, 
       aes(x = nmb_base, y = rank)) +
  geom_segment(
    aes(
      xend = nmb_lb,
      yend = rank,
    ),
    size = 7,
    lineend = "butt",
    col = "grey"
  ) +
  geom_segment(
    aes(
      xend = nmb_ub,
      yend = rank,
    ),
    size = 7,
    lineend = "butt",
    col = "grey"
  )+
  scale_y_discrete(name = "Scenarios",
                   limits = df_plot$scenario,
                   labels = v_plane_labs,
                   position = "left") +
  scale_x_continuous(name = "Net Monetary Benefit",
                     labels = scales::dollar_format(prefix = "£",
                                                    decimal.mark = ","),
                     limits = c(df_plot$nmb_base[1] - max(df_plot$diff)*1.05, df_plot$nmb_base[1] + max(df_plot$diff*1.05)))+
  
  theme_minimal()+
  theme(plot.margin = margin(10,10,10,10))+
  geom_vline(xintercept = df_plot$nmb_base[1])

ggsave(plot = plot_Tornado, 
       filename = "Figures/sensitivity_Tornado.png",
       scale = 0.7, 
       width = 8, 
       height = 12)
