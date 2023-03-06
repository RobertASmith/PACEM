# Combine the different per cycle results into a single matrix to be output into an array for PSA
combineResults <- function(m_D,   # a matrix of death status (1 = dead) with rows`for each cycle and columns for each individual`
                           a_prev , # an array of prevalence, with 1 = condition. Cycle as 1 dim, Individual as 2, Conditions as 3
                           a_hs,    # an array of health states (incident) with cycle as 1st dim, individual as 2nd dim, condition as 3rd dim
                           m_results, # matrix containing utilities and costs by cycle 
                           m_PA){  # matrix with physical activity, individuals as columns and cycles as rows
  # check that we have the same number of people and cycles
  assertthat::are_equal(nrow(m_D), nrow(m_results))
  assertthat::are_equal(dim(m_PA), dim(a_hs[,,1]))
  assertthat::are_equal(dim(m_D), dim(m_PA))
  # number of individuals are given by m_D rows
  n_i <- ncol(m_D)
  # calculate the proportion survival for each cycle 
  v_S    <- 1 - base::rowMeans(x = m_D, na.rm = 1)
  # Get a matrix of prevalence at each cycle (conditions as columns, cycles as rows)
  m_prev <- apply(
    X = a_prev,
    MARGIN = 3,
    FUN = function(m, m_D) {
      m[m_D == 1] <- NA
      rowMeans(m, na.rm = T)
    },
    m_D = m_D
  )
  # Get a matrix of incidence at each cycle (conditions as columns, cycles as rows)
  m_inc  <- apply(
    X = a_hs,
    MARGIN = 3,
    FUN = function(m, m_D) {
      m[m_D == 1] <- NA
      rowMeans(m, na.rm = T)
    },
    m_D = m_D
  )
  # get a matrix of physical actiivty percentiles, needs to be transposed.
  m_PA <- t(apply(X = m_PA[1:age_end, ], 
                  MARGIN = 1, 
                  FUN = quantile, 
                  na.rm = T))
  
  # calculate mean utility for those alive
  m_results[, "util"]  <- m_results[, "util"] / n_i / v_S
  m_results[, "costs"] <- m_results[, "costs"] / n_i 
  
  # combine these different metrics into a single matrix.
  m_out <- cbind(m_PA, 
                 v_S, 
                 m_prev, 
                 m_inc,
                 m_results)
  colnames(m_out) <- c(colnames(m_PA), 
                       "v_S", 
                       paste0("prev_", colnames(m_prev)), 
                       paste0("inc_", colnames(m_inc)),
                       colnames(m_results))
  return(m_out)
}


# run the deterministic model, 
runDeterm <- function(output,
                      n_i,
                      v_init_mets,
                      pa_decay_baseline,
                      pa_decay_intervention,
                      pa_additional_intervention,
                      df_PA_percentiles,
                      age_start,
                      age_adult,
                      age_end,
                      v_female,
                      effect_decay_type,
                      effect_decay_exp_rate,
                      effect_decay_linear_duration,
                      effect_decay_proportion,
                      v_constant,
                      v_ln,
                      df_inc,
                      v_durations,
                      df_cond_mort,
                      utility_method,
                      v_u_hp,
                      v_first_year_costs,
                      v_subsequent_year_costs,
                      v_HS_mult,
                      dr_costs,
                      dr_qalys,
                      HSseed){
  
  # Run the physical activity module
  l_PA_module_results <- runPAmodule(
    n_i = n_i,
    v_init_mets = v_init_mets,
    pa_decay_baseline = pa_decay_baseline,
    pa_decay_intervention = pa_decay_intervention,
    pa_additional_intervention = pa_additional_intervention,
    df_PA_percentiles = df_PA_percentiles,
    age_start = age_start,
    age_adult = age_adult,
    age_end = age_end,
    v_female = v_female,
    effect_decay_type = effect_decay_type,
    effect_decay_exp_rate = effect_decay_exp_rate,
    effect_decay_linear_duration = effect_decay_linear_duration,
    effect_decay_proportion = effect_decay_proportion
  )
  
  ## HS module ----
  set.seed(HSseed)
  
  l_EPI_baseline <- runEpiModel(age_start = age_start,
                                age_end = age_end,
                                v_female = v_female,
                                df_PA_traj_baseline = l_PA_module_results[["Baseline"]],
                                df_PA_traj_intervention = l_PA_module_results[["Baseline"]],
                                v_constant = v_constant,
                                v_ln = v_ln,
                                df_inc = df_inc,
                                v_durations = v_durations,
                                df_cond_mort = df_cond_mort)
  set.seed(HSseed)
  
  l_EPI_intervention <- runEpiModel(age_start = age_start,
                                    age_end = age_end,
                                    v_female = v_female,
                                    df_PA_traj_baseline = l_PA_module_results[["Baseline"]],
                                    df_PA_traj_intervention = l_PA_module_results[["Intervention"]],
                                    v_constant = v_constant,
                                    v_ln = v_ln,
                                    df_inc = df_inc,
                                    v_durations = v_durations,
                                    df_cond_mort = df_cond_mort)
  
  ## Econ module ----
  # calculate total utilities and costs for each intervention
  m_results_baseline <- runEconModel(
    v_first_year_costs = v_first_year_costs,
    v_subsequent_year_costs = v_subsequent_year_costs,
    a_hs = l_EPI_baseline[["a_hs"]],
    v_u_hp = v_u_hp,
    utility_method = utility_method,
    m_D = l_EPI_baseline[["m_D"]],
    a_prev = l_EPI_baseline[["a_prev"]],
    v_HS_mult = v_HS_mult,
    dr_costs = dr_costs,
    dr_qalys = dr_qalys)
  
  m_results_intervention <- runEconModel(
    v_first_year_costs = v_first_year_costs,
    v_subsequent_year_costs = v_subsequent_year_costs,
    a_hs = l_EPI_intervention[["a_hs"]],
    v_u_hp = v_u_hp,
    utility_method = utility_method,
    m_D = l_EPI_intervention[["m_D"]],
    a_prev = l_EPI_intervention[["a_prev"]],
    v_HS_mult = v_HS_mult,
    dr_costs = dr_costs,
    dr_qalys = dr_qalys)
  
  if(output == "econ"){
    return(list("base" = m_results_baseline, "int" = m_results_intervention ))
  } else if(output == "epi"){
    
    m_base <- combineResults(m_D = l_EPI_baseline[["m_D"]],
                             a_prev = l_EPI_baseline[["a_prev"]],
                             a_hs = l_EPI_baseline[["a_hs"]],
                             m_results = m_results_baseline,
                             m_PA = l_PA_module_results[["Baseline"]])
    
    m_int <- combineResults(m_D = l_EPI_intervention[["m_D"]],
                            a_prev = l_EPI_intervention[["a_prev"]],
                            a_hs = l_EPI_intervention[["a_hs"]],
                            m_results = m_results_intervention,
                            m_PA = l_PA_module_results[["Intervention"]])
    
    return(list("base" = m_base, "int" = m_int))
 
  }
  
}