
# FUNCTIONS ----

#' Create Physical Activity Percentile Distribution
#'
#' This function estimates the physical activity percentiles for males and females based on a baseline simulation.
#' @param v_mets A vector of physical activity levels, measured in mets.
#' @param v_female A vector indicating sex, with 1 for female and 0 for male.
#' @return A data frame with the physical activity percentiles for females and males.
#' @examples
#' createPAdist(v_mets = c(1.5, 2.5, 3.5), v_female = c(1, 0, 1))
createPAdist <- function(v_mets,
                         v_female) {
  # here we want percentiles (total of 100 from 0.01 to 1)
  quantiles <- seq(0.01, 1, 0.01)
  # use the baseline simulation to estimate the physical activity percentiles for males and females.
  df_mets_perc_18 <-  data.frame(
    fmle = quantile(x = v_mets[v_female == 1], quantiles),
    male = quantile(x = v_mets[v_female == 0], quantiles)
  )
  return(df_mets_perc_18)
}

#' Assign Physical Activity Percentile
#'
#' This function assigns a vector of individuals their percentile of physical activity within their sex.
#' @param v_mets A vector of physical activity levels, measured in mets.
#' @param v_female A vector indicating sex, with 1 for female and 0 for male.
#' @param df_PAdist A data frame with the physical activity percentiles for males and females.
#' @return A vector of physical activity percentiles for the individuals.
#' @examples
#' df_PAdist <- createPAdist(v_mets = c(1.5, 2.5, 3.5), v_female = c(1, 0, 1))
#' getActivityPercentile(v_mets = c(1.5, 2.5, 3.5), v_female = c(1, 0, 1), df_PAdist = df_PAdist)
#' getActivityPercentile(v_mets = c(1.5, 2.5), v_female = c(1, 0), df_PAdist = df_PAdist)
#' getActivityPercentile(v_mets = c(1.5), v_female = c(1), df_PAdist = df_PAdist)
assignPAcentile <- function(v_mets,
                            v_female,
                            df_PAdist){
  
  # create empty matix for percentiles
  v_met_percentiles <- rep(x = NA, length(v_mets))
  
  # female first, assign physical activity percentile closest to level of PA 
  v_met_percentiles[v_female==1] <- base::findInterval(x = v_mets[v_female==1], 
                                                       vec = df_PAdist$fmle,
                                                       all.inside = T) 
  # males second
  v_met_percentiles[v_female==0] <- base::findInterval(x = v_mets[v_female==0], 
                                                       vec = df_PAdist$male,
                                                       all.inside = T) 
  # return the vector of PA percentiles within age and sex.
  return(v_met_percentiles)
}



#' Simulate child physical activity using an assumption of linear physical activity decay
#'
#' @param v_pa_start A vector of initial physical activity levels
#' @param pa_decay A scalar representing the rate of decay of physical activity over time
#' @param v_pa_additional A vector of additional PA - taken from effect distribution...
#' @param age_start An integer representing the starting age of the child
#' @param age_adult An integer representing the age at which the child becomes an adult
#'
#' @return A matrix of simulated physical activity levels over time
#'
#' @examples
#' simChildPA(v_pa_start = c(60, 70, 80),
#'            pa_decay = 0.06,
#'            v_pa_additional = rnorm(n = 3, mean =  186, 10),
#'            age_start = 5,
#'            age_adult = 18) - 
#'   simChildPA(v_pa_start = c(60, 70, 80),
#'              pa_decay = 0.06,
#'              v_pa_additional = 0,
#'              age_start = 5,
#'              age_adult = 18)
#'
simChildPA <- function(v_pa_start,
                       pa_decay,
                       v_pa_additional,
                       age_start,
                       age_adult) {
  
  # convert the vector of initial PA levels to a matrix
  m_pa <- as.matrix(v_pa_start)
  # for each, extrapolate linear decay of physical activity
  # over the period
  m_pa <- apply(
    X = m_pa,
    MARGIN = 1,
    FUN = function(x) {
      x * (1 - pa_decay) ^ (0:(age_adult - age_start))
    }
  )
  # additional physical activity (e.g. for example due to unit increase
  # rather than a decay rate decrease)
  m_pa <- t(t(m_pa)+v_pa_additional)
  # physical activity rates lower than 0 recoded to keep sensible levels
  m_pa[m_pa < 0] <- 10
  # check done correctly, no errors
  assertthat::assert_that(min(m_pa) > 0,
                          msg = "Values of PA lower than zero")
  return(m_pa)
}

#' Estimate adult physical activity (PA) trajectory
#'
#' @param v_female A vector denoting whether each individual is female (1) or male (0).
#' @param v_percentile A vector denoting the assigned percentile at age 18 for each individual.
#' @param age_adult The start age for the adult PA trajectory.
#' @param age_end The end age for the adult PA trajectory.
#' @param df_PA_percentiles A data frame of PA percentiles by age from HSE14-17, percentile on row, age on column (with females starting at 101)
#' @return A matrix with one row for each cycle and one column for each person, containing their estimated PA trajectory over the specified period.
#' @export
assignAdultPA_traj <- function(v_female, 
                               v_percentile, 
                               age_adult, 
                               age_end, 
                               df_PA_percentiles){
  # Check that the sex and percentile vectors are the same length
  assertthat::are_equal(length(v_female), length(v_percentile))
  
  # Create an empty matrix with one row for each cycle and one column for each person
  m_adult_PA_traj <- matrix(
    data = NA,
    nrow = age_end - age_adult + 1,
    ncol = length(v_female)
  )
  
  # Loop through each person
  for(x in 1:ncol(m_adult_PA_traj)){
    # Physical activity trajectory for individual given sex and age
    m_adult_PA_traj[, x] <- as.numeric(x = df_PA_percentiles[ v_percentile[x], v_female[x] * 100 + age_adult:age_end ])
  }
  
  # Update any zero METmins to be 1 (for later calculations)
  m_adult_PA_traj[m_adult_PA_traj[] == 0] <- 1
  
  # Check that no values above 10000 or below zero
  assertthat::assert_that(all(m_adult_PA_traj <= 10000 & m_adult_PA_traj > 0))
  # Return the matrix of physical activity trajectories over the period for each person
  return(m_adult_PA_traj)
  
}

#' Decay intervention effect
#' 
#' Calculates the decay of an intervention effect on physical activity (PA) trajectory.
#' 
#' @param baseline_PA_trajectory numeric matrix of baseline PA trajectory.
#' @param intervention_PA_trajectory numeric matrix of PA trajectory after intervention.
#' @param effect_decay_type character string indicating the type of decay to use, either "exp" or "lin" (for exponential and linear).
#' @param effect_decay_exp_rate numeric value indicating the annual rate of exponential decay (if chosen).
#' @param effect_decay_linear_duration integer value indicating the number of periods over which to apply linear decay (if chosen).
#' @param effect_decay_proportion numeric value between 0 and 1 indicating the proportion of the population's intervention effect to decay.
#' @return decayed_PA_trajectory numeric matrix of PA trajectory after applying the decay of the intervention effect.
#' @examples
#' decayInteff(baseline_PA_trajectory = baseline_PA,
#'            intervention_PA_trajectory = intervention_PA,
#'            effect_decay_type = "exponential",
#'            effect_decay_exp_rate = 0.1,
#'            effect_decay_linear_duration = 5,
#'            effect_decay_proportion = 0.5)
decayInteff <- function(baseline_PA_trajectory,
                        intervention_PA_trajectory,
                        effect_decay_type = NULL,
                        effect_decay_exp_rate = NULL,
                        effect_decay_linear_duration = NULL,
                        effect_decay_proportion = NULL){
  
  # create the decayed PA trajectory which is the intervention one to start with
  decayed_PA_trajectory <- intervention_PA_trajectory
  # assign the decay multiplier for the remainder of the periods
  # based on the selection of decay type and then duration/annual rate
  if (effect_decay_type == "exp") {
    assertthat::assert_that(!is.null(x = effect_decay_exp_rate), msg = "Exp chosen, no rate provided")
    v_decay_multiplier <-
      (1 - effect_decay_exp_rate) ^ (0:(age_end - age_adult))
  } else if (effect_decay_type == "lin") {
    assertthat::assert_that(!is.null(x = effect_decay_linear_duration), msg = "Linear chosen, no duration provided")
    v_decay_multiplier <-
      c(seq(
        from = 1,
        to = 0,
        length.out = effect_decay_linear_duration + 1
      ),
      rep(0, (
        age_end - age_adult - effect_decay_linear_duration
      )))
  } else{
    message("Error - decay type does not match options")
  }
  
  # EFFECT DECAY PROPORTIONS 
  # If the proportion is stated then less than 100% may be assigned 
  # to have their effect decayed.
  if (is.null(x = effect_decay_proportion) | effect_decay_proportion == 1) {
    # assign into decayed PA trajectory matrix
    decayed_PA_trajectory[age_adult:age_end,] <-
      intervention_PA_trajectory[age_adult:age_end,] * v_decay_multiplier + 
      baseline_PA_trajectory[age_adult:age_end,] * (1 - v_decay_multiplier)
    return(decayed_PA_trajectory)
  } else if(effect_decay_proportion == 0) {
      decayed_PA_trajectory[age_adult:age_end,] <- intervention_PA_trajectory[age_adult:age_end,]
      return(decayed_PA_trajectory)
    } else{
    # basic checks, columns same size on trajectories.
    n_pop <- ncol(intervention_PA_trajectory)
    assertthat::are_equal(n_pop, ncol(baseline_PA_trajectory))
    
    # create individuals to change trajectories for:
    ind_decayed <- base::sample(x = 1:n_pop, 
                                replace = F, 
                                size = floor(n_pop * effect_decay_proportion))
    
    # decayed PA trajectory only for those with decay (some proportion of pop)
    decayed_PA_trajectory[age_adult:age_end, ind_decayed] <-
      intervention_PA_trajectory[age_adult:age_end,  ind_decayed] * v_decay_multiplier + 
      baseline_PA_trajectory[age_adult:age_end, ind_decayed] * (1 - v_decay_multiplier)
  return(decayed_PA_trajectory)
    }
}

#' Generates physical activity (PA) trajectories for individuals
#'
#' This function generates PA trajectories for individuals over time,
#' both for the baseline scenario and an intervention scenario. PA
#' is simulated for children until they reach adulthood, and then
#' assigned based on PA percentiles from age 18 and beyond. The
#' effect of the intervention on PA can be specified as either
#' exponentially or linearly decaying.
#'
#' @param n_i \code{numeric} number of individuals in the model.
#' @param v_init_mets \code{numeric vector} containing initial PA level
#' (in METs) at age 11.
#' @param pa_decay_baseline \code{numeric} annual physical activity
#' decay rate (%) for the baseline scenario.
#' @param pa_decay_intervention \code{numeric} annual physical activity
#' decay rate (%) for the intervention scenario.
#' @param pa_additional_intervention \code{numeric} vector of increased intervention
#' group physical activity trajectory in childhood - additive every year - taken from distribution.
#' @param df_PA_percentiles \code{data frame} of PA percentiles by
#' age from HSE14-17, with percentile on the rows and age on the
#' columns (with females starting at 101).
#' @param age_start \code{integer} age to start the PA trajectory.
#' @param age_adult \code{integer} age at which individuals are
#' considered adults.
#' @param age_end \code{integer} age to end the PA trajectory.
#' @param v_female \code{integer vector} with 1 for female and 0 for
#' male.
#' @param effect_decay_type \code{character} either "exp" or "lin"
#' indicating the type of decay to apply to the effect of the
#' intervention.
#' @param effect_decay_exp_rate \code{numeric} rate of decay for the
#' exponential decay of the intervention effect (between 0 and 1).
#' @param effect_decay_linear_duration \code{integer} duration (in
#' number of years) to apply linear decay of the intervention effect.
#' @param effect_decay_proportion \code{numeric} proportion (between
#' 0 and 1) of the PA percentile to apply decay to.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{m_PA_baseline}: matrix of PA (in METs) for the
#' baseline scenario.
#' \item \code{m_PA_intervention}: matrix of PA (in METs) for the
#' intervention scenario.
#' }
#'
runPAmodule <- function(n_i, # number of individuals in the model
                        v_init_mets, # numeric vector containing initial mets age 11
                        pa_decay_baseline, # numeric value Annual physical activity decay (%) - Baseline
                        pa_decay_intervention, # numeric value Annual physical activity decay (%) - Intervention
                        pa_additional_intervention, # numeric value giving value to add to intervention trajectory compared to baseline
                        df_PA_percentiles, # A data frame of PA percentiles by age from HSE14-17, percentile on row, age on column (with females starting at 101)
                        age_start, # integer value age to start the PA trajectory 
                        age_adult, # integer value age switch from child to adult 
                        age_end,   # integer value age to end the PA trajectory
                        v_female,  # integer vector 1 for female 0 for male
                        effect_decay_type, # character value either exp or lin
                        effect_decay_exp_rate, # numeric value between 0 - 1 denoting rate of decay
                        effect_decay_linear_duration, # integer value duration to apply linear decay of effect
                        effect_decay_proportion # numeric value between 0 -1 proportion to apply decay to
) {
  
  # Physical Activity Module matrix
  m_PA_baseline <- matrix(data = NA, 
                          ncol = n_i,
                          nrow = age_end)
  
  # simulate child physical activity until adulthood, baseline
  m_PA_baseline[age_start:age_adult,] <-
    simChildPA(
      v_pa_start = v_init_mets,
      pa_decay = pa_decay_baseline,
      v_pa_additional = 0,
      age_start = age_start,
      age_adult = age_adult
    )
  
  # simulate child physical activity until adulthood, intervention
  m_PA_intervention <- m_PA_baseline * NA
  
  m_PA_intervention[age_start:age_adult,] <-
    simChildPA(
      v_pa_start = v_init_mets,
      pa_decay = pa_decay_intervention,
      v_pa_additional = pa_additional_intervention,
      age_start = age_start,
      age_adult = age_adult
    )
  
  # Store PA at 18 in the baseline model
  df_PAdist <-  createPAdist(v_mets = m_PA_baseline[age_adult, ],
                             v_female = v_female)
  
  # Assign the percentile of physical activity at age 18,
  # store this as a fixed characteristic, i.e. one that follows the individual
  PAperc_baseline <-
    assignPAcentile(v_mets = m_PA_baseline[age_adult, ],
                    v_female = v_female,
                    df_PAdist = df_PAdist)
  
  # store this as a fixed characteristic, i.e. one that follows the individual
  PAperc_intervention <-
    assignPAcentile(v_mets = m_PA_intervention[age_adult, ],
                    v_female = v_female,
                    df_PAdist = df_PAdist)
  
  
  # Assign in the physical activity trajectories for adulthood - BASELINE
  m_PA_baseline[age_adult:age_end,] <- assignAdultPA_traj(
    v_female = v_female,
    v_percentile = PAperc_baseline,
    df_PA_percentiles = df_PA_percentiles,
    age_adult = age_adult,
    age_end = age_end
  )
  
  # Assign in the physical activity trajectories for adulthood - INTERVENTION
  # create new matrix just for the intervention group
  m_PA_intervention <- m_PA_baseline
  m_PA_intervention[age_adult:age_end,] <- NA
  # fill it in based on percentiles...
  m_PA_intervention[age_adult:age_end,] <- assignAdultPA_traj(
    v_female = v_female,
    v_percentile = PAperc_intervention,
    df_PA_percentiles = df_PA_percentiles,
    age_adult = age_adult,
    age_end = age_end
  )
  
  
  
  # if the decay is selected then we can include a decay into the model.
  m_PA_decayed <- decayInteff(
    baseline_PA_trajectory = m_PA_baseline,
    intervention_PA_trajectory = m_PA_intervention,
    effect_decay_type = effect_decay_type,
    effect_decay_exp_rate = effect_decay_exp_rate,
    effect_decay_linear_duration = effect_decay_linear_duration,
    effect_decay_proportion = effect_decay_proportion
  )
  
  return(list("Intervention" = m_PA_decayed,
              "Baseline" = m_PA_baseline))
}

