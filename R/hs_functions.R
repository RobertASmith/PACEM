# FUNCTIONS ----

#' Calculates individual relative risks for various health conditions
#'
#' This function calculates the relative risks for various health conditions 
#' for individual people, based on their physical activity level (in MET-minutes) 
#' and a log-linear model. It compares the physical activity level of each individual 
#' to a benchmark vector of physical activity levels (from the baseline vector study) 
#' and forces the mean risk in the group to be 1 for each condition. 
#' The function returns a matrix with the relative risks for each individual for each condition.
#'
#' @param v_mets_benchmark \code{vector} of physical activity levels (in MET-minutes) to compare to (e.g. baseline levels).
#' @param v_mets_population \code{vector} of physical activity levels (in MET-minutes) for the population of interest.
#' @param v_constant \code{vector} of constants for the log-linear model.
#' @param v_ln \code{vector} of coefficients for the log-linear model.
#'
#' @return A matrix with the relative risks for each individual for each condition. 
#' The matrix has column names 
#' "BC" (breast cancer), 
#' "CC" (colorectal cancer), 
#' "T2D" (type 2 diabetes), 
#' "IS" (ischemic stroke), 
#' "IHD" (ischemic heart disease).
#' 
f_rr_ind <- function(v_mets_benchmark, # vector of mets which to compare to - may be baseline
                     v_mets_population, # vector of population mets for the scenario
                     v_constant,       # vector of constants for the log linear model
                     v_ln){            # vector of coefficients for the log linear model
  # check lengths of vectors match
  n_i <- length(v_mets_benchmark)
  assertthat::are_equal(n_i, length(v_mets_population))
  n_hs <- length(v_constant)
  assertthat::are_equal(n_hs, length(v_ln))
  
  # create a matrix of the constant for each person
  m_c <- matrix(data = v_constant,
                nrow = n_i, 
                ncol = n_hs, 
                byrow = TRUE)
  
  # create a matrix of the variable for each person give mets
  m_rr_benchmark  <- m_c + log(v_mets_benchmark) %*% t(v_ln)
  m_rr_population <- m_c + log(v_mets_population) %*% t(v_ln)
  
  # create a matrix of the relative risk compared to baseline ...
  # i.e forcing the mean risk in the group to be 1 for each condition
  m_rr_ind <- m_rr_population / matrix(
    data = base::colMeans(m_rr_benchmark),
    nrow = n_i,
    ncol = n_hs,
    byrow = T
  )
  
  # if any of the column means deviate significantly from 1 then warn me
  if(any(abs(base::colMeans(m_rr_ind) - 1) > 0.01) & all(v_mets_population == v_mets_benchmark)){
    message("Warning: relative risk column means deviate significantly from 1,
            this should only occur for the baseline group")
  }
  
  # add in the column names
  colnames(m_rr_ind) <- names(v_constant)
  
  # return a matrix with the relative risk for each condition given met-mins in period.
  return(m_rr_ind)
}

#' @title Calculate cycle-specific risk for a population
#' 
#' @description 
#' This function calculates the cycle-specific risk of certain conditions for a 
#' population of interest, based on their age and sex, their met-mins, 
#' and a log-linear model of PA -> condition using the provided dataframe 
#' containing incidence of condition for each age-sex combination.
#' It first extracts the risk of each condition for each individual using 
#' `f_rr_ind` function and then extracts the risk of each condition for each 
#' individual given age and sex using data-frame df_inc.
#' Finally it returns matrix of risk values for each individual for each condition for the cycle.
#' NOTE: a column must exist for each constant in the log-linear model of PA -> condition risk 
#'
#' @param age numeric - age of the population this cycle
#' @param v_female vector - 1 for female, 0 for male
#' @param v_mets_population vector - met-mins for population of interest
#' @param v_mets_benchmark vector - met-mins for the baseline population
#' @param v_constant constant - in log-linear model of PA -> condition
#' @param v_ln coefficient - in log-linear model of PA -> condition
#' @param df_inc data-frame - containing incidence of condition for each age-sex combination
#'
#' @return matrix - of risk values, for each individual for each condition for the cycle
#'
#' @examples
#' df_example <- data.frame("BC" = c(0.03, 0.05), 
#'                          "CC" = c(0.02, 0.01), 
#'                          "T2D" = c(0.06, 0.02),
#'                          row.names = c("40_1", "40_0")) 
#'calcCycleRisk(
#'  age = 40,
#'  v_female = c(1, 1, 0, 1, 0),
#'  v_mets_population = c(2000, 1500, 1800, 2500, 3000),
#'  v_mets_benchmark = c(10, 10, 10, 10, 10),
#'  v_constant = c("BC" = 1.3, "CC" = 1.1, "T2D" = 1.2),
#'  v_ln = c("BC" = -0.03, "CC" = -0.02, "T2D" = -0.04),
#'  df_inc = df_example)
#' 
calcCycleRisk <- function(
    age, # numeric -  age of the population this cycle
    v_female, # vector - 1 for female, 0 for male
    v_mets_population, # vector of met-mins for population of interest
    v_mets_benchmark,  # vector of met-mins for the baseline population
    v_constant,        # constant in log-linear model of PA -> condition
    v_ln,              # coefficient in log-linear model of PA -> condition
    df_inc){           # data-frame containing incidence of condition each age-sex combination
  
  # extract risk of each condition for each individual
  # individual in row, condition in column...
  m_RR_cycle <- f_rr_ind(
    v_mets_population = v_mets_population,
    v_mets_benchmark = v_mets_benchmark,
    v_constant = v_constant,
    v_ln = v_ln
  )
  
  # extract risk of each condition for each individual given age and sex
  # individual in row, condition in column...
  m_risk <- df_inc[paste0(age, "_", v_female), colnames(m_RR_cycle)]
  
  # return the matrix of risk values, for each individual for each condition for the cycle
  return(m_risk * m_RR_cycle)
}

#' Generate a matrix with identical structure but random numbers from 0 to 1
#'
#' @param m_structure initial matrix on which to copy structure
#' @param lower lower bound for random numbers (default: 0)
#' @param upper upper bound for random numbers (default: 1)
#' @return a matrix with the same dimensions as m_structure, filled with random numbers between lower and upper
#' @examples
#' my_matrix <- matrix(1:9, 3, 3)
#' genRandMatrix(my_matrix)
#' genRandMatrix(my_matrix, lower = -1, upper = 0)
#' 
genRandMatrix <- function(m_structure, 
                          lower = 0, 
                          upper = 1){
  # create a copy
  m_rand <- m_structure
  # fill with random values
  m_rand[] <- runif(n = base::prod(dim(m_structure)),
                    min = lower, 
                    max = upper)
  # return random values
  return(m_rand)
}

#' Generates a matrix of health state for a given age and sex and metmins level
#'
#' @param age age relating to the cycle
#' @param v_female is a binary variable representing sex of the population
#' @param v_mets_population vector numeric amount of metabolic equivalent of task (METs) of physical activity participated by population 
#' @param v_mets_benchmark vector numeric amount of metabolic equivalent of task (METs) of physical activity participated by benchmark group
#' @param v_constant constant variable used in the risk calculation
#' @param v_ln variable in the risk calculation
#' @param df_inc dataframe of incidence of condition by age and sex
#' @return a logical matrix denoting whether an incident case of the condition has occurred in the period for each person
#' @examples
# df_example <- data.frame("BC" = c(0.5, 0.05),
#                        "CC" = c(0.02, 0.01),
#                        "T2D" = c(0.06, 0.02),
#                        row.names = c("40_1", "40_0"))
# getcycleHS(age = 40,
#           v_female = c(1, 1, 0, 1, 0),
#           v_mets_population = c(2000, 1500, 1800, 2500, 3000),
#           v_mets_benchmark = c(10, 10, 10, 10, 10),
#           v_constant = c("BC" = 1.3, "CC" = 1.1, "T2D" = 1.2),
#           v_ln = c("BC" = -0.03, "CC" = -0.02, "T2D" = -0.04),
#           df_inc = df_example)
getcycleHS <- function(age,
                       v_female, 
                       v_mets_population,
                       v_mets_benchmark,
                       v_constant,
                       v_ln,
                       df_inc){
  
  # calculate risk for that age, sex and mets combination
  m_risk <- calcCycleRisk(age = age,
                          v_female = v_female,
                          v_mets_population = v_mets_population,
                          v_mets_benchmark = v_mets_benchmark,
                          v_constant = v_constant,
                          v_ln = v_ln,
                          df_inc = df_inc)
  # assess that risk against a random matrix of values between 0 and 1
  m_hs <- m_risk > genRandMatrix(m_risk)
  
  # return logical matrix denoting whether an incident case of the condition has occurred in the period for each person 
  return(m_hs)
  
}

#' Run incidence model
#' 
#' Simulates the incidence of several conditions over a range of ages for a 
#' group of individuals based on sex and physical activity trajectories.
#' The function uses a loop to iterate through the ages in the given range, calling the getcycleHS() 
#' function to calculate the risk of the condition for each age, sex, and METs combination. 
#' It then compares the risk against a random matrix of values between 0 and 1, 
#' creating a logical matrix that denotes whether an incident case of the condition has occurred in the period for each person. 
#' The function then stores the matrix of incidence for each age in an array and returns the full array at the end of the loop.
#' 
#' @param age_start starting age for simulation
#' @param age_end ending age for simulation
#' @param v_female is a binary variable representing sex of the population
#' @param v_mets_population amount of metabolic equivalent of task (METs) of physical activity participated by population 
#' @param v_mets_benchmark amount of metabolic equivalent of task (METs) of physical activity participated by benchmark group
#' @param v_constant constant variable used in the risk calculation
#' @param v_ln variable in the risk calculation
#' @param df_inc dataframe of incidence of conditions by age and sex
#' @return an array with dimensions of (age_end, n_i, n_hs) which denotes whether an incident case of the condition has occurred in the period for each person
#' @example
runIncidenceModel <- function(
    age_start,
    age_end,
    v_female,
    df_PA_traj_baseline,
    df_PA_traj_population,
    v_constant,
    v_ln,
    df_inc){
  
  # check lengths of vectors match
  n_i <- ncol(df_PA_traj_baseline)
  assertthat::are_equal(n_i, ncol(df_PA_traj_population))
  n_hs <- length(v_constant)
  assertthat::are_equal(n_hs, length(v_ln))
  
  # create an array which has rows equal to max age, columns equal to the number
  # of individuals and a 3d based on the number of conditions in v_constant. 
  a_hs <- array(data = NA, 
                dim = c(age_end, n_i, n_hs),
                dimnames = list(1:age_end,
                                1:n_i,
                                names(v_ln))) # names from vector risk
  
  # loop through the cycles, adding the matrix of incidence for each cycle
  for(a in age_start:age_end){
    a_hs[a,,] <- 
      getcycleHS(age = a,
                 v_female = v_female,
                 v_mets_population = df_PA_traj_population[a, ],
                 v_mets_benchmark = df_PA_traj_baseline[a, ],
                 v_constant = v_constant,
                 v_ln = v_ln,
                 df_inc = df_inc)
  }
  
  # return the full array
  return(a_hs)
  
}

#' Propagate ones for duration d on the first existence of an incidence event.
#'
#'
#' The propagate_ones() function takes in a numeric vector of 1s and zeros as an input,
#' and assigns the next duration values after the first 1 with a 1 as well. 
#' The function starts by identifying the index where the value of x is 1.
#' Then, it assigns the next duration periods with the condition. 
#' If there is no 1's in the input it would return the same input vector. 
#' This function returns the updated numeric vector assuming length of duration. 
#' This function helps when you want to spread the effect of an incidence event for 
#' a certain duration, for example, in a medical study if someone is diagnosed 
#' with a certain condition the likelihood of the condition being present in the 
#' following days is higher, this function could be used to model that scenario.
#' @param v_x numeric vector of 1s and zeros
#' @param duration duration of the incident event
#' @return a numeric vector of 1s and zeros with updated values assuming length of duration
#' @examples
#' x <- rep(0,20)
#' x[5] <- 1
#' x[1] <- NA
#' propagate_ones(x, 10)
propagate_ones <- function(v_x, 
                           duration) {
  # get length of x
  n <- length(v_x)
  # identify index where value of x is 1
  v_one_index <- which(v_x == 1)
  # assign next 'duration' periods with condition
  if (length(v_one_index) > 0) {
    # identify minimum value
    v_one_min <- min(v_one_index)
    # prevalent_period 
    v_prev_period_max <- min(v_one_min + duration, n)
    # assign ones into this period
    v_x[v_one_min:v_prev_period_max] <- 1
  }
  return(v_x)
}

calcPrev <- function(a_hs, v_durations){
# health prevalence is then created as an array of NAs.
a_prev <- (a_hs[, , ] * NA)

# Loop through conditions, in each case estimating prevalence
for (cond in names(a_hs[1,1,])) {
  #print(v_durations[cond])
  a_prev[, , cond] <- apply(
    X = a_hs[, , cond],
    MARGIN = 2,
    FUN = propagate_ones,
    duration = v_durations[cond]
  )
}

return(a_prev)
}




#' Calculate matrix of mortality risk for a given age and population
#' 
#' This function calculates a matrix of mortality risk for a given age and population.
#' The function first extracts the names of the health states from the first row of the m_prev matrix. 
#' It then creates a vector of age and sex combinations by concatenating the age value and the v_female vector. 
#' Using the v_age_sex vector and the v_hs_names, it extracts the condition mortality and other cause mortality 
#' for the given cycle and population from the m_cond_mort matrix.
#' It then multiplies these values by the health state membership of each person, 
#' stored in m_prev, to get a matrix of survival probabilities for each health state. 
#' Finally, it calculates the total mortality risk for each individual by taking the 
#' product of the survival probabilities across all health states and multiplying it 
#' by 1 minus the other cause mortality rate. The function returns a vector of 
#' mortality rates by health state.
#' 
#' @param age numeric, the age (cycle) of the population
#' @param v_female logical, vector of logicals 1 for female (vs 0 for male)
#' @param m_cond_mort matrix of mortality rates for condition by age, a column for each condition and for OCM
#' @param m_prev matrix of logicals, a column for HS membership for each person (row)
#' @return a vector of mortality rates by health state
#' @export
# calcMortRisk_matrix(age = 25, 
#                     v_female = rbinom(n = 100, size = 1, prob = 0.5),
#                     m_cond_mort = matrix(data = runif(n = 60, min = 0.02, max = 0.2),
#                                          nrow = 10, 
#                                          ncol = 6, 
#                                          dimnames = list(paste0(rep(21:25, 2), "_", rep(0:1, each = 5)),
#                                                          c("BC", "CC", "T2D", "IS", "IHD", "OCM"))),
#                     m_prev = matrix(data = rbinom(n = 500, size = 1, prob = 0.5),
#                                     nrow = 100,
#                                     ncol = 5,
#                                     dimnames = list(1:100,
#                                                     c("BC", "CC", "T2D", "IS", "IHD"))))
calcMortRisk_matrix <- function(age, 
                                v_female,
                                m_cond_mort,
                                m_prev){ 
  # get names of health states
  v_hs_names <- names(m_prev[1,])
  # get age sex combinations
  v_age_sex <- paste0(age, "_", v_female)
  # get condition mortality and other cause mortality for cycle and population
  m_HS_MR_cycle <- m_cond_mort[v_age_sex, v_hs_names]
  v_OCM         <- m_cond_mort[v_age_sex, "OCM"]
  # multiply by health state membership
  m_pS <- as.matrix(x = 1 - m_HS_MR_cycle * m_prev[ ,v_hs_names])
  # estimate total mortality risk for each individual
  # 1 - HS_survival * OCM_survival
  v_MR <- (1 - matrixStats::rowProds(x = m_pS, rows = T)) + v_OCM
  # check that all values of MR are greater than 0
  assertthat::assert_that(all(v_MR >= 0 & v_MR <= 1, na.rm = T),
                          msg = "error, mortality rate not between 0 and 1")
  # return a vector of mortality rates by health state
  return(v_MR)
}

#' Remove health states other than dead from those who are dead.
#'
#' @param m_dead matrix containing, for each cycle (row) who (column) is dead.
#' @param a_prev the array of prevalence, first dimension is cycle, second person, third health state.
#'
#' @return Returns the adjusted prevalence array with health states other than dead removed from those who are dead.
#'
#' @examples
# removeHSfromDEAD(m_dead = matrix(c(0, 0, 1, 1), 
#                                  nrow = 2, 
#                                  ncol = 2), 
#                  a_prev = array(rbinom(n = 24, 
#                                        size = 1, 
#                                        prob = 0.5), 
#                                 dim = c(2, 2, 3)))
#'
removeHSfromDEAD <- function(m_dead, 
                             a_prev){
  
  # for each person multiply health states by survival status
  for (x in 1:dim(a_prev)[2]){
    a_prev[, x, ] <- (1 - m_dead[, x]) * a_prev[, x, ]
  }
  # return the adjusted prevalence array.
  return(a_prev)
}





#' @title runEpiModel
#' @description This function runs an epidemiological model with the given inputs, which include age range, gender, physical activity data, and data on incidence, duration, and mortality risk. It returns an array of health states, an array of prevalence, and a matrix of dead/alive.
#' @param age_start integer, start age of the population
#' @param age_end integer, end age of the population
#' @param v_female logical, vector indicating gender of each individual
#' @param df_PA_traj_baseline dataframe, physical activity data for the baseline group
#' @param df_PA_traj_intervention dataframe, physical activity data for the intervention group
#' @param v_constant numeric, constant term for incidence calculation
#' @param v_ln numeric, natural log term for incidence calculation
#' @param df_inc dataframe, incidence data
#' @param v_durations numeric, vector of durations for each condition
#' @param df_cond_mort dataframe, data on conditional mortality risk
#' @return list, an array of health states, an array of prevalence, and a matrix of dead/alive
runEpiModel <- function(age_start,
                        age_end,
                        v_female,
                        df_PA_traj_baseline,
                        df_PA_traj_intervention,
                        v_constant,
                        v_ln,
                        df_inc,
                        v_durations,
                        df_cond_mort) {
  # check that physical activity levels are above 0 in all cases.
  assertthat::assert_that(all(df_PA_traj_baseline > 0,na.rm = T),
                          msg = "Negative PA values in the baseline group")
  assertthat::assert_that(all(df_PA_traj_intervention > 0,na.rm = T),
                          msg = "Negative PA values in the baseline group")
  
  # Run the incidence model for the population
  a_hs <- runIncidenceModel(
    age_start = age_start,
    age_end = age_end,
    v_female = v_female,
    df_PA_traj_population = df_PA_traj_intervention,
    df_PA_traj_baseline = df_PA_traj_baseline,
    v_constant = v_constant,
    v_ln = v_ln,
    df_inc = df_inc
  )
  
  # calculate prevalence from incidence
  a_prev <- calcPrev(a_hs = a_hs, 
                     v_durations = v_durations)
  # health prevalence is then created as an array of NAs.
  #a_prev <- (a_hs[, , ] * NA)
  #
  ## Loop through conditions, in each case estimating prevalence
  #for (cond in names(a_hs[1,1,])) {
  #  #print(v_durations[cond])
  #  a_prev[, , cond] <- apply(
  #    X = a_hs[, , cond],
  #    MARGIN = 2,
  #    FUN = propagate_ones,
  #    duration = v_durations[cond]
  #  )
  #}
  
  # if prevalence is identical to incidence then print warning.
  if (all(a_prev == a_hs, na.rm = T)) {
    message("Warning: Prevalence and incidence are equal, was this intended?")
  }
  
  # create a matrix to store the mortality risk for each individual in each cycle.
  m_MR <- a_prev[,,1] * NA
  # loop through, filling in the matrix for baseline
  for(x in 1:nrow(m_MR)){
    m_MR[x, ] <- calcMortRisk_matrix(
      age = x,
      v_female = v_female,
      m_cond_mort = base::as.matrix(df_cond_mort),
      m_prev = a_prev[x, , ]
    )
  }
  # Generate random matrix
  m_rand <- genRandMatrix(m_structure = m_MR, 
                          lower = 0, 
                          upper = 1)
  # assign deaths based on matrix of deaths
  m_D <- m_MR > m_rand
  # propagate 1s through the model, death is an absorbing state.
  m_D <- apply(
    simplify = T,
    X = m_D * 1,
    MARGIN = 2,
    FUN = propagate_ones,
    duration = 100
  )
  # adjust health states for death - can't be dead and in a health state
  a_prev <- removeHSfromDEAD(m_dead = m_D, 
                             a_prev =  a_prev)
  # returns a list containing the array of health states, array of prevalence, and matrix of dead/alive
  return(list("a_hs" = a_hs,
              "a_prev" = a_prev,
              "m_D" = m_D))
  
}

