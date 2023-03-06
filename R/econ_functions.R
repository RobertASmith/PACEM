# FUNCTIONS ----

# ESTIMATE AGE-SPECIFIC UTILITIES...
#' @title Utility calculations for general and healthy population
#' @description The function takes in vectors of ages, utilities for general population and healthy population, 
#' the percentage of health state burden at age which is due to conditions modelled versus all conditions, and age range to interpolate to.
#' It creates a dataframe with two columns, 'gp' for general population and 'hp' for healthy population, 
#' linearly interpolates the utilities from data provided, and estimates the utility for those with none of the conditions modelled but may have other conditions. 
#' @param v_age a vector of ages for which utility data provided
#' @param v_U_gp a vector of general population utility for each age
#' @param v_U_H a vector of healthy (no conditions) utility for each age
#' @param HS_burden percentage of the health state burden at age which is due to conditions modelled vs all conditions
#' @param age_range age range to interpolate to.
#' @return A dataframe with columns 'gp', 'hp', and 'no_HS' for general population, healthy population, and population with none of the conditions modelled respectively.
#' @examples 
#' v_age <- c(15,20,30,40,50)
#' v_U_gp <- c(0.6,0.7,0.8,0.9,0.95)
#' v_U_H <- c(0.9,0.95,0.97,0.98,0.99)
#' HS_burden <- 0.2
#' age_range <- 1:70
#' calcPopUtils(v_age, v_U_gp, v_U_H, HS_burden, age_range)
calcPopUtils <- function(v_age, # a vector of ages for which utility data provided
                         v_U_gp, # a vector of general population utility for each age
                         v_U_H, # a vector of healthy (no conditions) utility for each age
                         HS_burden, # percentage of the health state burden at age which is due to conditions modelled vs all conditions
                         age_range){ # age range to interpolate to.
  # create a dataframe with two columns, 'gp' for general pop and hp for healthy pop
  # linearly interpolate from data provided
  df_U <- data.frame("gp" = approx(
    x = v_age,
    y = v_U_gp,
    xout = age_range,
    rule = 2,
    method = "linear"
  )$y,
  
  "hp" = approx(
    x = v_age,
    y = v_U_H,
    xout = age_range,
    rule = 2,
    method = "linear"
  )$y)
  # use the HS burden to estimate utility for those with none of conditions modelled
  # but may have other conditions
  df_U$no_HS  <- df_U$hp * HS_burden + df_U$gp * (1 - HS_burden)
  # check between 0 and 1
  assertthat::assert_that(all(df_U <= 1 & df_U >= 0))
  
  return(df_U)
}

#' @title Utility calculations for a single age cohort
#' @description The function takes in a method of calculation ("multiplicative" or "minimum"),
#' It multiplies the matrix of health states by the health state multipliers, replaces any zeros with ones,
#' and calculates utilities for each individual in the cohort using either multiplicative or minimum method.
#' @param method a string indicating the method of calculation, "multiplicative" or "minimum"
#' @param u_hp a scalar value of utility for healthy population
#' @param v_HS_mult a vector of health state multipliers
#' @param m_HS a matrix of health states, rows indicate individuals and columns indicate health states
#' @param v_D a vector indicating death status of each individual, 1 if dead and 0 if alive
#' @return A vector of health stateutilities for each individual in the cohort
#' @examples 
# method <- "multiplicative"
# u_hp <- 0.9
# v_HS_mult <- c(0.9,0.8,0.7,0.6)
# names(v_HS_mult) <- LETTERS[1:3]
# m_HS <- matrix(c(1,1,1,0,1,0,1,1,0,0,0,1), nrow = 4, ncol = 3, byrow = TRUE, dimnames = list(1:4,LETTERS[1:3]))
# v_D <- c(0,0,1,0)
# calcCycleUtils(method, u_hp, v_HS_mult, m_HS, v_D)

calcCycleHSUtils <- function(method,
                             u_hp, # utility of the baseline healthy population (i.e. those without the conditions modelled)
                             v_HS_mult, # vector of health state multipliers (don't vary with time, these are decrements for utility)
                             m_HS, # matrix of health states (rows individuals, columns HS)
                             v_D){ # vector with whether each individual is dead
  
  # multiply matrix of health states by health state multipliers...
  m_HS_util <- m_HS * v_HS_mult[colnames(m_HS)]
  # replace the zeros with ones
  m_HS_util[m_HS_util == 0] <- 1
  
  switch(method,
         "multiplicative" = {
           # if method is multiplicative find the product of each row.
           v_utils <-
             matrixStats::rowProds(m_HS_util) * u_hp * (1 - v_D) # this needs to be changed to be weighted by PSA.
           
         },
         "minimum" = {
           # if method is minimum find the minimum utility in each row
           v_utils <-
             pmin(matrixStats::rowMins(m_HS_util* u_hp), u_hp * (1 - v_D))
           
         },
         print("error"))
  # note death is assumed to have 0 utility in code
  return(v_utils)
}

#' Calculates the total cycle utilities for each individual.
#' This is used in calculating total QALYs, but can also be used to compare against general population utility values.
#' @param method either "multiplicative" or "minimum"
#' @param v_u_hp a vector of utilities for healthy person (i.e. no health states in the model)
#' @param m_D a matrix with individuals on columns and time on rows, binary alive/dead
#' @param a_prev an array for prevalence: time, person, health state
#' @param v_HS_mult a vector for health state multipliers on utility.
#' @return a matrix of total utility scores across all individuals
#' @examples
# calcTotalCycleUtils(method = "multiplicative",
#                     v_u_hp = rep(1, times = 70),
#                     m_D = matrix(data = 0, nrow = 70, ncol = 5),
#                     a_prev = array(data = 1, dim = c(70, 5, 2), dimnames = list(1:70, 1:5, c("condition1", "condition2"))),
#                     v_HS_mult = c(condition1 = 0.8, condition2 = 0.6))
#' @export
calcTotalCycleUtils <- function(method, # either multiplicative or minimum
                                v_u_hp,    # vector of utilities for healthy person (i.e. no health states in model)
                                m_D, # matrix with individuals on column and time on row, binary alive dead
                                a_prev, # array for prevalence: time, person, health state
                                v_HS_mult) { # vector for HS multipliers on utility.
  
  assertthat::are_equal(nrow(m_D), dim(a_prev)[1])
  
  # create empty array ready to insert cycle utils
  m_total_util <- array(data = NA, dim = dim(m_D))
  
  # loop through each age, adding in the individuals' (column) utility scores
  # in each cycle (age).
  for(age in 11:nrow(m_total_util)) {
    m_total_util[age,] <- calcCycleHSUtils(
      method = method,
      u_hp = v_u_hp[age],
      v_HS_mult = v_HS_mult,
      m_HS = a_prev[age, , ],
      v_D = m_D[age, ]
    )
  }
  # returns total utility scores across all individuals
  return(m_total_util)
}



#' Calculates the undiscounted per-cycle costs over the time horizon of the model.
#' Includes both first year costs - on the year of diagnosis - and the subsequent year costs.
#' @param v_first_year_costs a vector of first year health state costs
#' @param v_subsequent_year_costs a vector of subsequent year health state costs
#' @param a_hs an array of health state values (binary) for developing condition in year
#' @param a_prev an array of health state values for having condition
#' @return a matrix of total costs, one column for each condition one row for each cycle
#' @examples
# calcCycleCosts(v_first_year_costs = c(condition1 = 100, condition2 = 200),
#               v_subsequent_year_costs = c(condition1 = 50, condition2 = 75),
#               a_hs = array(data = 1, dim = c(5,5,2), dimnames = list(1:5, 1:5, c("condition1", "condition2"))),
#               a_prev = array(data = 1, dim = c(5,5,2), dimnames = list(1:5, 1:5, c("condition1", "condition2"))))
calcCycleCosts <- function(v_first_year_costs, # a vector of first year costs
                           v_subsequent_year_costs, # a vector of subsequent year health state costs
                           a_hs, # an array of health state values (binary) for developing condition in year
                           a_prev) { # an array of health state values for having condition
  
  # check names match
  assertthat::are_equal(names(v_first_year_costs) , names(v_subsequent_year_costs))
  assertthat::assert_that(all(names(v_first_year_costs) %in% dimnames(a_prev)[[3]]))
  
  # Calculate incident costs by year and condition (total for all individuals)
  m_incident_costs <- sapply(
    X = names(v_first_year_costs),
    FUN = function(x, a_hs, v_first_year_costs) {
      base::rowSums(a_hs[, , x] * v_first_year_costs[x], na.rm = T)
    },
    a_hs = a_hs,
    v_first_year_costs = v_first_year_costs
  )
  # calculate prevalent costs in same way
  m_prevalent_costs <- sapply(
    X = names(v_subsequent_year_costs),
    FUN = function(x, a_prev, v_subsequent_year_costs) {
      base::rowSums(a_prev[, , x] * v_subsequent_year_costs[x], na.rm = T)
    },
    a_prev = a_prev,
    v_subsequent_year_costs = v_subsequent_year_costs
  )
  # combine for total costs
  m_total_costs <- m_incident_costs + m_prevalent_costs
  # returns a matrix of total costs, one column for each condition one row for each cycle
  return(m_total_costs)
  
}



#' @title Run the economic part of the model to estimate total discounted costs and QALYs
#' @param v_u_hp Vector containing the health state utilities for those with none of the conditions included in the model
#' @param utility_method Character describing method used to combined health state utilities in context of decrements, "minimum" or "multiplicative"
#' @param m_D Matrix of death status - 1 is dead for cycle (row) and person (column)
#' @param a_prev A 3D array for existence of health condition: time, person, health state
#' @param v_HS_mult A vector of health state multipliers (must be between 0 and 1) - used to estimate total utility.
#' @param v_first_year_costs A vector of first year health state costs
#' @param v_subsequent_year_costs A vector of subsequent year health state costs
#' @param a_hs An array of health state values (binary) for developing condition in year
#' @param dr_costs Annual discount rate for costs
#' @param dr_qalys Annual discount rate for QALYs
#' @return A vector of results with the sum of discounted (and undiscounted) costs and QALYs
#' @export
#' @examples
runEconModel <- function(v_u_hp, 
                         utility_method = "multiplicative", 
                         m_D,   
                         a_prev, 
                         v_HS_mult, 
                         v_first_year_costs,  
                         v_subsequent_year_costs, 
                         a_hs,  
                         dr_costs, 
                         dr_qalys){ 
  
  # calculate total cycle utils for each individual
  m_util_cycle_ind <- calcTotalCycleUtils(
    v_u_hp = v_u_hp,
    method = utility_method,
    m_D = m_D,
    a_prev = a_prev,
    v_HS_mult = v_HS_mult)
  # take the rowsums to get total per cycle utils
  v_util_cycle     <- rowSums(m_util_cycle_ind, na.rm = T)
  
  
  # Get a matrix of cycle costs (rows) for each condition (columns).
  m_costs_cycle_condition <-
    calcCycleCosts(
      v_first_year_costs = v_first_year_costs,
      v_subsequent_year_costs = v_subsequent_year_costs,
      a_hs = a_hs,
      a_prev = a_prev
    )
  # take the rowsums to get total per cycle costs
  v_costs_cycle <- rowSums(x = m_costs_cycle_condition, na.rm = T)
  
  # create a set of results with the sum of discounted (and undiscounted) costs and qalys
  #v_results <-
  #  c(disc_costs   = sum(discVec(v_costs_util[11:80], dr = dr_costs)),
  #    undisc_costs = sum(discVec(v_costs_util[11:80], dr = 0)),
  #    disc_qalys   = sum(discVec(v_util_cycle[11:80], dr = dr_qalys)),
  #    undisc_qalys = sum(discVec(v_util_cycle[11:80], dr = 0)))
  
  # return this set of results
  return(cbind("util" = v_util_cycle, 
               "costs" = v_costs_cycle))
}

#' Creates a discount vector based on a discount rate and number of cycles.
#' @param dr a decimal representing the discount rate
#' @param cycles the number of time periods to create the discount vector for
#' @return a vector of discount weights, where the first element is 1 and each subsequent element is discounted by the rate
#' @examples
#' makeDiscountVec(dr = 0.03, cycles = 10)
#' @export
makeDiscountVec <- function(dr,
                            cycles) {
  # create a vector of time periods starting from
  # 0 and increasing to T-1
  v_t <- 0:(cycles - 1)
  # use the discount rate (a decimal) and time vector to create a discount vector
  v_dw <- 1 / ((1 + dr) ^ v_t)
  return(v_dw)
}

#' Creates a discounted vector of costs using a discount rate and the input costs vector
#' @param v_costs a vector of costs
#' @param dr a decimal representing the discount rate, defaults to 0.03
#' @return a vector of discounted costs
#' @examples
#' discVec(v_costs = rep(1,10), dr = 0.03)
#' @export
discVec <- function(v_costs,
                    dr = 0.03) {
  # generate a vector of discount weights.
  v_dw <- makeDiscountVec(dr = dr,
                          cycles = length(v_costs))
  # multiply the vector of costs by the vector of discount rates.
  temp <- v_costs * v_dw
  return(temp)
}

#' Calculate an ICER from costs and effects for the baseline and intervention arms
#'
#' @param c_base  costs in the base scenario
#' @param c_int   costs in the intervention scenario
#' @param e_base  effects in the base scenario
#' @param e_int   effects in the intervention scenario
#'
#' @return a single numeric value
#'
#' @examples calcICER(c_base = 1000, c_int = 2900, e_base = 10, e_int = 20)
calcICER <- function(c_base, 
                     c_int, 
                     e_base,
                     e_int){
  assertthat::assert_that(all(is.numeric(c(c_base, c_int, e_base, e_int))))
  # simple equation for the ICER
  ICER <- (c_int - c_base) / (e_int - e_base)
  return(as.numeric(ICER))
}





#' @title Discount PSA Outcome
#'
#' @description This function calculates the discounted PSA outcome by discounting a matrix of values over time.
#'
#' @param m_outcome matrix of outcome values, column for each PSA iteration, row for cycle.
#' @param dr discount rate.
#'
#' @return A vector of the net present value for each PSA iteration
#'
#' @examples
#' m_outcome <- matrix(runif(300,0,1), ncol = 3)
#' dr <- 0.1
#' discPSAoutcome(m_outcome, dr)
#' 
#' @export
#'
discPSAoutcome <- function(m_outcome, dr) {
  apply(
    X = m_outcome,
    MARGIN = 2,
    FUN = function(x)
      sum(discVec(x, dr = dr))
  )
}


mean_with_predecessor <- function(m) {
  # initialize output vector
  m_output <- m * NA
  
  # set predecessor of first value to 0
  m_output[1,] <- m[1,]/2
  
  # compute mean with predecessor for all other values
  m_output[2:nrow(m),] <- (m[1:(nrow(m)-1),] + m[2:nrow(m),])/2
  
  # return output
  return(m_output)
}

#' @title Calculate PSA Outcomes
#'
#' @description This function calculates the net present value of Costs and QALYs for cycles 
#'  from age 11 to 80 for both the baseline and intervention group.
#'
#' @param l_PSA_results a list of PSA results, must include arrays for 'base' and 'int' 
#' by age (rows) and costs and util as columns.
#'
#' @return A data-frame with four columns: 
#' costs for the baseline group, 
#' costs for the intervention group, 
#' QALYs for the baseline group, 
#' QALYs for the intervention group.
#'
#' @examples
#' l_PSA_results <- list(base = matrix(c(100, 90, 80, 70, 60), ncol = 3), 
#'                       int = matrix(c(100, 90, 80, 70, 60), ncol = 3), 
#'                       params = list(dr_costs = 0.1, dr_utils = 0.05))
#' calcPSAoutcomes(l_PSA_results)
#'
#' @export
#'
calcPSAoutcomes <- function(l_PSA_results) { # a list of PSA results, must include elements for 'base' and 'int'
  # Costs
  v_C_base <-
    discPSAoutcome(m_outcome = mean_with_predecessor(l_PSA_results$base[11:80, "costs", ]),
                   dr = l_PSA_results$params["dr_costs"])
  v_C_int <-
    discPSAoutcome(m_outcome = mean_with_predecessor(l_PSA_results$int[11:80, "costs", ]),
                   dr = l_PSA_results$params["dr_costs"])
  
  # Utilities
  v_U_base <-
    discPSAoutcome(m_outcome = mean_with_predecessor(l_PSA_results$base[11:80, "util", ]),
                   dr = l_PSA_results$params["dr_qalys"])
  v_U_int <-
    discPSAoutcome(m_outcome = mean_with_predecessor(l_PSA_results$int[11:80, "util", ]),
                   dr = l_PSA_results$params["dr_qalys"])
  
  # return a data-frame with costs & qaly for the baseline group and the intervetention group (4 columns)
  return(data.frame(v_C_base, v_C_int, v_U_base, v_U_int))
}





makeCEAC <- function(total_costs = example_TC,
                     total_qalys = example_TQ,
                     treatment = c("treat 1","notreat"),
                     lambda_min = 0,
                     lambda_max = 50000){
  
  all_names = colnames(total_costs)
  
  # assign each treatment a colour & give name based on column
  legend_colors = rainbow(n = ncol(total_costs))
  names(legend_colors) = colnames(total_costs)
  
  # assign each treatment a dash
  dash_numbers = 1:ncol(total_costs)
  names(dash_numbers) = colnames(total_costs)
  
  
  # assign colours and dashes based on treatment choices
  legend_colors = legend_colors[names(legend_colors) %in% treatment]
  legend_dash   = dash_numbers[names(dash_numbers) %in% treatment]
  
  # Take the appropriate columns from cost and qalys matrices/df
  v_TC = total_costs[,colnames(total_costs) %in% treatment]
  v_TQ = total_qalys[,colnames(total_qalys) %in% treatment]
  
  # Create a vector of 100 lambda values
  lambdas <- seq(from = lambda_min,
                 to = lambda_max,
                 length.out = 1000)
  
  # ProbCE for each Lambda - slow but reliable
  df_CEAC = c()
  
  for(l in lambdas){
    
    # get net benefit for each strategy
    nb = v_TQ * l - v_TC
    
    # idenfity whether each strategy is the maximum or not
    nb = apply(nb,1,function(x) x == max(x))
    
    # make data-frame of probability each option is CE
    nb = data.frame(Intervention = colnames(v_TC),
                    lambda = l,
                    value = apply(nb,1,mean))
    
    # bind to previous lambda data
    df_CEAC = rbind(df_CEAC,nb)
  }
  
  
  # MAKE PLOT
  
  ggplot2::ggplot(data = df_CEAC,
                  ggplot2::aes(x = lambda,y= value,col = Intervention)
  )+
    
    # theme
    ggplot2::theme_minimal() +
    
    # legend
    ggplot2::theme(legend.position = "top",
                   legend.text = ggplot2::element_text(size=11),
                   legend.title = ggplot2::element_text(size=11),
                   title = ggplot2::element_text(size=11)) +
    
    # lines
    ggplot2::geom_line(size=1.2) +
    
    # y axis
    ggplot2::scale_y_continuous(breaks=seq(0,1,0.25),limits = c(0,1),name = "Probability most cost-effective") +
    
    # labels
    ggplot2::xlab(label = "Willingness-to-pay (GBP)")+
    
    ggplot2::labs(title = "Cost Effectiveness Acceptability Curves",
                  subtitle = "The probability each preferred intervention is most cost effective against willingness to pay for each QALY threshold.") +
    
    # apply color scheme
    ggplot2::scale_color_manual(name="Treatment",values = legend_colors) +
    NULL
  
  
  
}
