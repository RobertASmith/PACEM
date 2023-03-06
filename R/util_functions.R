getRRmean <- function(data = df_rr,
                      condition = "BC",
                      value = "alpha"){
  if(value == "alpha"){
    
    data[data$cond == condition & data$bound == "m", "c"]
    
  }else if(value == "beta"){
    
    data[data$cond == condition & data$bound == "m", "ln"]
    
  } else{
    print("error")
  }
}

getRRsd <- function(data = df_rr,
                    condition = "BC",
                    value = "alpha"){
  
  if(value == "alpha"){ # for alpha values (i.e. c)
    
    # sd is the absolute difference between upper and lower bound divided by 3.92
    abs(data[data$cond == condition & data$bound == "h","c"] - 
          data[data$cond == condition & data$bound == "l","c"])/3.92
    
    
  }else if(value == "beta"){ # for beta values (i.e. LN)
    
    # sd is the absolute difference between upper and lower bound divided by 3.92
    abs(data[data$cond == condition & data$bound == "h","ln"] - 
          data[data$cond == condition & data$bound == "l","ln"])/3.92
    
    
  } else{
    print("error")
  }
}


#' Calculate Standard Deviation from a 95% Confidence Interval
#'
#' This function calculates the standard deviation from a 95% confidence interval.
#'
#' @param lwr a numeric value representing the lower bound of the confidence interval
#' @param upr a numeric value representing the upper bound of the confidence interval
#'
#' @return a numeric value representing the standard deviation.
#'
#' @note This function assumes that the sample size is large enough for the normal approximation to be valid. 
#' The `qt()` function from the stats package is used to calculate the quantile of the standard normal distribution.
#' This function returns the standard deviation of the population assuming the sample is normally distributed.
#'
#' @examples
#' calculate_std_from_ci(lwr = 2, upr = 8)
#' # Output: 2
#'
#' calculate_std_from_ci(lwr = -1.96, upr = 1.96)
#' # Output: 1
#'
#' @export
#'
calculate_std_from_ci <- function(lwr, upr) {
  (upr - lwr) / (2 * qt(p = 0.975, df = Inf))
}


#' Generate Quantiles of a Distribution
#' 
#' @param mean Mean of the distribution
#' @param sd Standard deviation of the distribution
#' @return Vector of 1000 values representing the quantiles of the distribution from 0.001 to 0.999
#' @export
#'
generate_quantiles <- function(mean, sd, length.out){
  prob <- seq(0.01, 0.99, length.out = length.out)
  quantiles <- qnorm(prob, mean = mean, sd = sd)
  return(quantiles)
}

