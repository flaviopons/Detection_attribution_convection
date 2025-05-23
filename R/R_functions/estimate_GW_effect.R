source('/home/estimr3/fpons/R_functions/matrix_to_long.R')
source('/home/estimr3/fpons/R_functions/run_regression.R')
source('/home/estimr3/fpons/R_functions/extract_results_regression.R')

# var_ana = hazard variable in correspondence to analogs, 3D array (nlon x nlat x T)
# X0 = (monthly) time series of independent variable
# time_X0 = dates in correspondence of values of X0
# X = time series of independent variable (GMST, RMST...) corresponding to analogs
# time_X = dates in correspondence of values of X
# date_event = Date of the event under consideration in yyyy-mm-dd format
# year1, year2 = years corresponding to two chosen levels of GMST
# X1, X2 = values of X corresponding to year1, year2
# If year1, year2 are provided, X1 and X2 are retrieved from data.
# If X1, X2 are provided, year1, year2 are retrieved from data.

estimate_GW_effect <- function(var_ana, X0, time_X0, X, time_X, date_event, X1=NULL, X2=NULL, 
                               year1=NULL, year2=NULL, type = 'linear', adj.level = 0.05,
                               adj.method = 'BH'){
  
  # Linear by default. It is a parameter of 'run_regression.R'.
  
  # Check if at least one between reference temperature or reference years is provided
  if (is.null(X1) && is.null(X2) && (is.null(year1) || is.null(year2))) {
    stop("Either (X1 and X2) OR (year1 and year2) must be provided.")
  }
  

  # Case 1: If year1/year2 are provided, extract X1/X2 from X
  if (!is.null(year1) && !is.null(year2)) {
    # Extract years from time_X0
    years_X0 <- as.numeric(format(time_X0, "%Y"))
    
    # Find indices matching year1/year2
    idx1 <- which(month(time_X0)==month(date_event) & year(time_X0)==year1)
    idx2 <- which(month(time_X0)==month(date_event) & year(time_X0)==year2)
    
    # Check if years exist in time_X
    if (length(idx1) == 0 || length(idx2) == 0) {
      stop("year1 or year2 not found in time_X0.")
    }
    
    # Extract X1/X2 (assuming single match per year; modify if needed)
    X1 <- X0[idx1]  # Takes first match if multiple exist
    X2 <- X0[idx2]
  }
  
  T_diff_year1 <- (X-X1)
  T_diff_year1_exp <- sweep(array(1, dim = c(dim(var_ana)[1], dim(var_ana)[2], dim(var_ana)[3])), 
                           3, T_diff_year1, "*")  # Broadcast to MxNxT
  T_diff_year2 <- (X-X2)
  T_diff_year2_exp <- sweep(array(1, dim = c(dim(var_ana)[1], dim(var_ana)[2], dim(var_ana)[3])), 
                           3, T_diff_year2, "*")  # Broadcast to MxNxT
  
  # Reshape the var_ana array
  # Convert the MxNxTT array into a matrix where each row is a time series at a grid point
  var_ana_matrix <- matrix(var_ana, nrow = dim(var_ana)[1] * dim(var_ana)[2], ncol = dim(var_ana)[3])
  
  # Apply the regression function to each row of the matrix
  regression_models <- apply((var_ana_matrix), 1, run_regression, X, type)
  regression_results <- extract_results_regression(regression_models)

    
  # Reshape the results back into an MxNx2 array (2 coefficients: intercept and slope)
  #coefficients_array <- array(regression_results$coefficients, dim = c(2, dim(var_ana)[1], dim(var_ana)[2]))
  #pvals_array <- array(regression_results$pvalues, dim = c(2, dim(var_ana)[1], dim(var_ana)[2]))
  
  alpha <- matrix(regression_results$coefficients[,1],dim(var_ana)[1],dim(var_ana)[2])  # contains the intercepts
  alpha_exp <- array(alpha, dim = c(dim(var_ana)[1], dim(var_ana)[2], dim(var_ana)[3]))  # Replicate beta along the time dimension
  beta <- matrix(regression_results$coefficients[,2],dim(var_ana)[1],dim(var_ana)[2]) # contains the slopes
  beta_exp <- array(beta, dim = c(dim(var_ana)[1], dim(var_ana)[2], dim(var_ana)[3]))  # Replicate beta along the time dimension

  pvals_alpha <- matrix(regression_results$pvalues[,1],dim(var_ana)[1],dim(var_ana)[2]) 
  pvals_beta <- matrix(regression_results$pvalues[,2],dim(var_ana)[1],dim(var_ana)[2]) 
  
  test_alpha <- ifelse(pvals_alpha <= adj.level, 1, NA)
  test_beta <- ifelse(pvals_beta <= adj.level, 1, NA)
  
  pvals_alpha_adj <- matrix(p.adjust(as.vector(pvals_alpha), method = adj.method), nrow = nrow(pvals_alpha))
  pvals_beta_adj <- matrix(p.adjust(as.vector(pvals_beta), method = adj.method), nrow = nrow(pvals_beta))
  
  test_alpha_adj <- ifelse(pvals_alpha_adj <= adj.level, 1, NA)
  test_beta_adj <- ifelse(pvals_beta_adj <= adj.level, 1, NA)
  
  if (type == 'linear'){
  var_ana_detrended <- var_ana - beta_exp*X
  var_ana_year1 <- var_ana - beta_exp*T_diff_year1_exp
  var_ana_year2 <- var_ana - beta_exp*T_diff_year2_exp
  }
  
  if (type == 'log-lin'){
    var_ana_detrended <- var_ana*exp(-beta_exp*X)
    var_ana_year1 <- var_ana*(exp(-beta_exp*T_diff_year1_exp))
    var_ana_year2 <- var_ana*(exp(-beta_exp*T_diff_year2_exp))
  }
  
  if (type == 'tweedie'){
    var_ana_detrended <- var_ana*exp(-beta_exp*X)
    var_ana_year1 <- var_ana*(exp(-beta_exp*T_diff_year1_exp))
    var_ana_year2 <- var_ana*(exp(-beta_exp*T_diff_year2_exp))
  }
  
  if (type == 'sinh'){
    var_ana_detrended <- sinh(asinh(var_ana) - beta_exp * X)
    var_ana_year1 <- sinh(asinh(var_ana) - beta_exp * T_diff_year1_exp)
    var_ana_year2 <- sinh(asinh(var_ana) - beta_exp * T_diff_year2_exp)
    
  }

  GW_effect_results <- list(c(X1, X2), c(year1, year2), alpha_exp, beta_exp, 
                            pvals_alpha, pvals_beta, test_alpha, test_beta, 
                            pvals_alpha_adj, pvals_beta_adj,
                            test_alpha_adj, test_beta_adj, var_ana_detrended, var_ana_year1, 
                            var_ana_year2)
  
  names(GW_effect_results) <- c('reference temperatures', 'reference years', 'alpha_coef',
                                'beta_coef', 'pvals_alpha', 'pvals_beta', 'test_alpha',
                                'test_beta',
                                'pvals_alpha_adj', 'pvals_beta_adj', 'test_alpha_adj',
                                'test_beta_adj', 'detrended hazard', 'hazard temperature 1', 
                                'hazard temperature 2')
  return(GW_effect_results)
    }


