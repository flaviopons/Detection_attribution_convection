library(fitdistrplus) # Gamma fit
library(SuppDists) # inverse gaussian
library(sn) #skew normal
library(GeneralizedHyperbolic) #NIG distribution
library(statmod)  # inverse gaussian

fit_prob_hazard <- function(Y, Y_ref, distr){
  
  # this function fits a distribution to the hazard variable at each grid point
  # and compute the probability of the hazard for the event 
  
Y_2d <- matrix(Y, nrow = dim(Y)[1]*dim(Y)[2], ncol = dim(Y)[3])

# for Gamma e Inverse Gaussian: the fit is only on positive values, zeros must be handled
# using a delta in zero.

if (distr == 'Gamma'){
  fitted_params_Y <- apply(Y_2d, 1, function(series) {
    if (length(series[series>0]) > 0){
      fit <- fitdist(series[series>0], distr = "gamma", method = "mme")
      return(fit$estimate)  # Extract the fitted parameters
    }
    else return(c(NA,NA))
  })
  # Reshape the results back into a grid format
  fitted_params_Y_grid <- array(fitted_params_Y, dim = c(2, dim(Y)[1], dim(Y)[2]))
    
  prob_zero <- apply(Y_2d, 1, function(series) {
    length(series[series==0])/length(series)
  })
    
  prob_pos <- sapply(1:dim(fitted_params_Y)[2], function(k) {
    params <- fitted_params_Y[, k]
    pgamma(as.vector(Y_ref)[k], params[1], params[2], lower.tail = FALSE)
  })
  prob_event_vec <- prob_zero + (1-prob_zero)*prob_pos
}

if (distr == 'IG'){
  fitted_params_Y <- apply(Y_2d, 1, function(series) {
    if (length(series[series>0])>0) fit <- glm(series[series>0] ~ 1, family = inverse.gaussian(link = "log"))
    return(c(exp(coef(fit)[1]), summary(fit)$dispersion))  # Extract the fitted parameters
  })
  
  fitted_params_Y_grid <- array(fitted_params_Y, dim = c(2, dim(Y)[1], dim(Y)[2]))
  
  prob_zero <- apply(Y_2d, 1, function(series) {
    length(series[series==0])/length(series)
  })
  
  prob_pos <- sapply(1:dim(fitted_params_Y)[2], function(k) {
    params <- fitted_params_Y[, k]
    pinvGauss(as.vector(Y_ref)[k], params[1], params[2], lower.tail = FALSE)
  })
  
  prob_event_vec <- prob_zero + (1-prob_zero)*prob_pos
}

if (distr == 'NIG'){
  fitted_params_Y <- apply(Y_2d, 1, function(series) {
    fit <- nigFit(series, method = "Nelder-Mead")
    return(fit$param)
  })
  fitted_params_Y_grid <- array(fitted_params_Y, dim = c(4, dim(Y)[1], dim(Y)[2]))

  prob_event_vec <- sapply(1:dim(fitted_params_Y)[2], function(k) {
    params <- fitted_params_Y[, k]
    pnig(as.vector(Y_ref)[k], params[1], params[2], params[3], params[4], lower.tail = FALSE)
  })
}

if (distr == 'SN'){
  print('Probabilities for the Skew-Normal distributio are given in the lower tail.')
  fitted_params_Y <- apply(Y_2d, 1, function(series) {
    fit <- selm(series ~ 1, family = "SN")    
    return(fit@param$dp) 
  })
  fitted_params_Y_grid <- array(fitted_params_Y, dim = c(3, dim(Y)[1], dim(Y)[2]))
   prob_event_vec <- sapply(1:dim(fitted_params_Y)[2], function(k) {
     params <- fitted_params_Y[,k]
     psn(as.vector(Y_ref)[k], params[1], params[2], params[3])
   })
}

if (distr == 'Normal'){
  fitted_params_Y <- apply(Y_2d, 1, function(series) {
    mm <- mean(series, na.rm = TRUE)   
    sd <- stats::sd(series, na.rm = TRUE)   
    return(c(mm, sd)) 
  })
  fitted_params_Y_grid <- array(fitted_params_Y, dim = c(2, dim(Y)[1], dim(Y)[2]))
  prob_event_vec <- sapply(1:dim(fitted_params_Y)[2], function(k) {
    params <- fitted_params_Y[,k]
    pnorm(as.vector(Y_ref)[k], params[1], params[2], lower.tail = FALSE)
  })
}

prob_event <- matrix(prob_event_vec, nrow = dim(Y)[1], ncol = dim(Y)[2])

fit_result <- list(fitted_params_Y_grid, prob_event)
names(fit_result) <- c('fitted parameters', 'event hazard probability')
return(fit_result)
}
