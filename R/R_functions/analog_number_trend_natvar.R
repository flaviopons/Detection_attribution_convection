library(pscl)
library(MASS)
library(lubridate)

analog_number_trend_natvar <- function(indices, time_tele, months_season, date_ana, crit = 'AIC'){

  df_indices <- data.frame(scale(indices))
  names(df_indices) <- names(indices)
  df_indices$time_tele <- time_tele
  
  ind_months_tele <- which(month(time_tele) %in% months_season)
  df_indices <- df_indices[ind_months_tele,]
  
  # Extract year-month (first day of the month) from the sequence of dates
  ymon_dates <- floor_date(date_ana, "month")
  # Ensure time_tele is also in year-month format
  ymon_tele <- floor_date(df_indices$time_tele, "month")
  # Count the occurrences of dates in each year-month of time_tele
  counts_month <- sapply(ymon_tele, function(month) sum(ymon_dates == month))
  
  # Fit Poisson and Negative Binomial without and with zero inflation
  # Standard models
  pois <- glm(counts_month ~ GST+AMO+PDO+ENSO, data = df_indices, family = poisson)
  nb <- glm.nb(counts_month ~ GST+AMO+PDO+ENSO, data = df_indices)
  # Zero-inflated models (adjust predictors after "|" for zero-inflation component)
  zip <- zeroinfl(counts_month ~ GST+AMO+PDO+ENSO | GST+AMO+PDO+ENSO, data = df_indices,  dist = "poisson")
  zinb <- zeroinfl(counts_month ~ GST+AMO+PDO+ENSO | GST+AMO+PDO+ENSO, data = df_indices, dist = "negbin")
  
  # a possible dispersion test is dividing residual deviance and residual degrees of freedom
  # if larger than 1 negative binomial could be considered
  dispersion_test <- sum(residuals(pois, type = "pearson")^2) / pois$df.residual
  
  # compare AIC and BIC
  aic_models <- c(AIC(pois), AIC(zip), AIC(nb), AIC(zinb))
  bic_models <- c(BIC(pois), BIC(zip), BIC(nb), BIC(zinb))
  
  mods <- list(Poisson = pois, ZIP = zip, NB = nb, ZINB = zinb)
  if(crit == 'AIC') best_model <- names(mods)[which.min(aic_models)]
  if(crit == 'BIC') best_model <- names(mods)[which.min(bic_models)]
  
  result <- list(mods, dispersion_test, aic_models, bic_models, best_model)
  names(result) <- c('models', 'dispersion_test', 'aic', 'bic', 'best_model')
  return(result)
}
