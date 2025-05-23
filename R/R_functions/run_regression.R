library(mgcv)

run_regression <- function(Y_series, X_series, type = 'linear') {
  
  if (all(is.na(Y_series))) {
    return(NULL)
  }
  
  # Use tryCatch to avoid errors when the series in NA (e.g. for SST on land points)
  model <- tryCatch(
    {
      if (type == 'linear') {
        lm(Y_series ~ X_series)
      } else if (type == 'tweedie') {
        gam(Y_series ~ X_series, family = tw(link = "log"))
      } else if (type == 'log-lin') {
        mm <- min(Y_series[Y_series > 0], na.rm = TRUE)
        Y_series[Y_series == 0] <- mm
        lm(log(Y_series) ~ X_series)
      } else if (type == 'sinh') {
        lm(asinh(Y_series) ~ X_series)
      }
    },
    error = function(e) {
      warning("Regression failed for a row: ", e$message)
      NULL  # Restituisce NULL in caso di errore
    }
  )
  
  return(model)
}
