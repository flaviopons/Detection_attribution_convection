extract_results_regression <- function(model_list){
  coef_matrix <- matrix(NA, nrow = length(model_list), ncol = 2)
  pvalue_matrix <- matrix(NA, nrow = length(model_list), ncol = 2)
  
  for(i in seq_along(model_list)) {
    model <- model_list[[i]]
    
    if(inherits(model, "lm") || inherits(model, "gam")) {
      smry <- tryCatch(summary(model), error = function(e) NULL)
      
      if(!is.null(smry)) {
        coefs <- smry$coefficients
        if(nrow(coefs) >= 2) {
          coef_matrix[i, ] <- coefs[1:2, "Estimate"]
          pvalue_matrix[i, ] <- coefs[1:2, "Pr(>|t|)"]
        }
      }
    }
  }
  
  colnames(coef_matrix) <- c("alpha", "beta")
  colnames(pvalue_matrix) <- c("p_alpha", "p_beta")
  
  return(list(coefficients = coef_matrix, pvalues = pvalue_matrix))
}
