
summarize_flexsurv <- function(fit, compute_HR = TRUE) {
  # 1. Extract res matrix
  res <- fit$res
  if(is.null(res)) stop("No results found in fit$res")
  
  # 2. Convert to data frame
  df <- as.data.frame(res)
  df$term <- rownames(df)
  rownames(df) <- NULL
  
  # 3. Filter only regression parameters (exclude distribution parameters)
  df <- df[!df$term %in% c("mu","sigma","Q"), ]
  
  # 4. Compute Wald p-values
  df$p_value <- 2 * (1 - pnorm(abs(df$est / df$se)))
  
  # 5. Optional: exponentiated estimates (time ratios)
  if(compute_HR) df$HR <- exp(df$est)
  
  # 6. Round numeric columns
  num_cols <- c("est","L95%","U95%","se","p_value","HR")
  for(col in num_cols) if(col %in% colnames(df)) df[[col]] <- round(df[[col]],3)
  
  # 7. Reorder columns nicely
  df <- df %>% select(term, est, se, `L95%`, `U95%`, HR, p_value)
  
  return(df)
}


lrt <- 2 * (fit_full$loglik - fit_reduced$loglik)
df <- fit_full$npars - fit_reduced$npars
p_value <- 1 - pchisq(lrt, df)

p_value