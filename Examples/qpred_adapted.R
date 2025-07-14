qpred.idr_2 <- function(predictions, quantiles) {
  
  # Check input
  if (!is.vector(quantiles, "numeric") || min(quantiles) < 0 ||
      max(quantiles) > 1) 
    stop("quantiles must be a numeric vector with entries in [0,1]")
  q0 <- function(data) {
    # Evaluate quantile function (stepfun) at given quantiles 
    stats::stepfun(x = data$cdf,
                   y = c(data$points, data$points[nrow(data)]), right = FALSE)(quantiles)
  }
  qVals <- lapply(predictions, q0)
  do.call(rbind, qVals)
}