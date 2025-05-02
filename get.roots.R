get.roots <- function(order, beta, type_interp = "linear") {
  if(!(order %in% c(1,2,3,4))) {
    stop("order must be one of the values 1,2,3,4.")
  }
  if (beta > 2) {
    beta <- beta - floor(beta - 1)
  }
  mt <- get(paste0("m", order, "table"))
  rb <- rep(0, order + 1)
  rc <- rep(0, order)
  if(type_interp == "linear"){
    if(order == 1) {
      rc = approx(mt$beta, mt[[paste0("rc")]], beta)$y
    } else {
      rc = sapply(1:order, function(i) {
        approx(mt$beta, mt[[paste0("rc.", i)]], beta)$y
      })
    }
    rb = sapply(1:(order+1), function(i) {
      approx(mt$beta, mt[[paste0("rb.", i)]], xout = beta)$y
    })
    factor = approx(mt$beta, mt$factor, xout = beta)$y
  } else if(type_interp == "spline") {
    if(order == 1) {
      rc = spline(mt$beta, mt[[paste0("rc")]], xout = beta)$y
    } else {
      rc = sapply(1:order, function(i) {
        spline(mt$beta, mt[[paste0("rc.", i)]], xout = beta)$y
      })
    }
    rb = sapply(1:(order+1), function(i) {
      spline(mt$beta, mt[[paste0("rb.", i)]], xout = beta)$y
    })
    factor = spline(mt$beta, mt$factor, xout = beta)$y
  } else {
    stop("invalid type. The options are 'linear' and 'spline'.")
  }
  return(list(rb = rb, rc = rc, factor = factor))
}

m1table <- rSPDE:::m1table
m2table <- rSPDE:::m2table
m3table <- rSPDE:::m3table
m4table <- rSPDE:::m4table

roots <- get.roots(order = 2, beta = 1.5, type_interp = "linear")
