# For each m and beta, this function returns c_m/b_{m+1} and the roots of rb and rc
my.get.roots <- function(order, beta) {
  mt <- get(paste0("m", order, "table"))
  rb <- rep(0, order + 1)
  rc <- rep(0, order)
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
  return(list(rb = rb, rc = rc, factor = factor))
}

# Given the roots, return the polynomial coefficients in increasing order like a+bx+cx^2+...
# So if the output is c(6, -5, 1), that means 6 - 5x + x^2
poly_from_roots <- function(roots) {
  coef <- 1
  for (r in roots) {
    coef <- convolve(coef, c(1, -r), type = "open")
  }
  return(coef) # returns in increasing order like a+bx+cx^2+...
}

# Given a set of complex roots, this function groups complex roots into conjugate pairs
get_conjugate_pairs <- function(roots) {
  used <- rep(FALSE, length(roots))
  pairs <- list()
  for (i in seq_along(roots)) {
    if (!used[i]) {
      conj_root <- Conj(roots[i])
      # Find its conjugate (within tolerance)
      j <- which(!used & abs(Re(roots) - Re(conj_root)) < 1e-8 & abs(Im(roots) - Im(conj_root)) < 1e-8)
      j <- j[j != i]
      if (length(j) > 0) {
        used[c(i, j[1])] <- TRUE
        pairs[[length(pairs) + 1]] <- c(roots[i], roots[j[1]])
      }
    }
  }
  return(pairs)
}

# This function builds the polynomial c_m\b_{m+1} prod_{i=1}^m (1-r_{1i}x) + tau \prod_{j=1}^{m+1} (1-r_{2j}x)
# For example compute_sum_poly(0, rep(0,2), c(2,-3,1),1) = c(6,-7,0,1) because (1-2x)(1+3x)(1-x) = 6x^3 - 7x^2 + 0x + 1
compute_sum_poly <- function(factor, pr_roots, pl_roots, cte) {
  pr_coef <- c(0, poly_from_roots(pr_roots)) # in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
  pl_coef <- poly_from_roots(pl_roots) # in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
  pr_plus_pl_coef <- factor * pr_coef + cte * pl_coef
  return(pr_plus_pl_coef) # in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
}

compute_real_roots_and_complex_coef <- function(pr_plus_pl_coef) {
  pr_plus_pl_roots <- polyroot(rev(pr_plus_pl_coef)) # polyroot expects the coefficients in increasing order
  real_roots <- Re(pr_plus_pl_roots[abs(Im(pr_plus_pl_roots)) < 1e-8])
  complex_roots <- pr_plus_pl_roots[abs(Im(pr_plus_pl_roots)) >= 1e-8]
  complex_poly_coefs <- list()
  if (length(complex_roots) > 0) {
    pairs <- get_conjugate_pairs(complex_roots)
    for (pair in pairs) {
      coef <- rev(Re(poly_from_roots(c(pair[1], pair[2]))))  # returns in x^2 + bx + c order
      complex_poly_coefs[[length(complex_poly_coefs) + 1]] <- round(coef, 12)
    }
  } else {
    complex_poly_coefs <- list()  # No complex roots
  }
  return(list(new_factor = pr_plus_pl_coef[1], 
              pr_plus_pl_roots = pr_plus_pl_roots, 
              real_roots = real_roots, 
              complex_poly_coefs = complex_poly_coefs))
}


m1table <- rSPDE:::m1table
m2table <- rSPDE:::m2table
m3table <- rSPDE:::m3table
m4table <- rSPDE:::m4table

v <- c()
ms <- c()
alphas <- c()
time_steps <- c()
for (time_step in rev(c(0.0001, 0.001, 0.01, 0.1,1))) {
  for (beta in m1table$beta[1:86]) {
    for (m in c(1,2,3,4)) {
      roots <- my.get.roots(m, beta)
      v <- c(length(compute_real_roots_and_complex_coef(compute_sum_poly(roots$factor, roots$rc, roots$rb, time_step))$complex_poly_coefs),v)
      ms <- c(m,ms)
      alphas <- c(2*beta,alphas)
      time_steps <- c(time_step,time_steps)
    }
  }
}

sum(v)

df <- data.frame(ms, alphas, time_steps, v)





