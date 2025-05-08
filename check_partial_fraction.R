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

m1table <- rSPDE:::m1table
m2table <- rSPDE:::m2table
m3table <- rSPDE:::m3table
m4table <- rSPDE:::m4table

poly_from_roots <- function(roots) {
  coef <- 1
  for (r in roots) {
    coef <- convolve(coef, c(1, -r), type = "open")
  }
  return(coef) # returns in increasing order like a+bx+cx^2+...
}

compute_partial_fraction_param <- function(factor, pr_roots, pl_roots, cte) {
  pr_coef <- c(0, poly_from_roots(pr_roots)) # in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
  pl_coef <- poly_from_roots(pl_roots) # in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
  factor_pr_coef <- factor * pr_coef
  pr_plus_pl_coef <- factor_pr_coef + cte * pl_coef
  res <- gsignal::residue(factor_pr_coef, pr_plus_pl_coef)
  return(list(factor_pr_coef = factor_pr_coef, pr_plus_pl_coef = pr_plus_pl_coef, r = res$r, poles = res$p, k = res$k)) # in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
}

alpha <- 1 # from 0.5 to 2
m = 4
beta <- alpha/2
roots <- my.get.roots(m, beta)
Pl.roots <- roots$rb
Pr.roots <- roots$rc
factor <- roots$factor

res <- compute_partial_fraction_param(factor, Pr.roots, Pl.roots, 1)
f_x <- function(x) {
  up <- x^4*res$factor_pr_coef[2] + x^3*res$factor_pr_coef[3] + x^2*res$factor_pr_coef[4] + x*res$factor_pr_coef[5] + res$factor_pr_coef[6]
  down <- x^5*res$pr_plus_pl_coef[1] + x^4*res$pr_plus_pl_coef[2] + x^3*res$pr_plus_pl_coef[3] + x^2*res$pr_plus_pl_coef[4] + x*res$pr_plus_pl_coef[5] + res$pr_plus_pl_coef[6]
  return(up/down)
}

g_x <- function(x){
  res$r[1]/(x - res$poles[1]) + res$r[2]/(x - res$poles[2]) + res$r[3]/(x - res$poles[3]) + res$r[4]/(x - res$poles[4]) + res$r[5]/(x - res$poles[5]) + ifelse(is.null(res$k), 0, res$k)
}

x <- seq(-10, 10, by = 0.1)
g_x(x) 


sum((f_x(x) - g_x(x))^2)  # Check if the two polynomials are equal
plot(x, f_x(x), type = "l", col = "blue", main = "Polynomials Comparison", ylab = "f(x) and g(x)", xlab = "x")
lines(x, g_x(x), col = "red")
lines(x, f_x(x)-g_x(x)+1, col = "black")  # Add x-axis



# library(gsignal)
# 
# # Define transfer function: 1 / (x-2)(x+3)
# b <- c(0, 0, 1)
# a <- c(1, 1, -6)
# 
# # computes b/a
# res <- gsignal::residue(b, a) # a and b in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
#   
# print(res$r)  # residues
# print(res$p)  # poles
# print(res$k)  # direct term
