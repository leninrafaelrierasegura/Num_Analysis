poly_from_roots <- function(roots) {
  coef <- 1
  for (r in roots) {
    coef <- convolve(coef, c(1, -r), type = "open")
  }
  return(coef) # returns in increasing order like a+bx+cx^2+...
}
 
# Group complex roots into conjugate pairs
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

compute_sum_poly <- function(factor, pr_roots, pl_roots, cte) {
  
  pr_coef <- c(0, poly_from_roots(pr_roots)) # in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
  pl_coef <- poly_from_roots(pl_roots) # in decreasing order like x^n+bx^(n-1)+cx^(n-2)+...
  pr_plus_pl_coef <- factor * pr_coef + cte * pl_coef
  return(pr_plus_pl_coef)
}

compute_real_roots_and_complex_coef <- function(pr_plus_pl_coef) {

  pr_plus_pl_roots <- polyroot(rev(pr_plus_pl_coef))
  
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

factor <- 2
pr_roots <- c(1, 2, -9)
pl_roots <- c(3, 4, -5, 7)
cte <- -5
result <- compute_real_roots_and_complex_coef(compute_sum_poly(factor, pr_roots, pl_roots, cte))

f_x <- function(x) {
  factor * (1 - pr_roots[1]*x) * (1 - pr_roots[2]*x) * (1 - pr_roots[3]*x) + cte * (1 - pl_roots[1]*x) * (1 - pl_roots[2]*x) * (1 - pl_roots[3]*x) * (1 - pl_roots[4]*x)
}

g_x <- function(x) {
  result$new_factor * (x - result$real_roots[1]) * (x - result$real_roots[2]) * (x^2 + result$complex_poly_coefs[[1]][2] * x + result$complex_poly_coefs[[1]][3])
}

x <- seq(-10, 10, by = 0.1)

sum((f_x(x) - g_x(x))^2)  # Check if the two polynomials are equal
plot(x, f_x(x), type = "l", col = "blue", main = "Polynomials Comparison", ylab = "f(x) and g(x)", xlab = "x")
lines(x, g_x(x), col = "red")
lines(x, f_x(x)-g_x(x), col = "black")  # Add x-axis




compute_real_roots_and_complex_coef(c(1,-1,3,-3,2,-2)) # (x^2+1)(x^2+2)(x-1)
compute_real_roots_and_complex_coef(c(1,-1,0,0,1,-1)) # (x^4+1)(x-1)
compute_real_roots_and_complex_coef(c(1,-6,5,-30)) # (x^2+5)(x-6)
a <- compute_real_roots_and_complex_coef(c(1,0,-7,6)) # (x-1)(x+3)(x-2)
a$complex_poly_coefs |> length()

for (i in 1:length(a$complex_poly_coefs)){
  print(i)
}
