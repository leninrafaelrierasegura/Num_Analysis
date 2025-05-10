
system.time({
  U_true <- EIGENFUN %*% outer(1:length(EIGENVAL_ALPHA), 
                               1:length(time_seq), 
                               function(i, j) coeff[i] * exp(-EIGENVAL_ALPHA[i] * time_seq[j]))
})
system.time({
  A <- EIGENFUN %*% sweep(exp(-outer(EIGENVAL_ALPHA, time_seq)), 1, coeff, "*")
})
system.time({
  A2 <- EIGENFUN %*% (coeff * exp(-outer(EIGENVAL_ALPHA, time_seq)))
})

sum(U_true- A)
sum(U_true- A2)
sum(A-A2)

U_0_true <- EIGENFUN %*% (coeff_U_0 * exp(-outer(EIGENVAL_ALPHA, time_seq))) # coeff_U_0 and EIGENVAL_ALPHA must have the same length
FF <- EIGENFUN %*% (coeff_FF * exp(-outer(EIGENVAL_ALPHA, time_seq)))
FF_true <- t(time_seq * t(FF))


FF_true <- matrix(NA, nrow = length(x), ncol = length(time_seq))
FF_sol_true <- matrix(NA, nrow = length(x), ncol = length(time_seq))
for (k in 1:length(time_seq)) {
  FF_true[, k] <- ff(time_seq[k]) # this is the right hand side function
  FF_sol_true[, k] <- time_seq[k]*FF_true[, k] # this is the second term in the solution
}

FF_approx <- matrix(0, nrow = length(x), ncol = length(time_seq))
for (i in 1:length(coeff_FF)) {
  FF_approx <- FF_approx + INT_BASIS_EIGEN[i,] %*% t(coeff_FF[i] * exp(-time_seq*EIGENVAL_ALPHA[i]))
}