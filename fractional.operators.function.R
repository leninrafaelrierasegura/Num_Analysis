fractional.operators <- function(L, # kappa^2C + G
                                 beta,
                                 C,
                                 scale.factor, # kappa^2
                                 m = 1,
                                 tau = 1) {

  
  C <- Matrix::Diagonal(dim(C)[1], rowSums(C)) # lumped
  Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C)) # lumped 
  I <- Matrix::Diagonal(dim(C)[1])
  L <- L / scale.factor # C + G/kappa^2
  CiL <- Ci %*% L
  LCi <- L %*% Ci
  roots <- get.roots(m, beta)
  Pl.roots <- roots$rb
  Pr.roots <- roots$rc
    # construct Pl
    Pl <- I - CiL * roots$rb[1]
    if (length(roots$rb) > 1) {
      for (i in 2:length(roots$rb)) {
        Pl <- Pl %*% (I - CiL * roots$rb[i])
      }
    }
    # So far Pl = \prod_{j=1}^{m+1} (I - CiL * rb[j])
    Pl.scaling <- scale.factor^beta / roots$factor # (kappa^2)^beta / factor
    Pl <- C %*% Pl # Pl becomes C %*% \prod_{j=1}^{m+1} (I - CiL * rb[j])
    
    # construct Pr
    Pr <- I - CiL * roots$rc[1]
    if (length(roots$rc) > 1) {
      for (i in 2:length(roots$rc)) {
        Pr <- Pr %*% (I - CiL * roots$rc[i])
      }
    }
    # So far Pr = \prod_{i=1}^{m} (I - CiL * rc[j])
  
  Pl <- Pl * Pl.scaling # Pl now is (kappa^2)^beta / factor * C %*% \prod_{j=1}^{m+1} (I - CiL * rb[j])
  output <- list(
    Pl = Pl,
    Pr = Pr,
    Ci = Ci,
    C = C,
    CiL = CiL,
    L = L,
    m = m,
    beta = beta,
    type = "fractional approximation"
  )
  class(output) <- "rSPDEobj"
  return(output)
}