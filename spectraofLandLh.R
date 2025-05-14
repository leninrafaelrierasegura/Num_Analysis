

gets.graph.tadpole <- function(h){
  edge1 <- rbind(c(0,0),c(1,0))
  theta <- seq(from=-pi,to=pi,length.out = 10000)
  edge2 <- cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi)
  edges <- list(edge1, edge2)
  graph <- metric_graph$new(edges = edges)
  graph$set_manual_edge_lengths(edge_lengths = c(1,2))
  graph$build_mesh(h = h)
  return(graph)
}

kappa <- 100

h = 0.01
graph <- gets.graph.tadpole(h = h)
graph$compute_fem()
G <- graph$mesh$G
C <- graph$mesh$C
Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C)) 
L <- (kappa^2*C + G)%*%Ci
scaled_L <- L/kappa^2

val <- eigen(L)$val
scaled_val <- eigen(scaled_L)$val

c(0, kappa^-2)*kappa^2
c(1/max(val), 1/min(val))*kappa^2
range(Re(eigen(solve(L))$val))*kappa^2

c(0, kappa^-2)
c(1/max(scaled_val), 1/min(scaled_val))
range(Re(eigen(solve(scaled_L))$val))
