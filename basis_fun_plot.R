library(MetricGraph)
#F
edge1 <- rbind(c(0,0), c(0,2))
edge2 <- rbind(c(0,2), c(0,4))
edge3 <- rbind(c(0,2), c(1,2))
edge4 <- rbind(c(0,4), c(2,4))
#U
edge5 <- rbind(c(2,4), c(2,1))
thetau <- seq(pi, 2*pi, length.out = 100)
edge6 <- cbind(3+1*cos(thetau), 1+1*sin(thetau))
edge7 <- rbind(c(4,1), c(4,4))
#N
edge8 <- rbind(c(4,1), c(4,0))
edge9 <- rbind(c(4,4), c(6,0))
edge10 <- rbind(c(6,0), c(6,4))
#C
thetac1 <- seq(pi, pi/2, length.out = 100)
edge11 <- cbind(8+2*cos(thetac1), 2+2*sin(thetac1))
thetac2 <- seq(pi, 3*pi/2, length.out = 100)
edge12 <- cbind(8+2*cos(thetac2), 2+2*sin(thetac2))
#T
edge13 <- rbind(c(8,4), c(10,4))
edge14 <- rbind(c(9,4), c(9,0))
#I
edge15 <- rbind(c(10,4), c(12,4))
edge16 <- rbind(c(10,0), c(12,0))
edge17 <- rbind(c(11,0), c(11,4))
#O
thetao1 <- seq(pi, 2*pi, length.out = 100)
edge18 <- cbind(13+1*cos(thetao1), 1+1*sin(thetao1))
thetao2 <- seq(0, pi, length.out = 100)
edge19 <- cbind(13+1*cos(thetao2), 3+1*sin(thetao2))
edge20 <- rbind(c(12,3), c(12,1))
edge21 <- rbind(c(14,1), c(14,4))
#N
edge22 <- rbind(c(14,4), c(16,0))
edge23 <- rbind(c(16,0), c(16,4))
edge56 <- rbind(c(14,0), c(14,1))
#S
edge24 <- rbind(c(16,0), c(17,0))
tethas1 <- seq(-pi/2, pi/2, length.out = 100)
edge25 <- cbind(17+1*cos(tethas1), 1+1*sin(tethas1))
thetas2 <- seq(3*pi/2, pi/2, length.out = 100)
edge26 <- cbind(17+1*cos(thetas2), 3+1*sin(thetas2))
edge27 <- rbind(c(17,4), c(18,4))
#H
edge28 <- rbind(c(0,4), c(0,6))
edge29 <- rbind(c(0,6), c(0,8))
edge30 <- rbind(c(0,6), c(2,6))
edge31 <- rbind(c(2,4), c(2,8))
#A
edge32 <- rbind(c(2,4), c(3,8))
edge33 <- rbind(c(3,8), c(4,4))
edge34 <- rbind(c(2.5,6), c(3.5,6))
#T
edge35 <- rbind(c(3,8), c(6,8))
edge36 <- rbind(c(5,8), c(5,4))
#B
edge37 <- rbind(c(8,8), c(8,4))
thetab1 <- seq(-pi/2, pi/2, length.out = 100)
edge38 <- cbind(9+1*cos(thetab1), 5+1*sin(thetab1))
edge39 <- rbind(c(8,6), c(9,6))
edge40 <- rbind(c(9,8), c(8,8))
thetab2 <- seq(-pi/2, pi/2, length.out = 100)
edge41 <- cbind(9+1*cos(thetab2), 7+1*sin(thetab2))
#A
edge42 <- rbind(c(10,4), c(11,8))
edge43 <- rbind(c(11,8), c(12,4))
edge44 <- rbind(c(10.5,6), c(11.5,6))
#S

#I
edge45 <- rbind(c(14,4), c(16,4))
edge46 <- rbind(c(15,4), c(15,8))
edge47 <- rbind(c(14,8), c(16,8))
#S
edge48 <- rbind(c(16,4), c(17,4))
edge49 <- cbind(17+1*cos(tethas1), 5+1*sin(tethas1))
edge50 <- cbind(17+1*cos(thetas2), 7+1*sin(thetas2))
edge51 <- rbind(c(17,8), c(18,8))
#S
edge52 <- rbind(c(13,8), c(14,8))
edge53 <-cbind(13+1*cos(thetas2), 7+1*sin(thetas2))
edge54 <- cbind(13+1*cos(tethas1), 5+1*sin(tethas1))
edge55 <- rbind(c(12,4), c(13,4))

edges <- list(edge1, edge2, edge3, edge4, edge5, edge6, edge7,
              edge8, edge9, edge10, edge11, edge12, edge13, edge14,
              edge15, edge16, edge17, edge18, edge19, edge20, edge21,
              edge22, edge23, edge24, edge25, edge26, edge27,
              edge28, edge29, edge30, edge31, edge32, edge33, edge34,
              edge35, edge36, edge37, edge38, edge39, edge40, edge41,
              edge42, edge43, edge44, edge45, edge46, edge47,
              edge48, edge49, edge50, edge51, edge52, edge53, edge54, edge55, edge56)
graph <- metric_graph$new(edges = edges, perform_merges = TRUE)

graph$plot(direction = TRUE)
