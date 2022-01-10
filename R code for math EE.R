library(statnet)
library(degreenet)
library(ergm)
library(igraph)
library(statnet.common)
library(ergm.rank)
library(latentnet)
library(relevent)
library(ergm.userterms)
library(networksis)
library(EpiModel)
library(coda)
library(intergraph)
library(network)
library(remotes)
library(ggplot2)
library(sna)

set.seed(1)
data("faux.mesa.high")
data("faux.magnolia.high")
data("faux.desert.high")
data("faux.dixon.high")
mesa <- faux.mesa.high
magnolia <- faux.magnolia.high
desert <- faux.desert.high
dixon <- faux.dixon.high

#Convert to graph project
ms <- intergraph::asIgraph(mesa)
mg <- intergraph::asIgraph(magnolia)
ds <- intergraph::asIgraph(desert)
dx <- intergraph::asIgraph(dixon)


ms.b <- intergraph::asIgraph(mesa.b.sim)
ms.e1 <- intergraph::asIgraph(mesa.e1.sim)
ms.e2 <- intergraph::asIgraph(mesa.e2.sim)
ms.e3 <- intergraph::asIgraph(mesa.e3.sim)

#Plotting
plot(mesa, main = "Faux Mesa High", vertex.cex = 0.7)
plot(magnolia, main = "Faux Magnolia High", vertex.cex = 0.7)
plot(desert, main = "Faux Desert High", vertex.cex = 0.7)
plot(dixon, main = "Faux Dixon High", vertex.cex = 0.7)

par(mfrow=c(1,1))
plot(mesa, vertex.col='Grade', main = "Faux Mesa High", vertex.cex = 0.7)
legend('bottomleft',fill=7:12,
       legend=paste('Grade',7:12),cex=0.75)

plot(mesa, vertex.col='Race', main = "Faux Mesa High", vertex.cex = 0.7)
legend('bottomleft',fill=1:5,
       legend=c("Black", "Hisp", "NatAm", "Other", "White"),cex=0.75)

plot(mesa, vertex.col='Sex', main = "Faux Mesa High")
legend('bottomleft',fill=1:5,
       legend=c("Male", "Female"),cex=0.75)

plot(dixon, displayisolates = FALSE, vertex.cex = 0.7, main = "Faux Dixon High")

#Degree distribution
dd.ms <- degree_distribution(ms)
dd.ms.ed <- degree_distribution(mesa.ed)
dd.ms.b <- degree_distribution(ms.b)
dd.ms.e1 <- degree_distribution(mesa.ed)
dd.ms.e2 <- degree_distribution(mesa.ed)
dd.ms.e3 <-degree_distribution(mesa.ed)

#Erdos-Renyi

set.seed(1)
mesa.ed <- sample_gnm(205,203)



#Bernoulli

set.seed(1)
mesa.b <- ergm(mesa ~ edges)
summary(mesa.b)  #Finding coefficients


#Examining the quality of model fit 
mesa.b.gof <- gof(mesa.b)
mesa.b.gof
par(mfrow=c(2,2))
plot(mesa.b.gof)


#Calculation & summary of structural effects of observed and simulated networks of Mesa High
set.seed(1)
mesa.b.sim <- simulate((mesa.b))
summary(mesa ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))
summary(mesa.b.sim ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))


par(mfrow=c(1,2))
plot(mesa, vertex.col='Grade', main = "real mesa")
plot(simulate(fauxmodel.01), vertex.col='Grade', main = "simulated mesa")
legend('bottomleft',fill=7:12,
       legend=paste('Grade',7:12),cex=0.75)

#ERGM 1

set.seed(1)
mesa.e1 <- ergm(mesa ~ edges + isolates)
summary(mesa.e1)

#Examining the quality of model fit 
mesa.e1.gof <- gof(mesa.e1)
mesa.e1.gof
par(mfrow=c(2,2))
plot(mesa.e1.gof)

set.seed(1)
mesa.e1.sim <- simulate(mesa.e1)
summary(mesa.e1.sim ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))

#ERGM 2

set.seed(1)
mesa.e2 <- ergm(mesa ~ edges + isolates + nodefactor("Grade") + nodefactor("Race") +
                  nodefactor("Sex") + nodematch("Grade",diff=TRUE) +
                  nodematch("Race",diff=TRUE) + nodematch("Sex",diff=FALSE))
summary(mesa.e2)

mesa.e2.gof <- gof(mesa.e2)
mesa.e2.gof
par(mfrow=c(2,2))
plot(mesa.e2.gof)


set.seed(1)
mesa.e2.sim <- simulate(mesa.e2)
summary(mesa.e2.sim ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))

#ERGM 3 

set.seed(1)
mesa.e3 <- ergm(mesa ~ edges + nodefactor("Grade") + nodefactor("Race") +
                  nodefactor("Sex") + nodematch("Grade",diff=TRUE) +
                  nodematch("Race",diff=TRUE) + nodematch("Sex",diff=FALSE) + 
                  gwdegree(1.0,fixed=TRUE) + gwesp(1.0,fixed=TRUE) +  gwdsp(1.0,fixed=TRUE))
summary(mesa.e3)
mesa.e3.gof <- gof(mesa.e3)
mesa.e3.gof
par(mfrow=c(2,2))
plot(mesa.e3.gof)

set.seed(1)
mesa.e3.sim <- simulate(mesa.e3)
summary(mesa.e3.sim ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))


#Degree Distribution comparison

plot(dd.ms, type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Frequency")
lines(dd.ms.b, lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

#b
plot(summary(mesa ~ degree(0:10)), type = "l", lty = 1, lwd = 2,
            xlab = "Degree", ylab = "Count")
lines(summary(mesa.b.sim ~ degree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
            lty = 1:2)

#e1
plot(summary(mesa ~ degree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(mesa.e1.sim ~ degree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

#e2
plot(summary(mesa ~ degree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(mesa.e2.sim ~ degree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)


################################################################################

#MAGNOLIA

#Convert to graph project

mg.b <- intergraph::asIgraph(mag.b.sim)
mg.e1 <- intergraph::asIgraph(mag.e1.sim)
mg.e2 <- intergraph::asIgraph(mag.e2.sim)
mg.e3 <- intergraph::asIgraph(mag.e3.sim)

#Plotting

plot(magnolia)


#Degree distribution
dd.mg <- degree
dd.mg.ed <-
dd.mg.b <-
dd.mg.e1 <- 
dd.mg.e2 <-
dd.mg.e3 <-
  
#Erdos-Renyi
  
set.seed(1)
mag.ed <- sample_gnm(1461,974)

summary(magnolia ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))


#Bernoulli

set.seed(1)
mag.b <- ergm(magnolia ~ edges)
summary(mag.b)  #Finding coefficients

#Examining the quality of model fit 
mag.b.gof <- gof(mag.b)
mag.b.gof
par(mfrow=c(2,2))
plot(mag.b.gof)


#Calculation & summary of structural effects of observed and simulated networks of Mesa High
set.seed(1)
mag.b.sim <- simulate(mag.b)
summary(mag.b.sim ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))


par(mfrow=c(1,2))
plot(mesa, vertex.col='Grade', main = "real mesa")
plot(simulate(), vertex.col='Grade', main = "simulated mesa")
legend('bottomleft',fill=7:12,
       legend=paste('Grade',7:12),cex=0.75)


#ERGM fit
set.seed(1)
mag.e1 <- ergm(magnolia ~ edges + isolates)
summary(mesa.e1)

#Examining the quality of model fit 
mag.e1.gof <- gof(mag.e1)
mag.e1.gof
par(mfrow=c(2,2))
plot(mag.e1.gof)

set.seed(1)
mag.e1.sim <- simulate(mag.e1)
summary(mag.e1.sim ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))

set.seed(1)
mag.e2 <- ergm(magnolia ~ edges + isolates + nodefactor("Grade") + nodefactor("Race") +
                  nodefactor("Sex") + nodematch("Grade",diff=TRUE) +
                  nodematch("Race",diff=TRUE) + nodematch("Sex",diff=FALSE))
summary(mag.e2)

mag.e2.gof <- gof(mag.e2)
mag.e2.gof
par(mfrow=c(2,2))
plot(mag.e2.gof)


set.seed(1)
mag.e2.sim <- simulate(mag.e2)
summary(mag.e2.sim ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))


set.seed(1)
mag.e3 <- ergm(magnolia ~ edges + nodefactor("Grade") + nodefactor("Race") +
                  nodefactor("Sex") + nodematch("Grade",diff=TRUE) +
                  nodematch("Race",diff=TRUE) + nodematch("Sex",diff=FALSE) + 
                  gwdegree(1.0,fixed=TRUE) + gwesp(1.0,fixed=TRUE) + gwdsp(1.0,fixed=TRUE))

mag.e3.gof <- gof(mag.e3)
mag.e3.gof
par(mfrow=c(2,2))
plot(mag.e3.gof)

set.seed(1)
mag.e3.sim <- simulate(mag.e3)
summary(mag.e3.sim ~ edges + triangles +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))

#Degree Distribution


#b
plot(summary(magnolia ~ degree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(mag.b.sim ~ degree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

#e1
plot(summary(magnolia ~ degree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(mag.e1.sim ~ degree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

#e2
plot(summary(magnolia ~ degree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(mag.e2.sim ~ degree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)
################################################################################

#DESERT

#Convert to graph project

ds.b <- intergraph::asIgraph(mesa.b.sim)
ds.e1 <- intergraph::asIgraph(mesa.e1.sim)
ds.e2 <- intergraph::asIgraph(mesa.e2.sim)
ds.e3 <- intergraph::asIgraph(mesa.e3.sim)

#Plotting

plot(desert)



#Degree distribution
dd.ds <- degree
dd.ds.ed <-
dd.ds.b <-
dd.ds.e1 <- 
dd.ds.e2 <-
dd.ds.e3 <-
  
  #Erdos-Renyi
  
set.seed(1)
des.ed <- sample_gnm(107,439, direct = TRUE)


#Bernoulli

set.seed(1)
des.b <- ergm(desert ~ edges)
summary(des.b)  #Finding coefficients

#Examining the quality of model fit 
des.b.gof <- gof(des.b)
des.b.gof
par(mfrow=c(2,2))
plot(des.b.gof)


#Calculation & summary of structural effects of observed and simulated networks of Mesa High
set.seed(1)
des.b.sim <- simulate((des.b))
summary(desert ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))
summary(des.b.sim ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))


par(mfrow=c(1,2))
plot(des, vertex.col='Grade', main = "real mesa")
plot(simulate(), vertex.col='Grade', main = "simulated mesa")
legend('bottomleft',fill=7:12,
       legend=paste('Grade',7:12),cex=0.75)



set.seed(1)
des.e1 <- ergm(desert ~ edges + mutual)
summary(des.e1)


#Examining the quality of model fit 
des.e1.gof <- gof(des.e1)
des.e1.gof
par(mfrow=c(2,2))
plot(des.e1.gof)

set.seed(1)
des.e1.sim <- simulate(des.e1)
summary(des.e1.sim ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))

set.seed(1)
des.e2 <- ergm(desert ~ edges + mutual + nodefactor("grade") + nodefactor("race") +
                  nodefactor("sex") + nodematch("grade",diff=TRUE) +
                  nodematch("race",diff=TRUE) + nodematch("sex",diff=FALSE))
summary(des.e2)

des.e2.gof <- gof(des.e2)
des.e2.gof
par(mfrow=c(2,2))
plot(des.e2.gof)


set.seed(1)
des.e2.sim <- simulate(des.e2)
summary(des.e2.sim ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))


set.seed(1)
des.e3 <- ergm(desert ~ edges + nodefactor("Grade") + nodefactor("Race") +
                  nodefactor("Sex") + nodematch("Grade",diff=TRUE) +
                  nodematch("Race",diff=TRUE) + nodematch("Sex",diff=FALSE) + 
                  gwdegree(1.0,fixed=TRUE) + gwesp(1.0,fixed=TRUE) + gwdsp(1.0,fixed=TRUE))

des.e3.gof <- gof(des.e3)
des.e3.gof
par(mfrow=c(2,2))
plot(des.e3.gof)

set.seed(1)
des.e3.sim <- simulate(des.e3)
summary(mes.e3.sim ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('Sex') + nodematch('Grade') + nodematch('Race'))

#Degree Distribution

#b
plot(summary(desert ~ idegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(des.b.sim ~ idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

plot(summary(desert ~ odegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(des.b.sim ~ odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

#e1
plot(summary(desert ~ idegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(des.e1.sim ~ idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

plot(summary(desert ~ odegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(des.e1.sim ~ odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

#e2
plot(summary(desert ~ idegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(des.e2.sim ~ idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

plot(summary(desert ~ odegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(des.e2.sim ~ odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

################################################################################

#DIXON

#Convert to graph project

dx.b <- intergraph::asIgraph(dix.b.sim)
dx.e1 <- intergraph::asIgraph(dix.e1.sim)
dx.e2 <- intergraph::asIgraph(dix.e2.sim)
dx.e3 <- intergraph::asIgraph(dix.e3.sim)

#Plotting

plot(dixon)



#Degree distribution
dd.dx <- degree
dd.dx.ed <-
  dd.dx.b <-
  dd.dx.e1 <- 
  dd.dx.e2 <-
  dd.dx.e3 <-
  
#Erdos-Renyi
  
set.seed(1)
dix.ed <- sample_gnm(248, 1197, direct = TRUE)



#Bernoulli

set.seed(1)
dix.b <- ergm(dixon ~ edges)
summary(dix.b)  #Finding coefficients

#Examining the quality of model fit 
dix.b.gof <- gof(dix.b)
dix.b.gof
par(mfrow=c(2,2))
plot(dix.b.gof)


#Calculation & summary of structural effects of observed and simulated networks of Mesa High
set.seed(1)
dix.b.sim <- simulate(dix.b)
summary(dixon ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))
summary(dix.b.sim ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))


par(mfrow=c(1,2))
plot(dixon, vertex.col='grade', main = "real dix")
plot(simulate(), vertex.col='Grade', main = "simulated mesa")
legend('bottomleft',fill=7:12,
       legend=paste('Grade',7:12),cex=0.75)



set.seed(1)
dix.e1 <- ergm(dixon ~ edges + mutual)
summary(dix.e1)


#Examining the quality of model fit 
dix.e1.gof <- gof(dix.e1)
dix.e1.gof
par(mfrow=c(2,2))
plot(dix.e1.gof)

set.seed(1)
dix.e1.sim <- simulate(dix.e1)
summary(dix.e1.sim ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))

set.seed(1)
dix.e2 <- ergm(dixon ~ edges + mutual + nodefactor("grade") + nodefactor("race") +
                 nodefactor("sex") + nodematch("grade",diff=TRUE) +
                 nodematch("race",diff=TRUE) + nodematch("sex",diff=FALSE))
summary(dix.e2)

dix.e2.gof <- gof(dix.e2)
dix.e2.gof
par(mfrow=c(2,2))
plot(dix.e2.gof)


set.seed(1)
dix.e2.sim <- simulate(dix.e2)
summary(dix.e2.sim ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))


set.seed(1)
dix.e3 <- ergm(dixon ~ edges + isolates + mutual +nnodefactor("grade") + nodefactor("race") +
                 nodefactor("sex") + nodematch("grade",diff=TRUE) +
                 nodematch("race",diff=TRUE) + nodematch("Sex",diff=FALSE) + 
                 gwdegree(1.0,fixed=TRUE) + gwesp(1.0,fixed=TRUE) + gwdsp(1.0,fixed=TRUE))

dix.e3.gof <- gof(dix.e3)
dix.e3.gof
par(mfrow=c(2,2))
plot(dix.e3.gof)

set.seed(1)
dix.e3.sim <- simulate(dix.e3)
summary(dix.e3.sim ~ edges + triangles + mutual +
          density + isolates + meandeg + 
          nodematch('sex') + nodematch('grade') + nodematch('race'))

#Degree Distribution

#b
plot(summary(dixon ~ idegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(dix.b.sim ~ idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

plot(summary(dixon ~ odegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(dix.b.sim ~ odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

#e1

plot(summary(dixon ~ idegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(dix.e1.sim ~ idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

plot(summary(dixon ~ odegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(dix.e1.sim ~ odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)



#e2

plot(summary(dixon ~ idegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(dix.e2.sim ~ idegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
       lty = 1:2)

plot(summary(dixon ~ odegree(0:10)), type = "l", lty = 1, lwd = 2,
     xlab = "Degree", ylab = "Count")
lines(summary(dix.e2.sim ~ odegree(0:10)), lty = 2, lwd = 3)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3,
      lty = 1:2)





