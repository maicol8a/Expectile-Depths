# Expectile-Depths
Expectile-Depths is a collection of R functions designed for the computation of bivariate expectile depth and regions, as well as for plotting BExPlots. These tools are invaluable for statistical analysis, enabling users to visualize and compare the depth of expectiles in bivariate data sets.

# Features
Expectile Depth Calculation: Algorithms to compute the depth of expectiles in bivariate data, helping identify patterns and anomalies.
BExPlot Generation: Functions to create BExPlots, which are visual representations of expectile depth, aiding in data distribution analysis.
Expectile Regions: Tools to plot contours of expectile regions, highlighting areas of interest within the data.
D-D Plots: Comparative depth-depth plots to assess similarities or differences between two data sets.
Example Usage
The repository includes examples demonstrating how to use the provided functions. Below is a sample workflow:

# Example 1
## Load sample data
data(hemophilia)

attach(hemophilia)

AHF.normal <- cbind(AHFactivity[gr=="normal"], AHFactivity.1[gr=="normal"])

AHF.carrier <- cbind(AHFactivity[gr=="carrier"], AHFactivity.1[gr=="carrier"])

## Generate BExPlots
par(mfrow=c(3,2))

BExPlot(AHF.normal)

BExPlot(AHF.carrier)

## Plot expectile regions
plot(AHF.carrier, pch=3)

points(x=mean(AHF.carrier[,1]), y=mean(AHF.carrier[,2]), pch=16)

lines(exactexp(AHF.carrier))

lines(exactexp(AHF.carrier, alpha=0.01))

lines(exactexp(AHF.carrier, alpha=0.05))

lines(exactexp(AHF.carrier, alpha=0.4))

# Example 2
## Generate D-D plot
X <- matrix(rnorm(100), ncol=2)

X <- rbind(X, matrix(rnorm(4, mean=5), ncol=2))

Y <- matrix(rnorm(100), ncol=2)

Y <- rbind(Y, matrix(rnorm(4, mean=-5), ncol=2))

sim.data <- rbind(X, Y)

e1 <- vector(length=nrow(sim.data))

for(i in 1:length(e1)) e1[i] <- expdepth(x=sim.data[i,], data=X)

e2 <- vector(length=nrow(sim.data))

for(i in 1:length(e2)) e2[i] <- expdepth(x=sim.data[i,], data=Y)

plot(x=e1, y=e2, main="D-D plot, expectile depth", xlab="X", ylab="Y")

abline(0, 1, lty=2)

# Authors
The code was created by I. Cascos (https://github.com/icascos) and M. Ochoa. The theoretical description of the procedures is available at [e-archivo.uc3m.es.](https://e-archivo.uc3m.es/handle/10016/28434)
