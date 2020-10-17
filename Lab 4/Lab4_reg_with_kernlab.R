library(kernlab)

# Code from earlier part of lab
library("mvtnorm")

# Covariance function
SquaredExpKernel <- function(x1,x2,sigmaF=1,l=3){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

# Gaussian Regression Model
# INPUTS:
#   X: Vector of training inputs.
#   y: Vector of training targets/outputs.
#   XStar: Vector of inputs where the posterior distribution is evaluated, i.e. X∗.
#   sigmaNoise: Noise standard deviation σn.
#   k: Covariance function or kernel. That is, the kernel should be a separate function
posteriorGP <- function(X, y, XStar, sigmaNoise, k, sigmaf, ell){
  
  # Compute covariance matrix K(X,X)
  covM = k(X, X, sigmaf, ell)
  
  # Identity matrix 
  I = diag(dim(covM)[1])
  
  # Transposing the output from the cholesky-function since the implementation in R returns 
  # an upper triangular matrix and we want a lower triangular matrix. 
  L = t(chol(covM + (sigmaNoise**2)*I))
  
  # Compute covariance matrix K(X,XStar)
  KStar = k(X, XStar, sigmaf, ell)
  # Predictive mean
  alpha = solve(t(L), solve(L, y))
  FStarMean = t(KStar)%*%alpha
  
  # Compute covariance matrix K(XStar,XStar)
  KStar2 = k(XStar, XStar, sigmaf, ell)
  # Predictive variance
  v = solve(L, KStar)
  V = diag(KStar2 - t(v)%*%v)
  
  # Compute log marginal likelihood
  #n = length(y)
  #logSum = 0
  #for(i in (1:n)) {
  #  logSum = logSum + log(L[i,i])
  #}
  
  # logProbYX = 0.5*t(y)*alpha - logSum - n/2*log(2*pi) # Not needed in this lab? Is a part of algorithm in 
  # Rasmussen and Williams' book
  logProbYX = 1
  
  return(list('posteriorMean' = FStarMean, 'posteriorVariance' = V, 'logProb' = logProbYX))
  
}

# Function for plotting posterior distributions with 95% probability band.
# INPUTS:   
#   means: Computed means for a number of points
#   variances: Computed variances for a number of points 
#   x: x-values for observations
#   y: y-values for observations
#   XStar: Vector of inputs where the posterior distribution is evaluated, i.e. X∗.
plotPosterior <- function(means, variances, x, y, XStar){
  var = sd(means) * variances + means
  upperBand = means + sqrt(var)*1.96
  lowerBand = means - sqrt(var)*1.96
  
  plot(x, y, type="l",
       xlim = c(min(x), max(x)), 
       ylim = c((min(min(lowerBand), min(y))), max(max(upperBand), max(y))))
  lines(XStar, means, col = "red", lwd = 2)
  lines(XStar, upperBand, col = "blue", lwd = 2)
  lines(XStar, lowerBand, col = "blue", lwd = 2)
  
}

# New code

data = read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/
Code/TempTullinge.csv", header=TRUE, sep=";")

# Creating variables with index for every fifth reading
tempTime = 1:2190
time = tempTime[seq(1,2190,5)]
dayNr = 1:365
tempDay = c()
for(i in 1:6){
  tempDay = c(tempDay, dayNr)
}
day = tempDay[seq(1,2190,5)]
temperature = unlist(data['temp'])
tempSubset = temperature[time]

# TASK 2.2.1

# Returns a square exponential kernel function
# INPUTS:
#   sigmaf: standard deviation for f
#   ell: Smoothing factor
squareExpKernel <- function(sigmaf = 1, ell = 1) 
{
  rval <- function(x, y = NULL) {
    r = crossprod(x-y);
    #return((sigmaf**2)*exp(-r/2*(ell**2)))
    return(sigmaf^2 * exp(-0.5 * r / ell^2))
  }
  class(rval) <- "kernel"
  return(rval)
  
}

# Creating kernel with sigmaf = 1, ell = 1
kernel = squareExpKernel(1,1)
# Evaluating in point x = 1, x' = 2
covM = kernel(1,2)

# Computing kernel matrix for created kernel and vectors X, XStar
X = c(1, 3, 4)
XStar = c(2, 3 ,4)
kernelMatrix = kernelMatrix(kernel = kernel, X, XStar)

# TASK 2.2.2 & 2.2.3

#
lmFit = lm(tempSubset ~ time + time**2) # vad är detta exakt?
noise = sd(lmFit$residuals) 

# Creating kernel with sigmaf = 20, ell = 0.2
sigmaf = 20
ell = 0.2
kernel = squareExpKernel(sigmaf = sigmaf, ell = ell)

#
fitGP = gausspr(x = time, 
                y = tempSubset, 
                kernel = kernel, 
                var = noise**2)
meanPredict = predict(fitGP, time)

#
posterior = posteriorGP(X = scale(time), 
                        y = scale(tempSubset), 
                        XStar = scale(time), 
                        sigmaNoise = noise, 
                        k = SquaredExpKernel, 
                        sigmaf = sigmaf, 
                        ell = ell)

# Plotting results
plotPosterior(means = meanPredict, 
              variances = posterior$posteriorVariance,
              x = time,
              y = tempSubset,
              XStar = time)

# TASK 2.2.4

#
lmFit = lm(tempSubset ~ day + day**2) # vad är detta exakt?
noise = sd(lmFit$residuals) 

# Creating kernel with sigmaf = 20, ell = 0.2
sigmaf = 20
ell = 0.2
kernel = squareExpKernel(sigmaf = sigmaf, ell = ell)

#
fitGP = gausspr(x = day, 
                y = tempSubset, 
                kernel = kernel, 
                var = noise**2)
meanPredictDay = predict(fitGP, day)

#
posteriorDay = posteriorGP(X = scale(day), 
                        y = scale(tempSubset), 
                        XStar = scale(time), 
                        sigmaNoise = noise, 
                        k = SquaredExpKernel, 
                        sigmaf = sigmaf, 
                        ell = ell)

# Plotting results
plotPosterior(means = meanPredict, 
              variances = posteriorDay$posteriorVariance,
              x = time,
              y = tempSubset,
              XStar = time)

plot(time, tempSubset, type = "l")
lines(time, meanPredict, col = "red", lwd = 2)
lines(time, meanPredictDay, col = "blue", lwd = 2)

# TASK 2.2.5

periodicKernel <- function(sigmaf = 1, ell1 = 1, ell2 = 1, d) 
{
  rval <- function(x, y = NULL) {
    r = crossprod(x-y);
    exp1 = exp(-2*(sin(pi*sqrt(r)/d)**2)/ell1**2)
    exp2 = exp(-0.5*r*(ell2**2))
    return((sigmaf**2)*exp1*exp2)
  }
  class(rval) <- "kernel"
  return(rval)
  
}

lmFit = lm(tempSubset ~ time + time**2) # vad är detta exakt?
noise = sd(lmFit$residuals) 

# Creating kernel with sigmaf = 20, ell = 0.2
sigmaf = 20
ell1 = 1
ell2 = 10
d = sd(time)
kernel = periodicKernel(sigmaf = sigmaf, ell1 = ell1, ell2 = ell2, d)

#
fitGP = gausspr(x = time, 
                y = tempSubset, 
                kernel = kernel, 
                var = noise**2)
meanPredictPer = predict(fitGP, time)

#
#posteriorPer = posteriorGP1(X = scale(time), 
#                           y = scale(tempSubset), 
#                           XStar = scale(time), 
#                           sigmaNoise = noise, 
#                           k = periodicKernel, 
#                           sigmaf = sigmaf, 
#                           ell1 = ell1,
#                           ell2 = ell2)

# Plotting results
#plotPosterior(means = meanPredictPer, 
#              variances = posteriorPer$posteriorVariance,
#              x = time,
#              y = tempSubset,
#              XStar = time)

plot(time, tempSubset, type = "l")
lines(time, meanPredict, col = "red", lwd = 2)
lines(time, meanPredictDay, col = "blue", lwd = 2)
lines(time, meanPredictPer, col = "green", lwd = 2)
