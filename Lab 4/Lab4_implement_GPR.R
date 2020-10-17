###################################
########### GIVEN CODE ############
###################################

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

###################################
############ MY CODE ##############
###################################


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
  upperBand = means + sqrt(variances)*1.96
  lowerBand = means - sqrt(variances)*1.96
  
  plot(x, y, type="p", 
       xlim=c(min(XStar),
              max(XStar)),
       ylim=c(min(min(y), min(lowerBand)),
              max(max(y), max(upperBand))))
  lines(XStar, means, col = "red", lwd = 3)
  lines(XStar, upperBand, col = "blue", lwd = 2)
  lines(XStar, lowerBand, col = "blue", lwd = 2)
  points(x, y, col='green')
  
}

# TASK 2.1.2

# sigma_f: ????
sigmaf = 1
#The length scale, i.e how highly correlated function values should be for close input points (smoothing factor): 
ell = 0.3
# Observation:
x1 = 0.4
y1 = 0.719
# Evaluation points in interval [-1,1]:
XStar = seq(-1, 1, 0.01)

posterior1 = posteriorGP(x1, y1 , XStar, sigmaNoise = 0.1, SquaredExpKernel, sigmaf = sigmaf, ell = ell)
plotPosterior(posterior1$posteriorMean, posterior1$posteriorVariance, x1, y1, XStar)

# TASK 2.1.3

# Another observation: 
x2 = −0.6
y2 = −0.044
X = c(x1, x2)
Y = c(y1, y2)

posterior2 = posteriorGP(X, Y, XStar, sigmaNoise = 0.1, SquaredExpKernel, sigmaf = sigmaf, ell = ell)
plotPosterior(posterior2$posteriorMean, posterior2$posteriorVariance, X, Y, XStar)

# TASK 2.1.4

# Observations: 
X = c(-1.0, -0.6, -0.2, 0.4, 0.8)
Y = c(0.768, -0.044, -0.940, 0.719, -0.664)

posterior2 = posteriorGP(X, Y, XStar, sigmaNoise = 0.1, SquaredExpKernel, sigmaf = sigmaf, ell = ell)
plotPosterior(posterior2$posteriorMean, posterior2$posteriorVariance, X, Y, XStar)

# TASK 2.1.5

X = c(-1.0, -0.6, -0.2, 0.4, 0.8)
Y = c(0.768, -0.044, -0.940, 0.719, -0.664)
ell = 1

posterior2 = posteriorGP(X, Y, XStar, sigmaNoise = 0.1, SquaredExpKernel, sigmaf = sigmaf, ell = ell)
plotPosterior(posterior2$posteriorMean, posterior2$posteriorVariance, X, Y, XStar)
