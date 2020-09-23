library("HMM")
library("entropy")

set.seed(1)

# Number of readings
n = 100

### Question 1

# The states are numbered from 0-9 and a reading is anotated with a "'"
states = c("0","1","2","3","4","5","6","7","8","9")
observable_readings = c("0'","1'","2'","3'","4'","5'","6'","7'","8'","9'")

# Transmission model where the probability to stay in a given state is the same as 
# moving to next state.
transmission_model = t(matrix(
  c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0.5, 0.5, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5,
    0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0.5
  ),10,10))

# A emisson model where the probability of reading a certain state is uniformly distrubuted
# over the two previous states, the true state and the two following states. 
emission_model = t(matrix(
  c(0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.2, 0.2,
    0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.2,
    0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0,
    0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0,
    0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0,
    0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0,
    0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0,
    0, 0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2,
    0.2, 0, 0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2,
    0.2, 0.2, 0, 0, 0, 0, 0, 0.2, 0.2, 0.2
  ),10,10))

#Probabilities to start in each state
start_prob = rep(0.1,10)

hmm = initHMM(states, 
              observable_readings,
              startProbs = start_prob,
              transProbs = transmission_model, 
              emissionProbs = emission_model)
#print(hmm)

### Question 2

simed_hmm = simHMM(hmm, n)

print(simmed_hmm)

### Question 3

filtering <- function(hmm, observations) {
  
  alpha = exp(forward(hmm, observations))
  beta = exp(backward(hmm, observations))
  
  # Creating empty matrix to be filled in following loop
  filtered = matrix(, nrow = 10, ncol = length(observations))
  
  # Computing the filtered distribution
  for (i in 1:length(observations)) {
    filtered[,i] = alpha[,i]/sum(alpha[,i])
  }
  return(filtered)
}
  
smoothing <- function(hmm, observations) {
  
  alpha = exp(forward(hmm, observations))
  beta = exp(backward(hmm, observations))
  
  # Creating empty matrix to be filled in following loop
  smoothed = matrix(, nrow = 10, ncol = length(observations))
  
  # Computing the smoothed distribution 
  for (i in 1:length(observations)) {
    smoothed[,i] = alpha[,i]*beta[,i]/sum(alpha[,1:i]*beta[,1:i]) # Double check if this is correct?
  }
  
  
  return(smoothed)
}

observations = simed_hmm$observation

filtered_distribution = filtering(hmm, observations)
smoothed_distribution = smoothing(hmm, observations)

# Computing the most probable path with the viterbi algorithm from the HMM-package
most_prob_path = viterbi(hmm, observations)
  
### Question 4

predict <- function(dist, states) {

  predictions = c()
  
  for (i in 1:dim(dist)[2]) {
  
  # Predicts the state that has the highest probability, if there are multiple states with
  # equal probability, the state with the lowest number/index is predicted i.e if ["1","2"] is 
  # two states with equal probability in a certain point, "1" is predicted. 
  predictions[i] = states[which(dist[,i] == max(dist[,i]))]
  
  }
  
  return(predictions)
}

calc_acc <- function(hmm_in, pred) {
  
  # Fetching ture hidden states and creating vectors to store predictions in
  #true_states = simed_hmm$states
  true_states = hmm_in$states
  
  
  acc = sum(pred == true_states)/n
  
  return(acc)
}

# Computing the accuracy for the filtering method
filtered_pred = predict(filtered_distribution, states)
filtered_acc = calc_acc(simed_hmm, filtered_pred)

# Computing the accuracy for the smoothing method
smoothing_pred = predict(smoothed_distribution, states)
smoothing_acc = calc_acc(simed_hmm, smoothing_pred)

# Computing the accuracy for the smoothing method
mpp_acc = calc_acc(simed_hmm, most_prob_path)

### Question 5

# Repeating earlier calculations with different seed
set.seed(12345)

simed_hmm = simHMM(hmm, n)
observations = simed_hmm$observation

filtered_distribution = filtering(hmm, observations)
smoothed_distribution = smoothing(hmm, observations)
most_prob_path = viterbi(hmm, observations)

# Computing the accuracy for the filtering method
filtered_pred = predict(filtered_distribution, states)
filtered_acc_seed2 = calc_acc(simed_hmm, filtered_pred)

# Computing the accuracy for the smoothing method
smoothing_pred = predict(smoothed_distribution, states)
smoothing_acc_seed2 = calc_acc(simed_hmm, smoothing_pred)

# Computing the accuracy for the smoothing method
mpp_acc_seed2 = calc_acc(simed_hmm, most_prob_path)

### Question 6
entropies = c()

for (i in 1:n) {
  entropies[i] = entropy.empirical(filtered_distribution[,i])
  print(entropies[i])
}

plot(entropies)

### Question 7

next_predicted_dist = filtered_distribution[,n]%*%transmission_model # Dubbelkolla
filtered_distribution[,n]
next_predicted_dist




