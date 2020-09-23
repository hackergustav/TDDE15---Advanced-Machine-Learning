library(bnlearn)
library(gRain)

data("asia")

# Structuring data
train_set = asia[1:4000,]
test_set = asia[4001:5000,]
test_labels = test_set["S"]
test_set = subset(test_set, select = -c(S) )
test_columns = colnames(test_set)

bn_structure_learning <- function(dataset, restarts) {
  
  # Structure learning
  model_structure = hc(dataset, restart = restarts)
  #plot(model_structure)
  
  return(model_structure)
}

bn_parameter_learning <- function(bn_structure, training_set) {
  
  # Parameter learning
  model = bn.fit(bn_structure, training_set)
  
  # Converting model to grain objects (avoiding NaN-values?)
  grained_model = as.grain(model)
  
  # Compiling model to an BN?
  compiled_model = compile(grained_model)
  
  return(compiled_model)
}

bn_precdict <- function(bn_model, testing_set, columns) {
  
  predictions = c()
  
  print(columns)
  
  # Predictions
  for (i in 1:dim(testing_set)[1]) {
    
    # Collecting ith obs states in a vector 
    states = NULL
    for (j in columns) {
      if (testing_set[i,j] == "yes") {
        states = c(states, "yes")
      }  else {
        states = c(states, "no")
      }
    }
    
    
    # Prediction of var "S" on the ith obs
    evidence = setEvidence(bn_model, nodes = columns, states = states)
    
    probabilities = querygrain(evidence, c("S"))
    
    if(probabilities$S["no"] >= 0.5) {
      predictions[i] = "no"
    } else {
      predictions[i] = "yes"
    }
  }
  return(predictions)
}

# Learning a BN from dataset
bn_structure = bn_structure_learning(asia, 5)
bn = bn_parameter_learning(bn_structure, train_set)
preds = bn_precdict(bn, test_set, mb(bn_structure, c("S")))

# Plotting result
confusion_matrix_mb = table(preds, test_labels$S, dnn=c("Predictions", "True labels"))
plot(confusion_matrix_mb)
confusion_matrix_mb


