library(bnlearn)
library(gRain)

data("asia")

# Structuring data
train_set = asia[1:4000,]
test_set = asia[4001:5000,]
test_labels = test_set["S"]
test_set = subset(test_set, select = -c(S) )
test_columns = colnames(test_set)

bn_structure_learning <- function() {
  
  model_structure = empty.graph(c("A", "S", "T", "L", "B", "E", "X", "D"))
  arc.set = matrix(c("A", "S", "T", "S", "L", "S", "B", "S", "E", "S", "X", "S", "D", "S"), ncol = 2, 
                   byrow = TRUE, dimnames = list(NULL, c("from", "to")))
  arcs(model_structure) = arc.set
  #plot(bn_nodes)
  
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
bn_naive_structure = bn_structure_learning()
plot(bn_naive_structure)
bn = bn_parameter_learning(bn_naive_structure, train_set)
preds = bn_precdict(bn, test_set, test_columns)

# Plotting result
confusion_matrix = table(preds, test_labels$S, dnn=c("Predictions", "True labels"))
plot(confusion_matrix)
confusion_matrix

# True BN network for asia dataset
true_asia_bn_struct = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
true_asia_bn = bn_parameter_learning(true_asia_bn_struct, train_set)
preds_asia = bn_precdict(true_asia_bn, test_set, test_columns)

# Plotting result, true BN
confusion_matrix_true_bn = table(preds_asia, test_labels$S, dnn=c("Predictions", "True labels"))
plot(confusion_matrix_true_bn)
confusion_matrix_true_bn
