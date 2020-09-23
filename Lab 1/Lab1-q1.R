library(bnlearn)

data("asia")
set.seed(1)

# Initial run
model = hc(asia, restart = 1)
plot(model, main= 0)

# Runs the hc-algorithm 'i' times and if the model changes, 
# the network gets plotted.
for (i in 1:20) {
  updated_model = hc(asia, start = model, restart = 1)

  
  if (all.equal(updated_model, model) != TRUE) {
    plot(updated_model, main = i)
  }
  
  model = updated_model
}
