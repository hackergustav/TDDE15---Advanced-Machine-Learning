# TASK 2.3.1

library(kernlab)


data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/
GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")

names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")

data[,5] <- as.factor(data[,5])

set.seed(111); 

SelectTraining <- sample(1:dim(data)[1], size = 1000, replace = FALSE)
SelectTest <- sample(1:dim(data)[1], size = 1000, replace = FALSE)

trainingData = data[SelectTraining,]
testingData = data[SelectTest,]

GPfit = gausspr(fraud ~  varWave + skewWave, data=trainingData)

x1 = seq(min(trainingData$varWave), max(trainingData$varWave), length=100)
x2 = seq(min(trainingData$skewWave), max(trainingData$skewWave), length=100)
gridXY = meshgrid(x1, x2)
gridXY <- data.frame(cbind(c(gridXY$x), c(gridXY$y)))
names(gridXY) <- c("varWave", "skewWave")
prediction = predict(GPfit, gridXY, type="probabilities")

# Plotting for Prob
contour(x1,x2,matrix(prediction[,1],100,byrow = TRUE), 20, xlab = "varWave", ylab = "skewWave")
points(trainingData[trainingData[, 5]==0, 1], trainingData[trainingData[, 5]==0, 2], col="red")
points(trainingData[trainingData[, 5]==1, 1], trainingData[trainingData[, 5]==1, 2], col="blue")

confusionMatrix = table(predict(GPfit,trainingData[,1:4]), trainingData[,5]) # confusion matrix
confusionMatrix
accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix)
accuracy

# TASK 2.3.2
prediction = predict(GPfit, testingData)
confusionMatrix = table(prediction, testingData[,5]) # confusion matrix
confusionMatrix
accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix)
accuracy

# TASK 2.3.3
GPfit = gausspr(fraud ~  varWave + skewWave + kurtWave + entropyWave, data=trainingData)

prediction = predict(GPfit, testingData)
confusionMatrix = table(prediction, testingData[,5]) # confusion matrix
confusionMatrix
accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix)
accuracy

