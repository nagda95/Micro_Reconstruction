################################################################################
#
# Copulas : Battery Demo
# 
# Presented 2021-09-10
#
# created by : August Brandberg augustbr@kth.se
# date : 2021-09-10
#
#
#
#
# Installation steps:
#
# 	1. Install R (https://cran.r-project.org/bin/)
#   2. Recommended but not required: Install RStudio (https://www.rstudio.com/products/rstudio/download/#download)
#   3. Run this script using R.
#
# created by : August Brandberg augustbr at kth dot se
# date : 2021-08-02
#
#

# Install packages from CRAN
list.of.packages <- c("VineCopula", "CDVineCopulaConditional", "fitdistrplus","devtools","")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Install my package directly from GITHUB.
library("devtools")
install_github("abrandberg/pulpDistr")

# Clear workspace environment
rm(list=ls())

# Import packages
library("VineCopula")               # Statistical inference of vine copulas
library("CDVineCopulaConditional")  # Automatic conditional sampling
library("fitdistrplus")             # Statistical inference of parametric distributions
library("pulpDistr")

# Import data
#data = as.data.frame(read.csv("rawData.csv",header = FALSE))
#data = as.data.frame(read.csv("rawData_pred_r.csv",header = FALSE))

data = as.data.frame(read.csv("raw_data_jitter_EXP_Copula_np.csv",header = FALSE))
predict_data = as.data.frame(read.csv("predict_data01_jitter_contacts_np.csv",header = FALSE))

data = as.data.frame(read.csv("stardist_RSA_Contacts_rawData.csv",header = FALSE))
predict_data = as.data.frame(read.csv("stardist_RSAJulia_Contacts_pred.csv",header = FALSE))

# Remove the nan values
data <-unstack(within(stack(data), values[is.nan(values)] <- 0))

# Reorder data so that the data which will be supplied during the sampling stage
# is at the end.
data <- data[ c(2,1,3)]

# Assign names to the data frame
names(data) <- c("numContacts","Radius","numNeighbours")
#cor(data,method = "kendall")

# Split into test and train 50:50
dataTrain = data[1:6725 ,]
dataTest = data[4001:6725,]
totalDataTest = data[4001:6725,]

# Generate pseudo-observations
udata = generatePseudoU(dataTrain)


# Inputs. See documentation.
condVars = 2
copType = "CVine"
selMethod = "AIC"
checkIndep = TRUE
#fitSet = c(0,1)
fitSet = NA


# Fit copula
RVM <- CDVineCondFit(udata,Nx=condVars, type=copType, selectioncrit=selMethod,
                     indeptest=checkIndep, level=0.05, familyset = fitSet, rotations = TRUE)

# Optionally generate data from the fit, without supplying conditional data
SimX <- RVineSim(1*dim(udata)[1],RVM)
#cor(SimX, method = "kendall")

print("******************** Margin fits ********************")
listOfModels = c("weibull","gamma","lnorm","norm","exp")
fitRadius      = fitMarginFunction(    dataTrain$Radius, listOfModels, SimX[,2])
fitNumNeighbours   = fitMarginFunction( dataTrain$numNeighbours, listOfModels, SimX[,3])
fitNumContact    = fitMarginFunction(  dataTrain$numContacts, listOfModels, SimX[,1])


##########################################
# Custom routine to fit number of contacts
# Matrix with two columns:
# 1. Sorted value
# 2. Percentile
matrixForNumContact = as.matrix(cbind(sort(dataTrain$numContacts), 1:length(dataTrain$numContacts)/(length(dataTrain$numContacts))))
matrixForRadius = as.matrix(cbind(sort(dataTrain$Radius), 1:length(dataTrain$Radius)/(length(dataTrain$Radius))))
matrixForNumNeighbours = as.matrix(cbind(sort(dataTrain$numNeighbours), 1:length(dataTrain$numNeighbours)/(length(dataTrain$numNeighbours))))


##########################################



# Split the data that will be supplied
dataX = dataTrain[,c(2,3)]

# Answers
dataAnswer = dataTest[,1]
dataTest = dataTest[,c(2,3)]
dataTest=predict_data
names(dataTest) <- c("Radius","numNeighbours")

################################################################################
# Conditional sampling

idxToRetain = length(dataTest[,1])
denomNormalization = length(dataX[,1])+1+1
CC = matrix(nrow = idxToRetain,ncol = 2)
for (j in 1:dim(dataTest)[2]){
  ap = dataX[,j]
  ap3 = sort(ap)
  for(i in 1:idxToRetain ){
    CC[i,j] = min(min(which(ap3 > dataTest[i,j] )),length(dataX[,1])+1)/denomNormalization
  }
}



xTemp = generatePseudoU(dataTest)
conditioningData = CC

d=dim(RVM$Matrix)[1]
condition = conditioningData[, cbind(RVM$Matrix[(d+1)-1,(d+1)-1] , RVM$Matrix[(d+1)-2,(d+1)-2] )-1]
Sim <- CDVineCondSim(RVM,condition)

#listOfModels = c("weibull","gamma","lnorm","norm","exp")
#fitRadius      = fitMarginFunction(    dataTest$Radius, listOfModels, xTemp[,1])
#fitNumNeighbours   = fitMarginFunction( dataTest$numNeighbours, listOfModels, xTemp[,2])

# Margin draw custom for fitRadius
tempMatrix = matrix(data = NA, nrow = length(Sim[,1]), ncol = 1)
for (aLoop in 1:length(Sim[,1])){
  tempMatrix[aLoop] <- matrixForRadius[which.min(Sim[aLoop,2] > matrixForRadius[,2]) , 1]
}
fitRadius$marginDraw = tempMatrix

# Margin draw custom for fitNumNeighbours
tempMatrix = matrix(data = NA, nrow = length(Sim[,1]), ncol = 1)
for (aLoop in 1:length(Sim[,1])){
  tempMatrix[aLoop] <- matrixForNumNeighbours[which.min(Sim[aLoop,3] > matrixForNumNeighbours[,2]) , 1]
}
fitNumNeighbours$marginDraw = tempMatrix

# Margin draw custom for fitNumContact
tempMatrix = matrix(data = NA, nrow = length(Sim[,1]), ncol = 1)
for (aLoop in 1:length(Sim[,1])){
    tempMatrix[aLoop] <- matrixForNumContact[which.min(Sim[aLoop,1] > matrixForNumContact[,2]) , 1]
}
fitNumContact$marginDraw = tempMatrix

bp = cbind(fitNumContact$marginDraw, fitRadius$marginDraw,fitNumNeighbours$marginDraw)

cor(dataTrain,method = "kendall")
cor(totalDataTest,method = "kendall")
cor(Sim,method = "kendall")
cor(bp,method = "kendall")




# Export step
bp = bp[,c(2,1,3)]
dataTrain = dataTrain[,c(2,1,3)]
totalDataTest = totalDataTest[,c(2,1,3)]

# Export the data for plotting in MATLAB
write.csv(bp, file = "bp_stardist_RSAJulia_Contacts.csv", quote = FALSE, row.names = FALSE)
write.csv(dataTrain, file = "trainData.csv", quote = FALSE, row.names = FALSE)
write.csv(totalDataTest, file = "testData.csv", quote = FALSE, row.names = FALSE)
