library(tidyverse)
library(corrr)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(naniar)
library(viridis)
library(forcats)
library(ggridges)
library(knitr)
library("car")        
library("EnvStats")   
library(VIM)
library(mice)
library(Hmisc)
library(reshape2)
library(MASS)
library(caret)
library(AppliedPredictiveModeling)
library(pls)
library(lars)
library(elasticnet)
library(kernlab)
library(GGally)
library(mlbench)
library(randomForest)
library(ROSE)
library(rpart)
library(earth)
library(lubridate)
library(Metrics)
library(Amelia)
library(xgboost)
library("Rtsne")
library(plotly)
library(MASS)
library(klaR) #For Kmodes
library(cba) # For ROCK
library(cluster)
library(factoextra)
library(psych)
library(devtools)
library(ggord)

# Import required dataset
trainData <- read.csv(file = 'welldata.csv')
glimpse(trainData)

# Selecting numeric and Factor variable
dataNumeric <- trainData %>% dplyr::select_if(is.numeric)
dataFactor <- trainData %>% dplyr::select_if(negate(is.numeric))

# Create a function for data quality report for Numerical Variable
Q1 <- function(x, na.rm= TRUE) {
  quantile(x,na.rm=na.rm)[2]
}
Q3 <- function(x, na.rm= TRUE) {
  quantile(x,na.rm=na.rm)[4]
}

theNumericSummary <- function(x){
  c(length(x), n_distinct(x), sum(is.na(x)), mean(x,na.rm=TRUE),
    min(x,na.rm=TRUE), Q1(x,na.rm=TRUE), median(x,na.rm=TRUE),
    Q3(x,na.rm=TRUE), max(x,na.rm=TRUE), sd(x,na.rm=TRUE))
}

numericSummary <- dataNumeric %>% summarise(round(across(everything(),theNumericSummary)))
numericSummary <- cbind(stat = c("n", "unique", "missing", "mean", "min",
                                 "Q1", "median", "Q3", "max", "sd"), numericSummary)

numericSummaryFinal <- numericSummary %>%
  pivot_longer("DEPTH":"MW", names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(missing_pct = 100*missing/n, unique_pct = 100*unique/n) %>%
  dplyr::select(variable, n, missing, missing_pct, unique, unique_pct, everything())

#Final Data Quality Report Table
options(digits=3)
options(scipen=99)
numericSummaryFinal

#Factor Data Quality Report
getmodes <- function(v,type=1) {
  tbl <- table(v)
  m1<-which.max(tbl)
  if (type==1) {
    return (names(m1)) #1st mode
  }
  else if (type==2) {
    return (names(which.max(tbl[-m1]))) #2nd mode
  }
  else if (type==-1) {
    return (names(which.min(tbl))) #least common mode
  }
  else {
    stop("Invalid type selected")
  }
}

getmodesCnt <- function(v,type=1) {
  tbl <- table(v)
  m1<-which.max(tbl)
  if (type==1) {
    return (max(tbl)) #1st mode freq
  }
  else if (type==2) {
    return (max(tbl[-m1])) #2nd mode freq
  }
  else if (type==-1) {
    return (min(tbl)) #least common freq
  }
  else {
    stop("Invalid type selected")
  }
}

#create a new function called myfactorsummary
myFactorSummary <- function(x){
  c(length(x), n_distinct(x), sum(is.na(x)), round(getmodesCnt(x, type=1)/getmodesCnt(x,type=2),2),
    getmodes(x, type=1), getmodesCnt(x, type=1), getmodes(x, type=2), getmodesCnt(x,type=2),
    getmodes(x, type=-1), getmodesCnt(x, type=-1))
}

#apply function to all columns of housing factor and summarize 
factorSummary <- dataFactor %>% summarise(across(everything(), myFactorSummary))


# add labels
factorSummary <- cbind(stat = c("n", "unique", "missing", "freqRatio", "1st mode",
                                "1st mode freq", "2nd mode", "2nd mode freq", 
                                "least common", "least common freq"), factorSummary)

# create factor report 2
factorSummaryFinal <- factorSummary %>%
  pivot_longer("Lithology":"Well_ID", names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(missing_pct = 100*as.double(missing)/as.double(n),
         unique_pct = 100*as.double(unique)/as.double(n)) %>%
  dplyr::select(variable, n, missing, missing_pct, unique, unique_pct, everything())

#visualize report 2
options(digits=3)
options(scipen=99)
factorSummaryFinal

test <- subset(trainData, select = -c(DEPTH, Well_ID,Lithology))
ggcorr(test, label=T)

# or "inclusive", or "linear" by default
fig <- plot_ly(y = trainData$ROP, type = "box", quartilemethod="linear") 
fig <- fig %>% layout(title = 'Outlier Identification in ROP',
                        xaxis = list(title = ''),
                        yaxis = list(title = 'Rate of Penetration (ft/hr)'))
fig

# or "inclusive", or "linear" by default
fig <- plot_ly(y = trainData$WOB, type = "box", quartilemethod="linear") 
fig <- fig %>% layout(title = 'Outlier Identification in WOB',
                        xaxis = list(title = ''),
                        yaxis = list(title = 'Weight on Bit (klbs)'))
fig

# or "inclusive", or "linear" by default
fig <- plot_ly(y = trainData$Flow_Rate, type = "box", quartilemethod="linear") 
fig <- fig %>% layout(title = 'Outlier Identification in Flow Rate',
                        xaxis = list(title = ''),
                        yaxis = list(title = 'Flow Rate (gpm)'))
fig

# or "inclusive", or "linear" by default
fig <- plot_ly(y = trainData$Pump_Pressure, type = "box", quartilemethod="linear") 
fig <- fig %>% layout(title = 'Outlier Identification in Pump Pressure',
                        xaxis = list(title = ''),
                        yaxis = list(title = 'Pump Pressure (psi)'))
fig

# or "inclusive", or "linear" by default
fig <- plot_ly(y = trainData$RPM, type = "box", quartilemethod="linear") 
fig <- fig %>% layout(title = 'Outlier Identification in Rotary Speed',
                        xaxis = list(title = ''),
                        yaxis = list(title = 'Rotary Speed (RPM)'))
fig

# or "inclusive", or "linear" by default
fig <- plot_ly(y = trainData$Torque, type = "box", quartilemethod="linear") 
fig <- fig %>% layout(title = 'Outlier Identification in Torque',
                        xaxis = list(title = ''),
                        yaxis = list(title = 'Torque (psi)'))
fig

# or "inclusive", or "linear" by default
fig <- plot_ly(y = trainData$MW, type = "box", quartilemethod="linear") 
fig <- fig %>% layout(title = 'Outlier Identification in Mud Weight',
                        xaxis = list(title = ''),
                        yaxis = list(title = 'Mud Weight (ppg)'))
fig

#chart to show revenue by hour
trainData %>%
  mutate(name = reorder(Lithology,-DEPTH, .desc = TRUE)) %>%
  ggplot( aes(x = reorder(Lithology, -DEPTH), y = DEPTH)) +
  stat_summary(geom = "bar") +
  coord_flip() +
  xlab("Lithology") +
  theme_bw()

#chart to show revenue by hour
trainData %>%
  mutate(name = reorder(Lithology,-ROP, .desc = TRUE)) %>%
  ggplot( aes(x = reorder(Lithology, -ROP), y = ROP)) +
  stat_summary(geom = "bar", fun.y = "mean") +
  coord_flip() +
  xlab("Lithology") +
  theme_bw()

#chart to show revenue by hour
trainData %>%
  mutate(name = reorder(Lithology,-Torque, .desc = TRUE)) %>%
  ggplot( aes(x = reorder(Lithology, -Torque), y = Torque)) +
  stat_summary(geom = "bar", fun.y = "mean") +
  coord_flip() +
  xlab("Lithology") +
  theme_bw()

# Check the amount of missingness in training data
missmap(trainData, margins = c(15,5))
colSums(is.na(trainData))

# Mud weight missing value imputation (MEAN)
trainData$MW[is.na(trainData$MW)] <- mean(trainData$MW, na.rm=T)

#convert character data to factor
trainData <- trainData %>% mutate_if(is.character, as.factor) 

# Outlier Treatment (Remove outlier using domain knowledge and capping method)
trainData <- subset(trainData, ROP < 750)
trainData <- subset(trainData, ROP >= 0) 
trainData <- subset(trainData, WOB > 0) 
trainData <- subset(trainData, Pump_Pressure < 4648) 
trainData <- subset(trainData, Flow_Rate <= 1100) 
trainData <- subset(trainData, Flow_Rate > 200) 
trainData <- subset(trainData, RPM > 0) 
trainData <- subset(trainData, RPM <= 150) 
trainData <- subset(trainData, Torque > 0) 
glimpse(trainData)

hist.data.frame(trainData)

levels(trainData$Lithology)

# Standard scaler based on wellID
scaleData <- trainData %>% 
  group_by(Well_ID) %>%
  mutate(ROP = scale(ROP, center = TRUE, scale = TRUE)) %>%
  mutate(WOB = scale(WOB, center = TRUE, scale = TRUE)) %>%
  mutate(Pump_Pressure = scale(Pump_Pressure, center = TRUE, scale = TRUE)) %>%
  mutate(Flow_Rate = scale(Flow_Rate, center = TRUE, scale = TRUE)) %>%
  mutate(RPM = scale(RPM, center = TRUE, scale = TRUE)) %>%
  mutate(Torque = scale(Torque, center = TRUE, scale = TRUE)) %>%
  mutate(MW = scale(MW, center = TRUE, scale = TRUE)) %>%
  ungroup()

hist.data.frame(scaleData)

# Create data samples for clustering analysis
Lith_cluster = subset(trainData, select = c(Lithology))
#make this example reproducible
set.seed(1)
# Define one-hot encoding function
dummy <- dummyVars("~.", data = Lith_cluster)
# Implement the one-hot encoding function
Lith_cluster_ohe <- data.frame(predict(dummy, newdata=Lith_cluster))
head(Lith_cluster_ohe)

fviz_nbclust(Lith_cluster_ohe, kmeans, method = "silhouette") +
geom_vline(xintercept = 0, linetype = 2)

# Compute k-means with k = 10
set.seed(123)
Lith_cluster_kmeans <- kmeans(Lith_cluster_ohe, 7, nstart = 25)

result_Lith_cluster_kmeans <- table(Lith_cluster$Lithology, Lith_cluster_kmeans$cluster)
Lith_purity.kmean <- sum(apply(result_Lith_cluster_kmeans, 2, max)) / nrow(Lith_cluster)
Lith_purity.kmean

result_Lith_cluster_kmeans

fviz_nbclust(Lith_cluster_ohe, kmodes, method = "silhouette") +
geom_vline(xintercept = 0, linetype = 2)

# Compute k-means with k = 9
set.seed(123)
Lith_cluster_kmodes <- kmodes(Lith_cluster_ohe, 9)

result_Lith_cluster_kmodes <- table(Lith_cluster$Lithology, Lith_cluster_kmodes$cluster)
Lith_purity.kmodes <- sum(apply(result_Lith_cluster_kmodes, 2, max)) / nrow(Lith_cluster)
Lith_purity.kmodes

result_Lith_cluster_kmodes

# Create data samples for clustering analysis
ROP_Lith_cluster = subset(scaleData, select = c(ROP, Lithology))
#make this example reproducible
set.seed(1)
# Define one-hot encoding function
dummy <- dummyVars("~.", data = ROP_Lith_cluster)
# Implement the one-hot encoding function
ROP_Lith_cluster_ohe <- data.frame(predict(dummy, newdata=ROP_Lith_cluster))
head(ROP_Lith_cluster_ohe)

fviz_nbclust(ROP_Lith_cluster_ohe, kmeans, method = "silhouette") +
geom_vline(xintercept = 0, linetype = 2)

# Compute k-means with k = 10
set.seed(123)
ROP_Lith_cluster_kmeans <- kmeans(ROP_Lith_cluster_ohe, 10, nstart = 25)

result_ROP_Lith_cluster_kmeans <- table(ROP_Lith_cluster$Lithology, ROP_Lith_cluster_kmeans$cluster)
ROP_Lith_purity.kmean <- sum(apply(result_ROP_Lith_cluster_kmeans, 2, max)) / nrow(ROP_Lith_cluster)
ROP_Lith_purity.kmean

result_ROP_Lith_cluster_kmeans

dataNumFact = subset(scaleData, select = -c(DEPTH, Well_ID))
#make this example reproducible
set.seed(1)
# Define one-hot encoding function
dummy <- dummyVars("~.", data = dataNumFact)
# Implement the one-hot encoding function
dataNumFact_ohe <- data.frame(predict(dummy, newdata=dataNumFact))
head(dataNumFact_ohe)

fviz_nbclust(dataNumFact_ohe, kmeans, method = "silhouette") +
geom_vline(xintercept = 0, linetype = 2)

# Compute k-means with k = 2
set.seed(123)
Num_Fact_cluster_kmeans <- kmeans(dataNumFact_ohe, 2, nstart = 25)

result_Num_Fact_cluster_kmeans <- table(dataNumFact$Lithology, Num_Fact_cluster_kmeans$cluster)
Num_Fact_purity.kmean <- sum(apply(result_Num_Fact_cluster_kmeans, 2, max)) / nrow(dataNumFact)
Num_Fact_purity.kmean

result_Num_Fact_cluster_kmeans

# Create numeric dataframe for PCA analysis
data_lda = subset(scaleData, select = -c(DEPTH, Well_ID))

# Create pair correlation plot 
pairs.panels(data_lda[1:7],
             gap = 0,
             bg = data_lda$Lithology,
             pch = 21)

#make this example reproducible
set.seed(1)

#use 80% of dataset as training set and 20% as test set
sample_lda <- sample(c(TRUE, FALSE), nrow(data_lda), replace=TRUE, prob=c(0.80,0.20))
train_lda  <- data_lda[sample_lda, ]
test_lda   <- data_lda[!sample_lda, ]

linear <- lda(Lithology~., train_lda)
linear

attributes(linear)

p <- predict(linear, train_lda)
par(mar = c(1, 5, 1, 5))
ldahist(data = p$x[,1], g = train_lda$Lithology)

ggord(linear, train_lda$Lithology, ylim = c(-10, 10))

# Partition plot This provides the classification of each and every combination in the training dataset.
partimat(Lithology~., data = train_lda, method = "lda")

p1 <- predict(linear, train_lda)$class
tab <- table(Predicted = p1, Actual = train_lda$Lithology)
tab

sum(diag(tab))/sum(tab)

p2 <- predict(linear, test_lda)$class
tab1 <- table(Predicted = p2, Actual = test_lda$Lithology)
tab1

sum(diag(tab1))/sum(tab1)

moh_scale <- trainData$Lithology
scaleData <- cbind(trainData, moh_scale)
glimpse(trainData)

levels(trainData$Lithology)

trainData$moh_scale <- as.character(trainData$moh_scale)
trainData$moh_scale[which(trainData$moh_scale == 'Alluvium' )] <- '3'
trainData$moh_scale[which(trainData$moh_scale) == 'Andesite' ] <- '6.17'
trainData$moh_scale[which(trainData$moh_scale) == 'Clay' ] <- '2.25'
trainData$moh_scale[which(trainData$moh_scale) == 'Dacite' ] <- '6.35'
trainData$moh_scale[which(trainData$moh_scale) == 'Diorite' ] <- '7'
trainData$moh_scale[which(trainData$moh_scale) == 'Felsic Dike' ] <- '6.5'
trainData$moh_scale[which(trainData$moh_scale) == 'Granite' ] <- '6.54'
trainData$moh_scale[which(trainData$moh_scale) == 'Granodiorite' ] <- '6.4'
trainData$moh_scale[which(trainData$moh_scale) == 'Monzonite' ] <- '6.5'
trainData$moh_scale[which(trainData$moh_scale) == 'Plutonic' ] <- '6.5'
trainData$moh_scale[which(trainData$moh_scale) == 'Quartz Diorite' ] <- '7'
trainData$moh_scale[which(trainData$moh_scale) == 'Quartz Monzonite' ] <- '7'
trainData$moh_scale[which(trainData$moh_scale) == 'Rhyolite' ] <- '6.55'
trainData$moh_scale[which(trainData$moh_scale) == 'Sandstone clay' ] <- '5'
trainData$moh_scale <- as.numeric(trainData$moh_scale)

summary(scaleData)

# Create target label
y_pca = subset(scaleData, select = c(Lithology))

# Create numeric dataframe for PCA analysis
x_pca = subset(scaleData, select = -c(DEPTH, Lithology, Well_ID))

# Analyze the proportion of variance from the different principal components
pca <- prcomp(x_pca, center = FALSE,scale. = FALSE)
summary(pca)

# Convert numerical data to principal components 
X_pca <- predict(pca,x_pca)

# Select the 7 principal components since it explains 94.6% data
X_pca <- data.frame(X_pca[, 1:7])

# Merge PCA and target data
data_pca <- cbind(X_pca,y_pca)


