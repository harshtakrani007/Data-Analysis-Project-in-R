#project data input

library(data.table)
library(bit)
install.packages("bit64")
library(bit64)
library(graphics)
#computing NPS
install.packages("gsubfn")
library(gsubfn)
install.packages("proto")
library(proto)
library(RSQLite)
install.packages("sqldf")
library(sqldf)
install.packages("arules")
library(arules)
install.packages("grid")
library(grid)
install.packages("arulesViz")
library(arulesViz)

#detach("package:arules", unload=TRUE)


FebData <- fread(file= "FebData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
MarchData <- fread(file= "MarchData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
AprilData <- fread(file = "AprilData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE   )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
MayData <- fread(file = "MayData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE   )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
JuneData <- fread(file = "JuneData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
JulyData <- fread(file = "JulyData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
AugustData <- fread(file = "AugustData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
SeptemberData <- fread(file = "SeptData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
OctoberData <- fread(file = "OctData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
NovemberData <- fread(file = "NovData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
DecemberData <- fread(file = "DecData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]
JanuaryData <- fread(file = "JanData.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE    )[,c(11,12,23,26,28,30,32,65,89,137:145,147,168,169,171,179,202,203,204,205,210,213,214,215,216,218,223,232)]

str(FebData)

#merging whole data
FullData <- rbind(FebData,MarchData,AprilData,MayData,JuneData,JulyData,AugustData,SeptemberData,OctoberData,NovemberData,DecemberData,JanuaryData)
str(FullData)
#removing row names
rownames(FullData) <- NULL
dim(FullData)
head(FullData)

#removing NA's based on NPS and likelihood to recommend
FullData <- FullData[!(FullData$NPS_Type == "" | is.na(FullData$Likelihood_Recommend_H)), ]
#checking
head(FullData$Likelihood_Recommend_H,20)
head(FullData$NPS_Type,20)
head(FullData)
#rm(FullData)
#removing columns with more than 60% NA values
#FullData <-  FullData[, -which(colMeans(is.na(FullData)) > 0.4)]


View(colMeans(is.na(FullData)) > 0.6)

#adding a likelihood_recommend_h_char column in full data
FullData$Likelihood_Recommend_H_Char <- FullData$Likelihood_Recommend_H
FullData$Likelihood_Recommend_H_Char <- as.factor(FullData$Likelihood_Recommend_H_Char)
FullData$Likelihood_Recommend_H_Char <- gsub(pattern = "9|10", replacement = "high" ,FullData$Likelihood_Recommend_H_Char)
FullData$Likelihood_Recommend_H_Char <- gsub(pattern = "7|8", replacement = "medium" ,FullData$Likelihood_Recommend_H_Char)
FullData$Likelihood_Recommend_H_Char <- gsub(pattern = "1|2|3|4|5|6", replacement = "low" ,FullData$Likelihood_Recommend_H_Char)

str(FullData)
#########################################################################################################

#adding a numeric type NPS_Type
FullData$NPS_Type_num <- FullData$NPS_Type
FullData$NPS_Type_num <- as.factor(FullData$NPS_Type_num)
FullData$NPS_Type_num <- gsub(pattern = "Promoter",replacement = "10",FullData$NPS_Type_num)
FullData$NPS_Type_num <- gsub(pattern = "Detractor",replacement = "6",FullData$NPS_Type_num)
FullData$NPS_Type_num <- gsub(pattern = "Passive",replacement = "8",FullData$NPS_Type_num)




############################################################################################################

#NPS for whole dataset

promoters <- length(which(FullData$NPS_Type == "Promoter"))
View(promoters)
detractors <- length(which(FullData$NPS_Type == "Detractor"))
View(detractors)
passives <- length(which(FullData$NPS_Type == "Passive"))
View(passives)

total <- passives + promoters + detractors
View(total)
NPS <- ((promoters-detractors)*100)/total
View(NPS)
#output
#57.17
############################################################################
memory.limit(100000)

# %of detractors
percentage_Detractors <- function(countrydata){
  promoter <- length(which(countrydata$NPS_Type == "Promoter"))
  detractor <- length(which(countrydata$NPS_Type == "Detractor"))
  passives <- length(which(countrydata$NPS_Type == "Passive"))
  total <- promoter + detractor + passives
  percent_detractor <- (detractor/total)*100
  return(paste( "percentage of detractors are",percent_detractor))
}

# NPS function
percentage_DetractorsAndNPS <- function(countrydata){
  promoter <- length(which(countrydata$NPS_Type == "Promoter"))
  detractor <- length(which(countrydata$NPS_Type == "Detractor"))
  passives <- length(which(countrydata$NPS_Type == "Passive"))
  total <- promoter + detractor + passives
  percent_detractor <- (detractor/total)*100
  NPS <- ((promoter - detractor)/total)* 100
  return(paste( "percentage of detractors are", percent_detractor, "and NPS score is", NPS))
}
View(FullData$Guest_Room_H)
##################################################################################################
#adding character columns for each metric
FullData$Guest_Room_H_char <- FullData$Guest_Room_H
FullData$Guest_Room_H_char <- as.factor(FullData$Guest_Room_H_char)
FullData$Guest_Room_H_char <- gsub("10|9","high",FullData$Guest_Room_H_char)
FullData$Guest_Room_H_char <- gsub("8|7","medium",FullData$Guest_Room_H_char)
FullData$Guest_Room_H_char <- gsub("0|1|2|3|4|5|6","low",FullData$Guest_Room_H_char)



FullData$Condition_Hotel_H_char <- FullData$Condition_Hotel_H
FullData$Condition_Hotel_H_char <- as.factor(FullData$Condition_Hotel_H_char)
FullData$Condition_Hotel_H_char <- gsub("10|9","high",FullData$Condition_Hotel_H_char)
FullData$Condition_Hotel_H_char <- gsub("8|7","medium",FullData$Condition_Hotel_H_char)
FullData$Condition_Hotel_H_char <- gsub("0|1|2|3|4|5|6","low",FullData$Condition_Hotel_H_char)

FullData$Customer_SVC_H_char <- FullData$Customer_SVC_H
FullData$Customer_SVC_H_char <- as.factor(FullData$Customer_SVC_H_char)
FullData$Customer_SVC_H_char <- gsub("10|9","high",FullData$Customer_SVC_H_char)
FullData$Customer_SVC_H_char <- gsub("8|7","medium",FullData$Customer_SVC_H_char)
FullData$Customer_SVC_H_char <- gsub("0|1|2|3|4|5|6","low",FullData$Customer_SVC_H_char)

FullData$Staff_Cared_H_char <- FullData$Staff_Cared_H
FullData$Staff_Cared_H_char <- as.factor(FullData$Staff_Cared_H_char)
FullData$Staff_Cared_H_char <- gsub("10|9","high",FullData$Staff_Cared_H_char)
FullData$Staff_Cared_H_char <- gsub("8|7","medium",FullData$Staff_Cared_H_char)
FullData$Staff_Cared_H_char <- gsub("0|1|2|3|4|5|6","low",FullData$Staff_Cared_H_char)

FullData$Internet_Sat_H_char <- FullData$Internet_Sat_H
FullData$Internet_Sat_H_char <- as.factor(FullData$Internet_Sat_H_char)
FullData$Internet_Sat_H_char <- gsub("10|9","high",FullData$Internet_Sat_H_char)
FullData$Internet_Sat_H_char <- gsub("8|7","medium",FullData$Internet_Sat_H_char)
FullData$Internet_Sat_H_char <- gsub("0|1|2|3|4|5|6","low",FullData$Internet_Sat_H_char)

FullData$Check_In_H_char <- FullData$Check_In_H
FullData$Check_In_H_char <- as.factor(FullData$Check_In_H_char)
FullData$Check_In_H_char <- gsub("10|9","high",FullData$Check_In_H_char)
FullData$Check_In_H_char <- gsub("8|7","medium",FullData$Check_In_H_char)
FullData$Check_In_H_char <- gsub("0|1|2|3|4|5|6","low",FullData$Check_In_H_char)



View(colMeans(is.na(FullData)) > 0.6)
FullData$Overall_Sat_H <- NULL
FullData$Tranquility_H <- NULL
str(FullData)

colnames(FullData)

#################################################################################################################################
FranceData <- subset(FullData, FullData$Country_PL=="France")
str(FranceData)

#FranceData
################################################################################################################################
#linear modelling


FranceData <- na.omit(FranceData)
View(FranceData)
model1 <- lm(data = FranceData[1:100000,], formula = Likelihood_Recommend_H ~ Guest_Room_H + Condition_Hotel_H + Customer_SVC_H + Staff_Cared_H + Check_In_H)
summary(model1)
library(ggplot2)

model2 <- lm(data = FranceData, formula = PMS_TOTAL_REV_USD_C ~ PMS_ROOM_REV_USD_C)
summary(model2)
#Multiple R-squared:  0.9461,	Adjusted R-squared:  0.9461 
#here we can see that total revenue is based on heavily relied on room revenue in hyatt hotels in france

model3 <- lm(data = FranceData, formula = PMS_TOTAL_REV_USD_C ~ PMS_FOOD_BEVERAGE_REV_USD_C)
summary(model3)
#Multiple R-squared:  0.4276,	Adjusted R-squared:  0.4275 
#can improve food and beverage services to improve revenue

model4 <- lm(data = FranceData, formula = PMS_TOTAL_REV_USD_C ~ PMS_OTHER_REV_USD_C + REVENUE_USD_R )
summary(model4)
#Multiple R-squared:  0.7828,	Adjusted R-squared:  0.7828 


ggplot(FranceData, aes(x=FranceData$Guest_Room_H, y=FranceData$Likelihood_Recommend_H  , color=FranceData$NPS_Type)) + 
  geom_smooth(method = "lm") +  ylab("LTR Rating") +  xlab(" Guest Room satisfaction metric") +  
  ggtitle(" Effect of Guest Room Rating rating on Net Promoter Score")

ggplot(FranceData, aes(x=FranceData$Condition_Hotel_H, y=FranceData$Likelihood_Recommend_H  , color=FranceData$NPS_Type)) + 
  geom_smooth(method = "lm") +  ylab("LTR Rating") +  xlab(" Condition of hotel metric") +  
  ggtitle(" Effect of condition of hotel rating on Net Promoter Score")

ggplot(FranceData, aes(x=FranceData$Customer_SVC_H, y=FranceData$Likelihood_Recommend_H  , color=FranceData$NPS_Type)) + 
  geom_smooth(method = "lm") +  ylab("LTR Rating") +  xlab(" customer service satisfaction metric") +  
  ggtitle(" Effect of customer service rating on Net Promoter Score")

ggplot(FranceData, aes(x=FranceData$Check_In_H, y=FranceData$Likelihood_Recommend_H  , color=FranceData$NPS_Type)) + 
  geom_smooth(method = "lm") +  ylab("LTR Rating") +  xlab(" check in satisfaction metric") +  
  ggtitle(" Effect of check in rating on Net Promoter Score")




####################################################################################################################
#KVSM anD SVM
rm(data1)

colnames(FranceData)
data1 <- FranceData[1:100000,c(10:17)]
View(data1)
View(data1)


randIndx <- sample(1:dim(data1)[1])

#taking 2/3rd of random data for training dataset and remaining 1/3rd for testing dataset

point <- floor(2*dim(data1)[1]/3)

traindta <- data1[randIndx[1:point],]
View(traindta)
testdta <- data1[randIndx[(point+1):dim(data1)[1]],]
View(testdta)
str(traindta)



traindta <- na.omit(traindta)
testdta <- na.omit(testdta)
library(kernlab)
#step3
KSVM <- ksvm(Likelihood_Recommend_H ~.,data = traindta, kernel="rbfdot",kpar="automatic", C=5, cross=3, prob.model=TRUE)
KSVM 

#Objective Function Value : -1624.44 
#Training error : 0.120246 
#Cross validation error : 1.66121 
#Laplace distr. width : 2.633128


#here we can see that the training and cross validation error is very low
library(e1071)


SVM <- svm(data= traindta,Likelihood_Recommend_H~.)
SVM  



predict <- predict(SVM, testdta,type="votes")
View(predict)

compTable <- data.frame(predict,testdta[,1])
View(compTable)

#renaming columns
colnames(compTable) <- c("Pred","test")
#calculating RMSerror
error <- sqrt(mean(compTable$test - compTable$Pred)^2)
error
#0.07239554

#calculating error
compTable$error <- abs(compTable$test - compTable$Pred)

#creatting data.frame to plot
data <- data.frame(testdta$Likelihood_Recommend_H,testdta$Guest_Room_H,compTable$error)

library(ggplot2)
#plotting
ScatPlot <- ggplot(data,aes(x=testdta.Likelihood_Recommend_H,y=testdta.Guest_Room_H)) +
  geom_point(aes(size=compTable.error,color=compTable.error))
ScatPlot

##################################################################################################################################
#arules
colnames(FranceData)
str(FranceData)

rm(FranceData_Arules)

FranceData_Arules <- FranceData[,c (1,3,8,35)]
#FranceData_Arules <- gsub("","NA",FranceData_Arules)
FranceData_Arules <- na.omit(FranceData_Arules)
str(FranceData_Arules)
View(FranceData_Arules)

FranceData_Arules <-  replace(FranceData_Arules,TRUE,lapply(FranceData_Arules, factor))
FranceData_Arules <- sapply(FranceData_Arules,as.factor)

str(FranceData_Arules)
length(FranceData$POV_CODE_C)

rules <- apriori(FranceData_Arules, parameter = list(support = 0.09, confidence = 0.5))

inspect(rules)
summary(rules)
plot(rules,method= "graph",measure="support",shading="lift",engine= "interactive")

goodrules <- rules[quality(rules)$lift>1]
inspect(goodrules)
summary(goodrules)
plot(goodrules,method= "graph",measure="support",shading="lift",engine= "interactive")

#3
x <- is.maximal(goodrules)
inspect(goodrules[x])
#lhs                                     rhs                                  support confidence     lift count
#[1] {}                                   => {Likelihood_Recommend_H_Char=high} 0.5521922  0.5521922 1.000000  5655
#[2] {Likelihood_Recommend_H_Char=low}    => {POV_CODE_C=BUSINESS}              0.1915829  0.9091752 1.023284  1962
#[3] {Likelihood_Recommend_H_Char=medium} => {POV_CODE_C=BUSINESS}              0.2126745  0.8970346 1.009620  2178
#[4] {GUEST_COUNTRY_R=UNITED STATES}      => {Likelihood_Recommend_H_Char=high} 0.1794747  0.6313981 1.143439  1838
#[5] {POV_CODE_C=BUSINESS,                                                                                         
#GUEST_COUNTRY_R=UNITED STATES}      => {Likelihood_Recommend_H_Char=high} 0.1559418  0.6218847 1.126211  1597

#association rules 2
colnames(FranceData)
FranceData_Arules2 <- subset(FranceData[,c(3,24:35,36)])
FranceData_Arules2 <- na.omit(FranceData_Arules2)
str(FranceData_Arules2)
View(FranceData_Arules2)

FranceData_Arules2 <-  replace(FranceData_Arules2,TRUE,lapply(FranceData_Arules2, factor))
rules2 <- apriori(FranceData_Arules2,parameter = list(support=0.78, confidence= 0.6))
summary(rules2)
plot(rules2,method= "graph",measure="support",shading="lift",engine= "interactive")
inspect(rules2)

goodrules2 <- rules2[quality(rules2)$lift>1]
inspect(goodrules2)
summary(goodrules2)
plot(goodrules2,method= "graph",measure="support",shading="lift",engine= "interactive")

#3
x <- is.maximal(goodrules2)
inspect(goodrules2[x])

##################################################################################################################
###################################################################################################################
colnames(FranceData)
#random forest
install.packages("randomForest")
library(randomForest)
RandomforestData <- subset(FranceData[,c(10:17)])
str(RandomforestData)
View(RandomforestData)
RandomforestData <- na.omit(RandomforestData)

RandomforestData$Likelihood_Recommend_H <-  ifelse(RandomforestData$Likelihood_Recommend_H >= 6,1,0)
RandomforestData$`F&B_Overall_Experience_H` <- as.integer(RandomforestData$`F&B_Overall_Experience_H`)

model4 <- randomForest(Likelihood_Recommend_H ~ Guest_Room_H + Condition_Hotel_H + Staff_Cared_H + Internet_Sat_H  +  Customer_SVC_H + Staff_Cared_H , ntree=4 ,mtry=2,data = RandomforestData)
model4
#% Var explained: 71.02
plot(model4)


####################################################################################################################################
#USdata

USdata <- subset(FullData, FullData$Country_PL=="United States")
str(USdata)


###########################################################################################
#linear modelling
model5 <- lm(data= USdata, Likelihood_Recommend_H ~ Guest_Room_H + Condition_Hotel_H + Customer_SVC_H)
summary(model5)
#Multiple R-squared:  0.6829,	Adjusted R-squared:  0.6829


model8 <- lm(data = USdata, Likelihood_Recommend_H~ NPS_Type_num)
summary(model8)
#Multiple R-squared:  0.8656,	Adjusted R-squared:  0.8656 




###################################################################################################################
#A rules
colnames(USdata)


USdata_arules <- subset(USdata[,c(3,22:34)])
USdata_arules <- na.omit(USdata_arules)
str(USdata_arules)
View(USdata_arules)

USdata_arules <-  replace(USdata_arules,TRUE,lapply(USdata_arules, factor))
rules3 <- apriori(USdata_arules,parameter = list(support=0.9, confidence= 0.8))
summary(rules3)
plot(rules3)
plot(rules3,method= "graph",measure="support",shading="lift",engine= "interactive")
inspect(rules3)

goodrules3 <- rules3[quality(rules3)$lift>1]
inspect(goodrules3)
summary(goodrules3)
plot(goodrules3)

#3
x3 <- is.maximal(goodrules3)
inspect(goodrules3[x3])
##########################################################################################

rm(USleisureData)

USleisureData <- subset(USdata, USdata$POV_CODE_C=="LEISURE")
View(USleisureData)
unique(USleisureData$POV_CODE_C)
colnames(USleisureData)

USleisureData_arules <- subset(USleisureData[,c(24:34,36)])
USleisureData_arules <-  replace(USleisureData_arules,TRUE,lapply(USleisureData_arules, factor))
str(USleisureData_arules)

rules4 <- apriori(USleisureData_arules,parameter = list(support=0.6, confidence= 0.6))
summary(rules4)
plot(rules4)
inspect(rules4)

goodrules4 <- rules4[quality(rules4)$lift>1]
inspect(goodrules4)
summary(goodrules4)
plot(goodrules4)

#3
x4 <- is.maximal(goodrules4)
inspect(goodrules4[x4])


###########################################################################
# A rules for metrics
colnames(USdata)
data2 <- subset(USdata[,c(34:41)])
data2 <-  replace(data2,TRUE,lapply(data2, factor))
str(data2)

rules4 <- apriori(data2,parameter = list(support=0.6, confidence= 0.6))
summary(rules4)
plot(rules4)
inspect(rules4)

goodrules4 <- rules4[quality(rules4)$lift>1]
inspect(goodrules4)
summary(goodrules4)
plot(goodrules4,method= "graph",measure="support",shading="lift",engine= "interactive")

#3
x4 <- is.maximal(goodrules4)
inspect(goodrules4[x4])


############################################################################
#plotting
library(sp)
install.packages("rworldmap")
library(rworldmap)

df <- data.frame(countryNames,Detractors,NPS)
df$countryNames <- gsub("([a-z])([A-Z])", "\\1 \\2",df$countryNames)
Worldmap <- joinCountryData2Map(df, joinCode="NAME",
                                nameJoinColumn="countryNames")
mapCountryData(Worldmap, nameColumnToPlot="NPS",colourPalette = "topo",catMethod="fixedWidth",addLegend = TRUE,borderCol = "grey",mapTitle = "NPS distribution across countries")
?mapCountryData


Worldmap2 <- joinCountryData2Map(df, joinCode="NAME",
                                 nameJoinColumn="countryNames")
mapCountryData(Worldmap2, nameColumnToPlot="Detractors",colourPalette = "heat",catMethod="fixedWidth",addLegend = TRUE,borderCol = "black",mapTitle = "Percentage of Detractors across countries")
