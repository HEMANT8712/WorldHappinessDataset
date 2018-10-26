
########################################################## 
# 
# Final Project--STAT5703-HEMANT GUPTA-101062246 
#                         ANURAG DAS-101089268 
#                         Manoj Kakarla-101071887
# 
# 
########################################################## 
#packages: 
install.packages("rggobi") 
install.packages("cluster") 
install.packages("stats") 
install.packages("fpc") 
install.packages("flexclust") 
install.packages("plyr") 
install.packages("fastICA") 

install.packages("rpart")
install.packages("caret")
install.packages("rpart.plot")
install.packages("MASS")
install.packages("DAAG")
install.packages("tree")


## Setting the Path of Directory

## Setting the Path of Directory

drive="C:"
path.upto <- paste("STAT5703-Final Project-Hemant Anurag Manoj", sep="/" )
code.dir <- paste(drive, path.upto,"Code", sep="/")
data.dir <- paste(drive, path.upto,"Data", sep="/")
work.dir <- paste(drive, path.upto,"Work", sep="/")
setwd(work.dir)

##Reading the CSV format of Cars Data
happiness.file<- paste(data.dir,"AllCombined.csv",sep="/")

happiness.data <- read.csv(happiness.file, header=TRUE)

head(happiness.data)

##########################################################
#
#               FUNCTIONS
#
##########################################################

#
##Function: Not Contain in the Set 
#

"%w/o%"<- function(x,y) x[!x %in% y]

#
## Function Set the indices for the training/test sets
#
get.train <- function (data.sz, train.sz)
{
  set.seed(123)
  # Take subsets of data for training/test samples
  # Return the indices
  train.ind <- sample(data.sz, train.sz)
  test.ind <- (data.sz) %w/o% train.ind
  list(train=train.ind, test=test.ind)
}

#
## Function Set the indices for the  sets
#
get.subset <- function (data, size)
{
  set.seed(123)
  # Take subsets of data for training/test samples
  # Return the indices
  data_subset <- sample(data, size)
  
}
#
#Function:========DBI Function=======#
#
Davies.Bouldin <- function(A, SS, m) {
  # A  - the centres of the clusters
  # SS - the within sum of squares
  # m  - the sizes of the clusters
  N <- nrow(A)   # number of clusters
  # intercluster distance
  S <- sqrt(SS/m)
  # Get the distances between centres
  M <- as.matrix(dist(A))
  # Get the ratio of intercluster/centre.dist
  R <- matrix(0, N, N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      R[i,j] <- (S[i] + S[j])/M[i,j]
      R[j,i] <- R[i,j]
    }
  }
  return(mean(apply(R, 1, max)))
}

#
#Function for Data Standardization
#

#################################
# Standardize data
####################################
f.data.std <- function(data) {
  data <- as.matrix(data)
  bar <- apply(data, 2, mean)
  s <- apply(data, 2, sd)
  t((t(data) - bar)/s)
}

#####################################
#WhitenedData: Centre and sphere data
#########################################
Sphere.Data <- function(data) {
  data <- as.matrix(data)
  data <- t(t(data) - apply(data, 2, mean))
  data.svd <- svd(var(data))
  sphere.mat <- t(data.svd$v %*% (t(data.svd$u) * (1/sqrt(data.svd$d))))
  return(data %*% sphere.mat)
}


##############################
#
#FUNCTION: error_cal- For calculating wrong classification in Cluster
#
############################
error_cal <- function(tbl,cluster_size)
{
  wrong_data <- 0
  for(clust in 1:cluster_size)
  {
    wrong_data <- wrong_data + (sum(tbl[,clust])-max(tbl[,clust]))
  }
  return(wrong_data)
}
####################################################
#Function: Clustering_euclidean
# Cluster Size Varies 2 to 15
# SSE and DBI is determined at each iteration
####################################################

clustering_euclidean <- function(data_set,data_set.orig, limit) 
{
  # set.seed(654321)
  oldpar <- par(mfrow = c(4,4))
  par(mar=c(2,1,2,1))
  
  errs <- rep(0, 7)
  DBI <- rep(0, 7)
  ##Package(cluster)
  library(cluster)
  library(stats)
  library(fpc)
  library(flexclust)
  
  
  #Loop for different CLuster Size
  for (i in limit)
  {
    min_error <- 450
    min_error_km <- 0
    best.seed <- 0
    #Loop for Seed
    for (j in 2:1000)
    {  
      
      set.seed(j)
      #Clustering Using K means
      KM <- kmeans((data_set[,]), i, 25)
      ct.km <- table(KM$cluster, KM$cluster)
      
      #Calculating toal wrong data for each seed         
      error <- error_cal(ct.km,i)
      if(min_error > error)
      {
        #Storing Error count and Kmeans output and best seed for min error
        min_error <- error
        min_error_km <-KM
        best.seed <- j
      }  
      
    }     
    print(paste("Best Seed for Cluster Size " , i ,"is " , best.seed))
    
    print(paste("Total Wrong in Cluster Size " , i ,"is " , min_error))
    
    print(paste("Centroids for Cluster Size " , i ,"are :"))
    
    print(min_error_km$centers)
    
    print(min_error_km)
    
    #Plotting the CLuster
    plotcluster(data_set, col=min_error_km$cluster,min_error_km$cluster, main=paste(i,"clusters-Eucli"))
    
    if(length(limit) > 1)
    {  
      #CLuster Analysis
      errs[i-1] <- sum(min_error_km$withinss)
      DBI[i-1] <- Davies.Bouldin(min_error_km$centers, min_error_km$withinss, min_error_km$size)
      
     
    }
    
  }
  if(length(limit) > 1)
  {
    plot(2:15, errs, main = "SSE")
    lines(2:15, errs)
    #
    plot(2:15, DBI, main = "Davies-Bouldin")
    lines(2:15, DBI)
    #
  }
  else
  {  
    return(min_error_km)
  }
  return(errs)
}



########################################################################################
   #                 Training and Test Dataset

########################################################################################
#Training Dataset of year 2015 and 2016
happiness.train <- happiness.data[happiness.data$Year!=2017,]

#Test Dataset
##Reading the CSV format of Cars Data
happiness.test <- happiness.data[happiness.data$Year==2017,]


out = c (nrow(happiness.train), nrow(happiness.test))

##Setting X axis name
x.names=c("Train","Test")

barplot(out, main="Data Spliting",xaxt="n")
axis(1,at = 1:2,labels=x.names)

########################################################################################

#######                    DATA VISUALIZATION

########################################################################################
happiness.15.16<- happiness.train
happiness.15.16$Region = as.factor(happiness.15.16$Region)

happiness.15.16$Region<-as.numeric(happiness.train$Region)

#happiness.15.16a <- as.data.frame(happiness.15.16)

#happiness.15.16[is.na(happiness.15.16$Region),2] <- 0
## Factoring the Car Data ORigin Column


happiness.15.16 <- happiness.15.16[order(happiness.15.16$Region),]


happiness.col <- happiness.15.16$Region

## Cleaning the data for Column Year
#for (j in 1:length(happiness.15.16$Region)){
#  if(setequal(happiness.15.16[j,2],"North America")){
#    happiness.15.16[j,2] <- 2
#  }
#  else
#  {}
#}

#> unique(happiness.15.16$Region)
#[1] Australia and New Zealand       Central and Eastern Europe     
#[3] Eastern Asia                    Latin America and Caribbean    
#[5] Middle East and Northern Africa North America                  
#[7] Southeastern Asia               Southern Asia                  
#[9] Sub-Saharan Africa              Western Europe              
## Ggobi PLot : Scatterplot and Parallel COordinate gives overall picture of 
###              data but its hard to scale it.

library(rggobi)
g <- ggobi(happiness.15.16)

display(g[1], "Scatterplot Matrix")

display(g[1], "Parallel Coordinates Display")

## Plotting the Cars.dat with respect to origin Pairs Plot as it also contain scale for each parameter


### CO-PLOT :-To get the concentation of points with yearly for each origin

par(mfrow=c(2,4))

plot(happiness.data$Year,happiness.data$Economy..GDP.per.Capita., col= happiness.data$Year)
plot(happiness.data$Year,happiness.data$Family, col= happiness.data$Year)
plot(happiness.data$Year,happiness.data$Health..Life.Expectancy., col= happiness.data$Year)
plot(happiness.data$Year,happiness.data$Freedom, col= happiness.data$Year)
plot(happiness.data$Year,happiness.data$Trust..Government.Corruption., col= happiness.data$Year)
plot(happiness.data$Year,happiness.data$Generosity, col= happiness.data$Year)
plot(happiness.data$Year,happiness.data$Dystopia.Residual, col= happiness.data$Year)


#pairs(happiness.15.16[,3:12], col=happiness.col)


#Categorize data by Region
#Region - Australia and Newzealand
happiness.AN <- happiness.15.16[happiness.15.16$Region==1,-1]
#Region- Central_Eastern_Europe
happiness.CEE <-happiness.15.16[happiness.15.16$Region==2,-1]
#Region-Eastern_Asia
happiness.EA <-happiness.15.16[happiness.15.16$Region==3,-1]
#Region-Latin_America_Caribbean
happiness.LAC <- happiness.15.16[happiness.15.16$Region==4,-1]
#Region-Middle_East_Northern_Africa
happiness.MENA <- happiness.15.16[happiness.15.16$Region==5,-1]
#Region-North_America\
happiness.NA <- happiness.15.16[happiness.15.16$Region==7,-1]
#Region-Sourtheastern_Asia
happiness.SEA <- happiness.15.16[happiness.15.16$Region==8,-1]
#Region-Sourthern_Asia
happiness.SA <- happiness.15.16[happiness.15.16$Region==9,-1]
#Region-Sub_Saharan_Africa
happiness.SSA <- happiness.15.16[happiness.15.16$Region==10,-1]
#Region-Western_Europe
happiness.WE <- happiness.15.16[happiness.15.16$Region==11,-1]


par(mfrow=c(2,5))

##Setting X axis name
x.names=c("AN","CEE","EA","LAC","MENA","NA","SEA","SA","SSA","WE")

## PLotting the Box Plot of All 6 parameter based on Origin and keeping na.rm =TRUE to avoid NA values
for(i in 2:10){
  plot.title=colnames(happiness.15.16)[i+1]
  boxplot(happiness.AN[,i],happiness.CEE[,i],happiness.EA[,i],happiness.LAC[,i],happiness.MENA[,i],happiness.NA[,i],happiness.SEA[,i],happiness.SA[,i],happiness.SSA[,i],happiness.WE[,i],main=plot.title,xaxt="n", col= heat.colors(10), na.rm=TRUE)
  axis(1,at = 1:10,labels=x.names)
}



#####We use Box plot to get the Aveage value of the data per year
#####Coplot shows the distribition of data point but clearly didnot mention
#####average avlue concentartion . This helps in deduction Accurately.

#Categorize data by Year
#Region - Australia and Newzealand
happiness.AN.15=happiness.AN[happiness.AN$Year==2015,-1]
happiness.AN.16=happiness.AN[happiness.AN$Year==2016,-1]

##BoxPlot of data based on Yearly variation and keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.AN)[i+1], "A&N",sep = '_')
  boxplot(happiness.AN.15[,i],happiness.AN.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}


#Region - Central Eastern Europe
happiness.CEE.15=happiness.CEE[happiness.CEE$Year==2015,-1]
happiness.CEE.16=happiness.CEE[happiness.CEE$Year==2016,-1]

##BoxPlot of data based on Yearly variation CEEd keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.CEE)[i+1], "CEE",sep = '_')
  boxplot(happiness.CEE.15[,i],happiness.CEE.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}

#Region - Eastern ASia
happiness.EA.15=happiness.EA[happiness.EA$Year==2015,-1]
happiness.EA.16=happiness.EA[happiness.EA$Year==2016,-1]

##BoxPlot of data based on Yearly variation EAd keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.EA)[i+1], "EA",sep = '_')
  boxplot(happiness.EA.15[,i],happiness.EA.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}

#Region - Latin America and Carribean
happiness.LAC.15=happiness.LAC[happiness.LAC$Year==2015,-1]
happiness.LAC.16=happiness.LAC[happiness.LAC$Year==2016,-1]

##BoxPlot of data based on Yearly variation ead keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.LAC)[i+1], "LA&C",sep = '_')
  boxplot(happiness.LAC.15[,i],happiness.LAC.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}

#Region - Middle East and NOrthern Africa
happiness.MENA.15=happiness.MENA[happiness.MENA$Year==2015,-1]
happiness.MENA.16=happiness.MENA[happiness.MENA$Year==2016,-1]

##BoxPlot of data based on Yearly variation and keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.MENA)[i+1], "ME& NA",sep = '_')
  boxplot(happiness.MENA.15[,i],happiness.MENA.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}


#Region - North America
happiness.NA.15=happiness.NA[happiness.NA$Year==2015,-1]
happiness.NA.16=happiness.NA[happiness.NA$Year==2016,-1]

##BoxPlot of data based on Yearly variation and keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.NA)[i+1], "NA",sep = '_')
  boxplot(happiness.NA.15[,i],happiness.NA.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}

#Region - South Eastern Asia
happiness.SEA.15=happiness.SEA[happiness.SEA$Year==2015,-1]
happiness.SEA.16=happiness.SEA[happiness.SEA$Year==2016,-1]

##BoxPlot of data based on Yearly variation and keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.SEA)[i+1], "SEA",sep = '_')
  boxplot(happiness.SEA.15[,i],happiness.SEA.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}

#Region - Southern Asia
happiness.SA.15=happiness.SA[happiness.SA$Year==2015,-1]
happiness.SA.16=happiness.SA[happiness.SA$Year==2016,-1]

##BoxPlot of data based on Yearly variation and keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.SA)[i+1], "SA",sep = '_')
  boxplot(happiness.SA.15[,i],happiness.SA.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}

#Region - Sub Saharan Africa
happiness.SSA.15=happiness.SSA[happiness.SSA$Year==2015,-1]
happiness.SSA.16=happiness.SSA[happiness.SSA$Year==2016,-1]

##BoxPlot of data based on Yearly variation and keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.SSA)[i+1], "SSA",sep = '_')
  boxplot(happiness.SSA.15[,i],happiness.SSA.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}

#Region - Western Europe
happiness.WE.15=happiness.WE[happiness.WE$Year==2015,-1]
happiness.WE.16=happiness.WE[happiness.WE$Year==2016,-1]

##BoxPlot of data based on Yearly variation and keeping na.rm =TRUE to avoid NA values
par(mfrow=c(3,3))
x.names=c("2015","2016")
for(i in 1:9){
  
  plot.title <- paste(colnames(happiness.WE)[i+1], "WE",sep = '_')
  boxplot(happiness.WE.15[,i],happiness.WE.16[,i],main=plot.title, xaxt="n", col= heat.colors(2),na.rm=TRUE)
  axis(1,at = 1:2,labels=x.names)
}


#####################################################################################
# DIMENSION REDUCTION-PCA

#####################################################################################


happiness.std <- f.data.std(happiness.data[5:11])
# Get principal component vectors using prcomp 
pc.happiness <- prcomp(happiness.std[,])
summary(pc.happiness)
plot(pc.happiness)
# First  principal components
happiness.pc <- data.frame(pc.happiness$x[,1:3])
head(happiness.pc)



happiness.train.std<- f.data.std(happiness.train[5:11])
# Get principal component vectors using prcomp 
pc.train.std <- prcomp(happiness.train.std)
summary(pc.train.std)
plot(pc.train.std)

# First  principal components
happiness.train.pc <- data.frame(pc.train.std$x[,1:3])
head(happiness.train.pc)


happiness.test.std<- f.data.std(happiness.test[5:11])
# Get principal component vectors using prcomp 
pc.test.std <- prcomp(happiness.test.std)
summary(pc.test.std)
plot(pc.test.std)

# First  principal components
happiness.test.pc <- data.frame(pc.test.std$x[,1:3])
head(happiness.test.pc)





#########################################################################################

##########################DIMENSION REDUCTION-ICA########################################################

#########################################################################################
library(fastICA)


#Whitening the whole happiness dataset
happiness.white <- Sphere.Data(happiness.data[,5:11])

##Taking number of components as 3 as want to compare it with result of PCA 3 components

happiness.white.ica <- fastICA(happiness.white, 3, alg.typ = "parallel", fun = "logcosh", alpha = 1,
                    method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose =
                      TRUE)

happiness.ica<-happiness.white.ica$S
#Estimated Source Matrix
head(happiness.ica)



#Whitening the Train Happiness dataset
happiness.train.white <- Sphere.Data(happiness.train[,5:11])

##Taking number of components as 3 as want to compare it with result of PCA 3 components
happiness.train.white.ica <- fastICA(happiness.train.white, 3, alg.typ = "parallel", fun = "logcosh", alpha = 1,
                         method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose =
                           TRUE)
#Estimated Source Matrix
happiness.train.ica<-happiness.train.white.ica$S
head(happiness.train.ica)




#Whitening the Test Happiness dataset
happiness.test.white <- Sphere.Data(happiness.test[,5:11])

##Taking number of components as 3 as want to compare it with result of PCA 3 components
happiness.test.white.ica <- fastICA(happiness.test.white, 3, alg.typ = "parallel", fun = "logcosh", alpha = 1,
                               method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose =
                                 TRUE)
#Estimated Source Matrix
happiness.test.ica<-happiness.test.white.ica$S
head(happiness.test.ica)



########################################################################################

#  DATA REDUCTION

#######################################################################################

#################################RAW DATASET#################################
cluster_range <- 2:15

##Clustering on train and test Raw DataSet 
happiness.ret <- clustering_euclidean(happiness.data[5:11],happiness.data, cluster_range)

#Optimal Cluster Size is 3
# Create a K-Means cluster with 3 groups based on cols 5:11
# (GDP, Family, Life Expectancy, Freedom, Generosity, Corruption, Dystopia Residual)
km <- clustering_euclidean(happiness.data[5:11],happiness.data, 3)

g1<-happiness.data[km$cluster == 1,]$Happiness.Score
g2<-happiness.data[km$cluster == 2,]$Happiness.Score
g3<-happiness.data[km$cluster == 3,]$Happiness.Score
# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of Raw Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5)) )


#Who is in the Happiest and least happiest Group?
top <- which.max(c(mean(g1),mean(g2), mean(g3))) # which is the top group
bottom <- which.min(c(mean(g1),mean(g2), mean(g3))) # which is the top group


happiest <- happiness.data[km$cluster == top, c(1,4)]
least_happiest <- happiness.data[km$cluster == bottom, c(1,4)]

#Cluster Size
print(km$size)

print(paste("RAW-Happiest Group is g", top, sep=""))

print(paste("RAW-Least Happiest Group is g", bottom, sep=""))

#Countries with highest and least Happiness score in clusters
head(happiest[order(happiest$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest[order(least_happiest$Happiness.Score, decreasing=TRUE), ])

#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries

print("Median Happiness Score for Whole Dataset-")
(median(happiness.data$Happiness.Score))

head(least_happiest[least_happiest$Happiness.Score > median(happiness.data$Happiness.Score), ])
head(happiest[happiest$Happiness.Score < median(happiness.data$Happiness.Score), ])


print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest[happiest$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest[happiest$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest[least_happiest$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest[least_happiest$Happiness.Score < median(happiness.data$Happiness.Score), ])


#################################Standard DATASET#################################
cluster_range <- 2:15

##Clustering on train and test Raw DataSet 
happiness.std.ret <- clustering_euclidean(happiness.std,happiness.data, cluster_range)

#Optimal Cluster Size is 3
# Create a K-Means cluster with 3 groups based on cols 5:11
# (GDP, Family, Life Expectancy, Freedom, Generosity, Corruption, Dystopia Residual)
km.std <- clustering_euclidean(happiness.std,happiness.data, 3)

g1.std<-happiness.data[km.std$cluster == 1,]$Happiness.Score
g2.std<-happiness.data[km.std$cluster == 2,]$Happiness.Score
g3.std<-happiness.data[km.std$cluster == 3,]$Happiness.Score
# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.std, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of STD. Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.std, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.std, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5)) )


#Who is in the Happiest and least happiest Group?
top <- which.max(c(mean(g1.std),mean(g2.std), mean(g3.std))) # which is the top group
bottom <- which.min(c(mean(g1.std),mean(g2.std), mean(g3.std))) # which is the top group


happiest.std <- happiness.data[km.std$cluster == top, c(1,4)]
least_happiest.std <- happiness.data[km.std$cluster == bottom, c(1,4)]

#Cluster Size
print(km.std$size)

print(paste("Standard-Happiest Group is g", top, sep=""))

print(paste("Standard-Least Happiest Group is g", bottom, sep=""))

#Countries with highest and least Happiness score in clusters
head(happiest.std[order(happiest.std$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.std[order(least_happiest.std$Happiness.Score, decreasing=TRUE), ])

#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries

print("Median Happiness Score for Whole Dataset-")
(median(happiness.data$Happiness.Score))

head(least_happiest.std[least_happiest.std$Happiness.Score > median(happiness.data$Happiness.Score), ])
head(happiest.std[happiest.std$Happiness.Score < median(happiness.data$Happiness.Score), ])


print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.std[happiest.std$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.std[happiest.std$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.std[least_happiest.std$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.std[least_happiest.std$Happiness.Score < median(happiness.data$Happiness.Score), ])

#################################Whitened DATASET#################################
cluster_range <- 2:15

##Clustering on train and test Raw DataSet 
happiness.white.ret <- clustering_euclidean(happiness.white,happiness.data, cluster_range)

#Optimal Cluster Size is 6
# Create a K-Means cluster with 3 groups based on cols 5:11
# (GDP, Family, Life Expectancy, Freedom, Generosity, Corruption, Dystopia Residual)
km.white <- clustering_euclidean(happiness.white,happiness.data, 6)

g1.white<-happiness.data[km.white$cluster == 1,]$Happiness.Score
g2.white<-happiness.data[km.white$cluster == 2,]$Happiness.Score
g3.white<-happiness.data[km.white$cluster == 3,]$Happiness.Score
g4.white<-happiness.data[km.white$cluster == 4,]$Happiness.Score
g5.white<-happiness.data[km.white$cluster == 5,]$Happiness.Score
g6.white<-happiness.data[km.white$cluster == 6,]$Happiness.Score

# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.white, xlim=c(0,10), col=rgb(0.5,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of white. Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.white, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.white, xlim=c(0,10), col=rgb(0,0.5,1,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g4.white, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g5.white, xlim=c(0,10), col=rgb(0,0,0.5,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g6.white, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)

legend("topleft", c("Group1", "Group2", "Group3","Group4","Group5","Group6")
       , fill=c(rgb(0.5,0,0,0.5),rgb(1,0,0,0.5),rgb(0,0.5,0,0.5),rgb(0,1,0,0.5),rgb(0,0,0.5,0.5),rgb(0,0,1,0.5)) )


#Who is in the Happiest and least happiest Group?
top <- which.max(c(mean(g1.white),mean(g2.white), mean(g3.white),mean(g4.white),mean(g5.white),mean(g6.white))) # which is the top group
bottom <- which.min(c(mean(g1.white),mean(g2.white), mean(g3.white),mean(g4.white),mean(g5.white),mean(g6.white))) # which is the top group


happiest.white <- happiness.data[km.white$cluster == top, c(1,4)]
least_happiest.white <- happiness.data[km.white$cluster == bottom, c(1,4)]

#Cluster Size
print(km.white$size)

print(paste("Whitened-Happiest Group is g", top, sep=""))

print(paste("Whitened-Least Happiest Group is g", bottom, sep=""))

#Countries with highest and least Happiness score in clusters
head(happiest.white[order(happiest.white$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.white[order(least_happiest.white$Happiness.Score, decreasing=TRUE), ])

#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries

print("Median Happiness Score for Whole Dataset-")
(median(happiness.data$Happiness.Score))

head(least_happiest.white[least_happiest.white$Happiness.Score > median(happiness.data$Happiness.Score), ])
head(happiest.white[happiest.white$Happiness.Score < median(happiness.data$Happiness.Score), ])


print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.white[happiest.white$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.white[happiest.white$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.white[least_happiest.white$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.white[least_happiest.white$Happiness.Score < median(happiness.data$Happiness.Score), ])
#################Standard Data Reduction############
happiness.std.new <- happiness.data
happiness.std.new$cluster <- km.std$cluster

happiness.std.new <- happiness.data
happiness.std.new$cluster <- km.std$cluster

med1<-median(happiness.std.new$Happiness.Score[happiness.std.new$cluster == 1])
med2<-median(happiness.std.new$Happiness.Score[happiness.std.new$cluster == 2])
med3<-median(happiness.std.new$Happiness.Score[happiness.std.new$cluster == 3])
print(paste("Median Values-> Cluster1 = ", med1,", Cluster2 = ",med2,", Cluster3 = ",med3, sep=""))

T1.std.ind <- which(happiness.std.new$cluster == 1)
T2.std.ind <- which(happiness.std.new$cluster == 2)
T3.std.ind <- which(happiness.std.new$cluster == 3)

T1.std.size <- round((2*length(T1.std.ind))/3)
T2.std.size <- round((2*length(T2.std.ind))/3)
T3.std.size <- round((2*length(T3.std.ind))/3)

T1.std <- get.subset(T1.std.ind,T1.std.size)
T2.std <- get.subset(T2.std.ind,T2.std.size)
T3.std <- get.subset(T3.std.ind,T3.std.size)

happiness.std.reduced.ind <- c(T1.std,T2.std,T3.std) 

happiness.std.reduced.data <- happiness.data[happiness.std.reduced.ind,]

head(happiness.std.reduced.data)

##DIvide this Reduced Dataset in Training and  Test Based on Year and Its proper 2/3 train and 
##1/3 Test Ratio
nrow(happiness.std.reduced.data[happiness.std.reduced.data$Year!=2017,])
nrow(happiness.std.reduced.data[happiness.std.reduced.data$Year==2017,])

#Train Dataset with year 2015 and 2016- Total Rows 205
happiness.std.reduced.train <- happiness.std.reduced.data[happiness.std.reduced.data$Year!=2017,]
#Test Dataset with Year 2017 with Total Rows 108##Setting X axis name

happiness.std.reduced.test <- happiness.std.reduced.data[happiness.std.reduced.data$Year==2017,]

out.std.reduced = c (nrow(happiness.std.reduced.data),nrow(happiness.std.reduced.train), nrow(happiness.std.reduced.test))

x.names=c("Complete","Train","Test")

barplot(out.std.reduced,main="Reduced Std. Data Spliting",xaxt="n",width=c(1,1,1))
axis(1,at = 1:3,labels=x.names)




########################################################################################

# UnSupervised Learning

#######################################################################################
########################################################################################

# UnSupervised Learning

#######################################################################################
##################################################################################

## Train Raw data 

##################################################################################

cluster_range <- 2:15

##Clustering on train and test Raw DataSet 
happiness.train.ret <- clustering_euclidean(happiness.train[5:11],happiness.train, cluster_range)

#Optimal Cluster Size is 3
# Create a K-Means cluster with 3 groups based on cols 5:10
# (GDP, Social Support, Life Expectancy, Freedom, Generosity, Corruption)
km.train <- clustering_euclidean(happiness.train[, 5:11],happiness.train,3)

g1.train<-happiness.train[km.train$cluster == 1,]$Happiness.Score
g2.train<-happiness.train[km.train$cluster == 2,]$Happiness.Score
g3.train<-happiness.train[km.train$cluster == 3,]$Happiness.Score
# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.train, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of Train Raw Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.train, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.train, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5)) )

#Who is in the Happiest and Least happy Group?
top <- which.max(c(mean(g1.train),mean(g2.train), mean(g3.train))) # which is the top group
bottom <- which.min(c(mean(g1.train),mean(g2.train), mean(g3.train))) # which is the top group

happiest.train <- happiness.train[km.train$cluster == top, c(1,4)]
least_happiest.train <- happiness.train[km.train$cluster == bottom, c(1,4)]
#Cluster Size
print(km.train$size)

print(paste("RAW-Train Happiest Group is g", top, sep=""))

print(paste("RAW-Train-Least Happiest Group is g", bottom, sep=""))


head(happiest.train[order(happiest.train$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.train[order(least_happiest.train$Happiness.Score, decreasing=TRUE), ])

#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries

print("Median Happiness Score for training Dataset-")
(median(happiness.train$Happiness.Score))

head(least_happiest.train[least_happiest.train$Happiness.Score > median(happiness.train$Happiness.Score), ])
head(happiest.train[happiest.train$Happiness.Score < median(happiness.train$Happiness.Score), ])

print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.train[happiest.train$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.train[happiest.train$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.train[least_happiest.train$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.train[least_happiest.train$Happiness.Score < median(happiness.data$Happiness.Score), ])


###########TEST Raw Dataset######################################

#Optimal Cluster Size is 3
# Create a K-Means cluster with 3 groups based on cols 5:10
# (GDP, Social Support, Life Expectancy, Freedom, Generosity, Corruption)
km.test <- clustering_euclidean(happiness.test[, 5:11],happiness.test,3)

g1.test<-happiness.test[km.test$cluster == 1,]$Happiness.Score
g2.test<-happiness.test[km.test$cluster == 2,]$Happiness.Score
g3.test<-happiness.test[km.test$cluster == 3,]$Happiness.Score
# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.test, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of Test Raw Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.test, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.test, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5)) )

#Who is in the Happiest Group?
top <- which.max(c(mean(g1.test),mean(g2.test), mean(g3.test))) # which is the top group
bottom <- which.min(c(mean(g1.test),mean(g2.test), mean(g3.test))) # which is the top group

happiest.test <- happiness.test[km.test$cluster == top, c(1,4)]
least_happiest.test <- happiness.test[km.test$cluster == bottom, c(1,4)]


#Cluster Size
print(km.test$size)

print(paste("RAW-Test Happiest Group is g", top, sep=""))

print(paste("RAW-Test Least Happiest Group is g", bottom, sep=""))

head(happiest.test[order(happiest.test$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.test[order(least_happiest.test$Happiness.Score, decreasing=TRUE), ])

#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries

print("Median Happiness Score for testing Dataset-")
(median(happiness.test$Happiness.Score))

head(least_happiest.test[least_happiest.test$Happiness.Score > median(happiness.test$Happiness.Score), ])
head(happiest.test[happiest.test$Happiness.Score < median(happiness.test$Happiness.Score), ])

print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.test[happiest.test$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.test[happiest.test$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.test[least_happiest.test$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.test[least_happiest.test$Happiness.Score < median(happiness.data$Happiness.Score), ])

##############################PCA######################################################

cluster_range <- 2:15
##Clustering on train and test Raw pcSet 
happiness.pc.ret <- clustering_euclidean(happiness.pc,happiness.data, cluster_range)

#Optimal Cluster Size is 3
# Create a K-Means cluster with 3 groups based on cols 5:10
# (GDP, Family, Life Expectancy, Freedom, Generosity, Corruption, Dystopia Residual)
km.pc <- clustering_euclidean(happiness.pc,happiness.data, 3)

g1.pc<-happiness.data[km.pc$cluster == 1,]$Happiness.Score
g2.pc<-happiness.data[km.pc$cluster == 2,]$Happiness.Score
g3.pc<-happiness.data[km.pc$cluster == 3,]$Happiness.Score
# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.pc, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of PCA Whole Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.pc, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.pc, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5)) )

#Cluster Size
print(km.pc$size)

#Who is in the Happiest Group?
top <- which.max(c(mean(g1.pc),mean(g2.pc), mean(g3.pc))) # which is the top group
bottom <- which.min(c(mean(g1.pc),mean(g2.pc), mean(g3.pc))) # which is the top group

happiest.pc <- happiness.data[km.pc$cluster == top, c(1,4)]
least_happiest.pc <- happiness.data[km.pc$cluster == bottom, c(1,4)]

print(paste("PCA-Happiest Group is g", top, sep=""))

print(paste("PCA-Least Happiest Group is g", bottom, sep=""))

head(happiest.pc[order(happiest.pc$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.pc[order(least_happiest.pc$Happiness.Score, decreasing=TRUE), ])
#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries
print("Median Happiness Score for Whole Dataset-")
(median(happiness.data$Happiness.Score))

head(least_happiest.pc[least_happiest.pc$Happiness.Score > median(happiness.data$Happiness.Score), ])
head(happiest.pc[happiest.pc$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.pc[happiest.pc$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.pc[happiest.pc$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.pc[least_happiest.pc$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.pc[least_happiest.pc$Happiness.Score < median(happiness.data$Happiness.Score), ])


######################Train Dataset-After PCA###################


cluster_range <- 2:15
##Clustering on train and test Raw pcSet 
happiest.train.pc.ret <- clustering_euclidean(happiness.train.pc,happiness.train, cluster_range)

#Optimal Cluster Size is 3
# Create a K-Means cluster with 3 groups based on cols 5:10
# (GDP, Family, Life Expectancy, Freedom, Generosity, Corruption,Dystopia Residual)
km.train.pc <- clustering_euclidean(happiness.train.pc,happiness.train, 3)

g1.train.pc<-happiness.train[km.train.pc$cluster == 1,]$Happiness.Score
g2.train.pc<-happiness.train[km.train.pc$cluster == 2,]$Happiness.Score
g3.train.pc<-happiness.train[km.train.pc$cluster == 3,]$Happiness.Score
# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.train.pc, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of Train PCA-Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.train.pc, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.train.pc, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5)) )

#Who is in the Happiest and Least Happiest group Group?
top <- which.max(c(mean(g1.train.pc),mean(g2.train.pc), mean(g3.train.pc))) # which is the top group
bottom <- which.min(c(mean(g1.train.pc),mean(g2.train.pc), mean(g3.train.pc))) # which is the top group

happiest.train.pc <- happiness.train[km.train.pc$cluster == top, c(1,4)]
least_happiest.train.pc <- happiness.train[km.train.pc$cluster == bottom, c(1,4)]

#Cluster Size
print(km.train.pc$size)

print(paste("Train-PCA-Happiest Group is g", top, sep=""))

print(paste("Train-PCA-Least Happiest Group is g", bottom, sep=""))


head(happiest.train.pc[order(happiest.train.pc$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.train.pc[order(least_happiest.train.pc$Happiness.Score, decreasing=TRUE), ])
#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries
print("Median Happiness Score for testing Dataset-")
(median(happiness.train$Happiness.Score))

head(least_happiest.train.pc[least_happiest.train.pc$Happiness.Score > median(happiness.train$Happiness.Score), ])
head(happiest.train.pc[happiest.train.pc$Happiness.Score < median(happiness.train$Happiness.Score), ])


print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.train.pc[happiest.train.pc$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.train.pc[happiest.train.pc$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.train.pc[least_happiest.train.pc$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.train.pc[least_happiest.train.pc$Happiness.Score < median(happiness.data$Happiness.Score), ])

##########TEST DATASET AFTER PCA###########################################

#Optimal Cluster Size is 3
# Create a K-Means cluster with 3 groups based on cols 5:10
# (GDP, Social Support, Life Expectancy, Freedom, Generosity, Corruption)
km.test.pc <- clustering_euclidean(happiness.test.pc,happiness.test,3)

g1.test.pc<-happiness.test[km.test.pc$cluster == 1,]$Happiness.Score
g2.test.pc<-happiness.test[km.test.pc$cluster == 2,]$Happiness.Score
g3.test.pc<-happiness.test[km.test.pc$cluster == 3,]$Happiness.Score
# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.test.pc, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of Test PCA-Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.test.pc, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.test.pc, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5)) )

#Who is in the Happiest and least Happiest Group?
top <- which.max(c(mean(g1.test.pc),mean(g2.test.pc), mean(g3.test.pc))) # which is the top group
bottom <- which.min(c(mean(g1.test.pc),mean(g2.test.pc), mean(g3.test.pc))) # which is the top group

happiest.test.pc <- happiness.test[km.test.pc$cluster == top, c(1,4)]
least_happiest.test.pc <- happiness.test[km.test.pc$cluster == bottom, c(1,4)]

#Cluster Size
print(km.test.pc$size)

print(paste("Test-PCA-Happiest Group is g", top, sep=""))

print(paste("Test-PCA-Least Happiest Group is g", bottom, sep=""))

head(happiest.test.pc[order(happiest.test.pc$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.test.pc[order(least_happiest.test.pc$Happiness.Score, decreasing=TRUE), ])
#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries
print("Median Happiness Score for testing Dataset-")
(median(happiness.test$Happiness.Score))

head(least_happiest.test.pc[least_happiest.test.pc$Happiness.Score > median(happiness.test$Happiness.Score), ])
head(happiest.test.pc[happiest.test.pc$Happiness.Score < median(happiness.test$Happiness.Score), ])

print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.test.pc[happiest.test.pc$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.test.pc[happiest.test.pc$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.test.pc[least_happiest.test.pc$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.test.pc[least_happiest.test.pc$Happiness.Score < median(happiness.data$Happiness.Score), ])

###################################ICA######################################################

cluster_range <- 2:15
##Clustering on train and test Raw icaSet 
happiest.ica.ret <- clustering_euclidean(happiness.ica,happiness.data, cluster_range)

#Optimal Cluster Size is 4
# Create a K-Means cluster with 4 groups based on cols 5:10
# (GDP, Family, Life Expectancy, Freedom, Generosity, Corruption,Dystopia Residual)
km.ica <- clustering_euclidean(happiness.ica,happiness.data, 4)

g1.ica<-happiness.data[km.ica$cluster == 1,]$Happiness.Score
g2.ica<-happiness.data[km.ica$cluster == 2,]$Happiness.Score
g3.ica<-happiness.data[km.ica$cluster == 3,]$Happiness.Score
g4.ica<-happiness.data[km.ica$cluster == 4,]$Happiness.Score

# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.ica, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of ICA Whole Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.ica, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.ica, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g4.ica, xlim=c(0,10), col=rgb(1,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3","Group4")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,1,0.5)) )


#Who is in the Happiest and least happiest Group?
top <- which.max(c(mean(g1.ica),mean(g2.ica), mean(g3.ica),mean(g4.ica))) # which is the top group
bottom <- which.min(c(mean(g1.ica),mean(g2.ica), mean(g3.ica),mean(g4.ica))) # which is the top group

happiest.ica <- happiness.data[km.ica$cluster == top, c(1,4)]
least_happiest.ica <- happiness.data[km.ica$cluster == bottom, c(1,4)]

#Cluster Size
print(km.ica$size)

print(paste("ICA-Happiest Group is g", top, sep=""))

print(paste("ICA-Least Happiest Group is g", bottom, sep=""))

head(happiest.ica[order(happiest.ica$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.ica[order(least_happiest.ica$Happiness.Score, decreasing=TRUE), ])
#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries
print("Median Happiness Score for Whole Dataset-")
(median(happiness.data$Happiness.Score))

head(least_happiest.ica[least_happiest.ica$Happiness.Score > median(happiness.data$Happiness.Score), ])
head(happiest.ica[happiest.ica$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.ica[happiest.ica$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.ica[happiest.ica$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.ica[least_happiest.ica$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.ica[least_happiest.ica$Happiness.Score < median(happiness.data$Happiness.Score), ])


###############Train Dataset After ICA###############
###################################ICA######################################################

cluster_range <- 2:15
##Clustering on train and test Raw icaSet 
happiness.train.ica.ret <- clustering_euclidean(happiness.train.ica,happiness.train, cluster_range)

#Optimal Cluster Size is 4
# Create a K-Means cluster with 4 groups based on cols 5:10
# (GDP, Social Support, Life Expectancy, Freedom, Generosity, Corruption)
km.train.ica <- clustering_euclidean(happiness.train.ica,happiness.train, 4)

g1.train.ica<-happiness.train[km.train.ica$cluster == 1,]$Happiness.Score
g2.train.ica<-happiness.train[km.train.ica$cluster == 2,]$Happiness.Score
g3.train.ica<-happiness.train[km.train.ica$cluster == 3,]$Happiness.Score
g4.train.ica<-happiness.train[km.train.ica$cluster == 4,]$Happiness.Score

# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.train.ica, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of Train ICA Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.train.ica, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.train.ica, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g4.train.ica, xlim=c(0,10), col=rgb(1,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3","Group4")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,1,0.5)) )

#Who is in the Happiest Group?
top <- which.max(c(mean(g1.train.ica),mean(g2.train.ica), mean(g3.train.ica),mean(g4.train.ica))) # which is the top group
bottom <- which.min(c(mean(g1.train.ica),mean(g2.train.ica), mean(g3.train.ica),mean(g4.train.ica))) # which is the top group

happiest.train.ica <- happiness.train[km.train.ica$cluster == top, c(1,4)]
least_happiest.train.ica <- happiness.train[km.train.ica$cluster == bottom, c(1,4)]

#Cluster Size
print(km.train.ica$size)

print(paste("Train-ICA-Happiest Group is g", top, sep=""))

print(paste("Train-ICA-Least Happiest Group is g", bottom, sep=""))

head(happiest.train.ica[order(happiest.train.ica$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.train.ica[order(least_happiest.train.ica$Happiness.Score, decreasing=TRUE), ])
#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries
print("Median Happiness Score for train Dataset-")
(median(happiness.train$Happiness.Score))

head(least_happiest.train.ica[least_happiest.train.ica$Happiness.Score > median(happiness.train$Happiness.Score), ])
head(happiest.train.ica[happiest.train.ica$Happiness.Score < median(happiness.train$Happiness.Score), ])

print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.train.ica[happiest.train.ica$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.train.ica[happiest.train.ica$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.train.ica[least_happiest.train.ica$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.train.ica[least_happiest.train.ica$Happiness.Score < median(happiness.data$Happiness.Score), ])


########TEST DATASET AFTER ICA#################################


#Optimal Cluster Size is 4
# Create a K-Means cluster with 4 groups based on cols 5:10
# (GDP, Social Support, Life Expectancy, Freedom, Generosity, Corruption)
km.test.ica <- clustering_euclidean(happiness.test.ica, happiness.test, 4)

g1.test.ica<-happiness.test[km.test.ica$cluster == 1,]$Happiness.Score
g2.test.ica<-happiness.test[km.test.ica$cluster == 2,]$Happiness.Score
g3.test.ica<-happiness.test[km.test.ica$cluster == 3,]$Happiness.Score
g4.test.ica<-happiness.test[km.test.ica$cluster == 4,]$Happiness.Score

# plot option "col=rgb(x,x,x,0.5)"" gives fill transparency

par(mfrow=c(1,1))

hist(g1.test.ica, xlim=c(0,10), col=rgb(1,0,0,0.5), breaks=seq(0.25,10,0.25), main = "Histogram of Test ICA Dataset Happiness Score for 3 cluster-groups", xlab = "Country Happiness Score")
hist(g2.test.ica, xlim=c(0,10), col=rgb(0,1,0,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g3.test.ica, xlim=c(0,10), col=rgb(0,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
hist(g4.test.ica, xlim=c(0,10), col=rgb(1,0,1,0.5), breaks=seq(0.25,10,0.25), add=T)
legend("topleft", c("Group1", "Group2", "Group3","Group4")
       , fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,1,0.5)) )

#Who is in the Happiest Group?
top <- which.max(c(mean(g1.test.ica),mean(g2.test.ica), mean(g3.test.ica),mean(g4.test.ica))) # which is the top group
bottom <- which.min(c(mean(g1.test.ica),mean(g2.test.ica), mean(g3.test.ica),mean(g4.test.ica))) # which is the top group

happiest.test.ica <- happiness.test[km.test.ica$cluster == top, c(1,4)]
least_happiest.test.ica <- happiness.test[km.test.ica$cluster == bottom, c(1,4)]

#Cluster Size
print(km.test.ica$size)

print(paste("Test-ICA-Happiest Group is g", top, sep=""))

print(paste("Test-ICA-Least Happiest Group is g", bottom, sep=""))

head(happiest.test.ica[order(happiest.test.ica$Happiness.Score, decreasing=TRUE), ])

tail(least_happiest.test.ica[order(least_happiest.test.ica$Happiness.Score, decreasing=TRUE), ])
#Identify who in the top group has lower happiness than the 
#median happiness for the entire set of countries
print("Median Happiness Score for testing Dataset-")
(median(happiness.test$Happiness.Score))

head(least_happiest.test.ica[least_happiest.test.ica$Happiness.Score > median(happiness.test$Happiness.Score), ])
head(happiest.test.ica[happiest.test.ica$Happiness.Score < median(happiness.test$Happiness.Score), ])

print("Number of countries above and below Median Value in Happiest Group")

nrow(happiest.test.ica[happiest.test.ica$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(happiest.test.ica[happiest.test.ica$Happiness.Score < median(happiness.data$Happiness.Score), ])

print("Number of countries above and below Median Value in Least Happiest Group")

nrow(least_happiest.test.ica[least_happiest.test.ica$Happiness.Score > median(happiness.data$Happiness.Score), ])
nrow(least_happiest.test.ica[least_happiest.test.ica$Happiness.Score < median(happiness.data$Happiness.Score), ])




############################Supervised Learning####################################
############Classifications rparts#######
library(rpart)
library(caret)
library(rpart.plot)
library(MASS)
library(DAAG)
library(tree)

#####################TREE################

library(tree)

region.tree <- tree(Region~Happiness.Rank+Happiness.Score+Economy..GDP.per.Capita.+
                      Family+Health..Life.Expectancy.+Freedom+Trust..Government.Corruption.
                    +Generosity+Dystopia.Residual ,method="class",data = happiness.std.reduced.train)

plot(region.tree, type = "uniform")
text(region.tree, all= T, cex=0.7)

########Prediction and Accuracy#########
model.2<-region.tree
test<-happiness.std.reduced.test
test_pred <- predict(model.2,test)
head(test_pred)
head(test)


###############Score tree############
score.tree <- tree(Happiness.Score~Happiness.Rank+Region+Economy..GDP.per.Capita.+
                     Family+Health..Life.Expectancy.+Freedom+Trust..Government.Corruption.
                   +Generosity+Dystopia.Residual ,method="anova",data = happiness.std.reduced.train)

(score.tree)

plot(score.tree, type = "uniform")
text(score.tree, all= T, cex=0.9)
##########Residual Analysis#########
plot(residuals(score.tree), ylab="residuals")
r<-(residuals(score.tree))
sum(r^2)



########Prediction and Accuracy#########
model.2<-score.tree
test_pred <- predict(model.2,test)
Actual.Values<-round(test$Happiness.Score)
r.pred<-floor(test_pred)
Predicted.Values<-as.numeric(r.pred)
table(Predicted.Values,Actual.Values)
confusion(Predicted.Values, Actual.Values)






############# Score tree without rank###############

score.tree2 <- tree(Happiness.Score~Region+Economy..GDP.per.Capita.+
                      Family+Health..Life.Expectancy.+Freedom+Trust..Government.Corruption.
                    +Generosity+Dystopia.Residual ,method="anova",data = happiness.std.reduced.train)



(score.tree2)

plot(score.tree2, type = "uniform")
text(score.tree2, all= T, cex=0.5)



##########Residual Analysis#########
plot(residuals(score.tree2), ylab="residuals")
r<-(residuals(score.tree2))
sum(r^2)



########Prediction and Accuracy#########
model.2<-score.tree2
test_pred <- predict(model.2,test)
Actual.Values<-round(test$Happiness.Score)
r.pred<-floor(test_pred)
Predicted.Values<-as.numeric(r.pred)
table(Predicted.Values,Actual.Values)
confusion(Predicted.Values, Actual.Values)




##########Applying rpart###################
region.rp.std <- rpart(Region~Happiness.Rank+Happiness.Score+Economy..GDP.per.Capita.+
                         Family+Health..Life.Expectancy.+Freedom+Trust..Government.Corruption.
                       +Generosity+Dystopia.Residual ,method="class",data = happiness.std.reduced.train)



########Applying rpart############
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(123)
dtree_fit <- train(Happiness.Score ~., data = happiness.std.reduced.train, method = "rpart",
                   parms = list(split = "information"),
                   trControl=trctrl,
                   tuneLength = 10)



##########Plotting rpart############
prp(dtree_fit$finalModel, box.palette = "Blue", tweak = 1.2, varlen=0)


##########Residual Analysis#########
plot(residuals(dtree_fit), ylab="residuals")
r<-(residuals(dtree_fit))
sum(r^2)

########Prediction and Accuracy#########

model<-dtree_fit
test<-happiness.std.reduced.test
test_pred <- predict(model,test)
Actual.Values<-round(happiness.std.reduced.test$Happiness.Score)
r.pred<-round(test_pred)
Predicted.Values<-as.numeric(r.pred)
table(Predicted.Values,Actual.Values)
confusion(Predicted.Values, Actual.Values)






##########Analysis without rank##########
########Applying rpart############
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(123)
dtree_fit.2 <- train(Happiness.Score ~Economy..GDP.per.Capita.+
                       Family+Health..Life.Expectancy.+Freedom+Trust..Government.Corruption.
                     +Generosity, data = happiness.std.reduced.train, method = "rpart",
                     parms = list(split = "information"),
                     trControl=trctrl,
                     tuneLength = 10)

##########Plotting rpart############
prp(dtree_fit.2$finalModel, box.palette = "Green", tweak = 1.2, varlen=0)




##########Residual Analysis#########
plot(residuals(dtree_fit.2), ylab="residuals")
r<-(residuals(dtree_fit.2))
sum(r^2)


########Prediction and Accuracy#########
model.2<-dtree_fit.2
test<-happiness.std.reduced.test
test_pred <- predict(model.2,test)
Actual.Values<-round(test$Happiness.Score)
r.pred<-floor(test_pred)
Predicted.Values<-as.numeric(r.pred)
table(Predicted.Values,Actual.Values)
confusion(Predicted.Values, Actual.Values)

#######################END######################################