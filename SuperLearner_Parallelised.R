############# Super learner ###############
library(tidyverse)
library(dplyr)
library(SuperLearner)
library(foreign)
library(ggplot2)
library(ck37r) #Chris Kennedy package
library(ctmle)
#devtools::install_github("tlverse/sl3")
#library(sl3)
library(caret)
library("dbarts")
library("e1071") 
library("ranger") 
library("glmnet")
library("tmle")
library("gam")
library(parallel)
###note that the continuous SuperLearner selects the best weighted average of the included algos.

##functions ###

##computes the Q0 superlearner by exposedunexposed group
#using parallel SL which works in Windows/Mac

my.snowSL.pred=function(W,A,Y,library) {
  my.matrix <- data.frame(W)
  m0 <- snowSuperLearner(Y[A==0], my.matrix[A==0,], newX = my.matrix, SL.library = library, cluster = cluster,
                         family = binomial(),cvControl = list(shuffle = F, V=10),method="method.NNLS")
  m1 <- snowSuperLearner(Y[A==1], my.matrix[A==1,], newX = my.matrix, SL.library = library,cluster = cluster,
                         family = binomial(),cvControl = list(shuffle = F, V=10),method="method.NNLS")
  Yhat.0 <- m0$SL.predict
  Yhat.1 <- m1$SL.predict
  Yhat.A <- ifelse(A==0, Yhat.0,Yhat.1)
  
  Q.SL.object=cbind(Yhat.A,Yhat.0,Yhat.1)
  ate_SL=mean(Yhat.1-Yhat.0)
  return(list(m0=m0,m1=m1,Q.SL.object=Q.SL.object,ate_SL=ate_SL))
}


###SL  wrappers #####

##GAM
SL.gam.3 <- function (...,deg.gam = 3) {
  SL.gam (...,deg.gam = deg.gam )
}

SL.gam.2 <- function (...,deg.gam = 2) {
  SL.gam (...,deg.gam = deg.gam )
}

###RandomForest ##

#### Tuning for random Forest ####
#mtry = ## the mtry is usually set to the sqrt of the num of factors for classification
# p/3 for regression 
#ntree = 1000
#nodesize = number of  obs in terminal node, related to tree depth 
#one by one  
SL.rf.4 = function(...) {
    SL.randomForest(..., ntree = 3000, mtry=4)
  }

SL.rf.6 = function(...) {
  SL.randomForest(..., ntree = 3000, mtry=6)
}

SL.rf.8 = function(...) {
  SL.randomForest(..., ntree = 3000, mtry=8)
}

#better way using create.Learner
rf.learners = create.Learner("SL.randomForest", tune = list(ntree = c(1000,3000),mtry = c(4,6,8), detailed_names = T, name_prefix = "RF"))

#there is a tuning function in the e1071 package

#rf.tuning<-tune.randomForest(data_obs_x, y = myoutcome, data = data_obs_xy, nodesize = NULL, 
#                  mtry = NULL, ntree = NULL)
#summary(rf.tuning)
#plot(rf.tuning)
#best.randomForest(x, tunecontrol = tune.control(), ...)

#### Tuning for Ranger RF are faster than RandomForest ####
tune.ranger = list(num.trees = c(200,500),
                   mtry = c(5,10))

# Also shorten the name prefix.
rang.learners = create.Learner("SL.ranger",tune = tune.ranger, detailed_names = T, name_prefix = "Rg")

# 9 configurations
length(rang.learners$names)


#### tuning for xgboost ####
tune = list(ntrees = c(100,1000),
            max_depth = c(1,4),
            shrinkage = c(0.001, 0.1))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
xgb.learners = create.Learner("SL.xgboost",tune = tune, detailed_names = T, name_prefix = "xgb")

# 12 configurations
length(xgb.learners$names)


#############
my.library<-c(xgb.learners$names,rang.learners$names,"SL.dbarts", "SL.gam","SL.gam.2", "SL.gam.3","SL.glm","SL.glm.interaction","SL.glmnet")


# Fit the CV.SuperLearner. 




#fit using multicores#
# Setup parallel computation - use all cores on our computer.
num_cores = detectCores(all.tests = FALSE, logical = TRUE)

# How many cores does this computer have?
num_cores
options(mc.cores = num_cores)

# Check how many parallel workers we are using (on macOS/Linux).
getOption("mc.cores")

# We need to set a different type of seed that works across cores.
# Otherwise the other cores will go rogue and we won't get repeatable results.
# This version is for the "multicore" parallel system in R.
set.seed(1, "L'Ecuyer-CMRG")

# Windows doesn't support multicore system , which is easier
#  cv_sl = CV.SuperLearner(Y = mytx, X = mymatrix, family = binomial(), V = 10,
#                          parallel = "multicore",
#                          SL.library = my.library)


# Make a snow cluster
cluster = parallel::makeCluster(num_cores-1)
# Load the SuperLearner package on all workers so they can find
# SuperLearner::All(), the default screening function which keeps all variables.
parallel::clusterEvalQ(cluster, library(SuperLearner))

# We need to explictly export our custom learner functions to the workers.
parallel::clusterExport(cluster, xgb.learners$names)
parallel::clusterExport(cluster, rang.learners$names)
parallel::clusterExport(cluster, "SL.gam.2")
parallel::clusterExport(cluster, "SL.gam.3")


####data
# If your data is called  bdata 
# where the outcome is Y and 
# the exposure is tx
# create an outcome and treatment variables as follows

myoutcome <- bdata$Y
mytx  <- bdata$tx

#select all the potential confounders, excluding mediators, called potential_confounders_cols
#define a matrix with them
mymatrix <-bdata[,!colnames(bdata) %in% c(potential_confounders_cols)]

data_obs = as.data.frame(mymatrix)
## this does not contain tx or outcome
data_obs_x = as.data.frame(cbind(data_obs,mytx))
data_obs_xy = as.data.frame(cbind(data_obs,mytx,myoutcome))

Y<-myoutcome
A<-mytx

# set a  seed for SNOW parallelization.
parallel::clusterSetRNGStream(cluster, 11029)

  cv_sl_0 = CV.SuperLearner(Y = myoutcome[mytx==0], X = data_obs[mytx==0,], family = binomial(), V = 10,
                          parallel = cluster,
                          SL.library = my.library)

# Review results.
sum.discrete.SL0<-summary(cv_sl_0)
table.discrete.SL0<-table(simplify2array(cv_sl_0$whichDiscreteSL))
plot(cv_sl_0) + theme_bw()



##now for the exposed 
parallel::clusterSetRNGStream(cluster, 1)

  cv_sl_1 = CV.SuperLearner(Y = myoutcome[mytx==1], X = data_obs[mytx==1,], family = binomial(), V = 10,
                            parallel = cluster,
                            SL.library = my.library)


sum.dscrete.SL1<-summary(cv_sl_1)
table.discrete.SL1<-table(simplify2array(cv_sl_1$whichDiscreteSL))
coef.SL1<-coef(cv_sl_1)
plot(cv_sl_1) + theme_bw()


#run SL to obtain Q0 predictions using the function above
parallel::clusterSetRNGStream(cluster, 1)

SL.object=my.snowSL.pred(mymatrix,mytx,myoutcome,my.library)

#model for unexposed to A
SL.object$m0
#model for exposed to A
SL.object$m1
#ATE from Superlearner
SL.object$ate_SL

Q0<-as.data.frame(SL.object$Q.SL.object)
names(Q0)<-cbind("QAW","Q0W","Q1W")
#save for later use, for example for TMLE
save(Q0, file = "Q0.RData")
#save(SL.object,file="SL.Q0object")

## Now SuperLearner for the propensity score ####
parallel::clusterSetRNGStream(cluster, 1128)

PS <- snowSuperLearner(mytx, data_obs, newX = data_obs, SL.library = my.library, cluster = cluster,
                       family = binomial(),cvControl = list(shuffle = F, V=10),method="method.NNLS")
summary(PS)
PS$coef
g1W <- PS$SL.predict

save(g1W, file = "g0.RData")
#save(PS, file = "PS_SLmodel.RData")


#### Run with external CV to check performace of the SL 
parallel::clusterSetRNGStream(cluster, 19)

  cv_sl_ps = CV.SuperLearner(Y = mytx, X = data_obs, family = binomial(), V = 10,
                            parallel = cluster,
                            SL.library = my.library)



# Review results.
summary(cv_sl_ps)
plot(cv_sl_ps) + theme_bw()
#ggsave("SuperLearner_PSmodel.png")


# Stop the cluster 
parallel::stopCluster(cluster)
