##############################################
## Variable selection for birth data
##  Double lasso, and lasso-CTMLE
##  September 2018
##  Karla DiazOrdaz
#############################################
library(SuperLearner)
library(foreign)
library(glmnet)
library(caret)
library(tidyverse)
library(dplyr)
library(boot)
library(survey)

# If your data is called  bdata 
# where the outcome is Y and 
# the exposure is tx
# create an outcome and treatment variables as follows

myoutcome <- bdata$Y
mytx  <- bdata$tx

#select all the potential confounders, excluding mediators, called potential_confounders_cols
#define a matrix with them
mymatrix <-bdata[,!colnames(bdata) %in% c(potential_confounders_cols)]


#### double lasso ####
## we use glmnet which needs a matrix of predictors x
## this is mymatrix
##
## and a  vector of outcome variable myoutcome
set.seed (1)
train=sample (1: nrow(mymatrix), nrow(mymatrix)/2)
test=(- train )
y.test=myoutcome[test]


###models need to be linear !!
# STEP 1: select variables that predict outcomes using lasso

#Belloni 2014 does not choose the penalty by CV, but we do, following
#Dukes 2018 who found it to behave better 

cv.out =cv.glmnet(mymatrix[train ,],myoutcome[train],family="gaussian",alpha =1)
plot(cv.out)
#select best value for the penalty lambda using CV
bestlam =cv.out$lambda.min  
lasso.out = glmnet(mymatrix[train ,],myoutcome[train],family="gaussian",alpha =1, lambda=bestlam)


lasso.pred.out=predict(lasso.out ,s=bestlam ,newx=mymatrix[test ,])

out=glmnet(mymatrix,myoutcome,alpha =1, family="gaussian",lambda =bestlam)
lasso.coef.out=predict(out ,type ="coefficients",s=bestlam )
selected.vars.out<-t(data.frame(name = lasso.coef.out@Dimnames[[1]][lasso.coef.out@i + 1], coefficient = lasso.coef.out@x))[1,]

# STEP 2: select variables that predict treatment using lasso
tx.test=mytx[test]

set.seed (101)
cv.tx =cv.glmnet(mymatrix[train ,],mytx[train],family="gaussian",alpha =1)
plot(cv.tx)
#select best value for the penalty lambda using CV
bestlam.tx =cv.tx$lambda.min
lasso.tx =glmnet (mymatrix[train ,],mytx[train], family="gaussian", alpha =1, lambda =bestlam.tx)

lasso.pred.tx=predict(lasso.tx ,s=bestlam.tx ,newx=mymatrix[test ,])
mean((lasso.pred.tx-tx.test)^2) #test MSE

out.tx=glmnet(mymatrix,mytx,family="binomial",alpha =1, lambda =bestlam.tx)
lasso.coef.tx=predict(out.tx ,type ="coefficients",s=bestlam.tx )

selected.vars.tx<-t(data.frame(name = lasso.coef.tx@Dimnames[[1]][lasso.coef.tx@i + 1], coefficient = lasso.coef.tx@x))[1,]
length(t(data.frame(name = lasso.coef.tx@Dimnames[[1]][lasso.coef.tx@i + 1], coefficient = lasso.coef.tx@x))[1,])

# STEP 3: linear regression with both sets of variables, no mediator
use = union(selected.vars.out,selected.vars.tx)


X <-  select(as.data.frame(mymatrix), one_of(use[-1]))
#bind treatment and outcome with selected covarites
XAY<- as.data.frame(cbind(X,mytx,myoutcome))
regresors<- paste(names(X),collapse=" + ")

#outcome regression
y_form  = paste("myoutcome","~","mytx","+", paste(names(X),collapse=" + "), collapse = "") 
#propensity score model
a_form  = paste("mytx", "~", paste(names(X),collapse=" + "), collapse = "") 

fit.out.reg <- glm(as.formula(y_form), data=XAY, family=gaussian())
summary(fit.out.reg) # show results

#predict potential outcomes
Y1 = predict(fit.out.reg, newdata = data.frame(mytx = 1, X), type="response")
Y0 = predict(fit.out.reg, newdata = data.frame(mytx = 0, X), type="response")

ACE=mean(Y1)-mean(Y0) #we would have to bootstrap for SEs

###bootstrapping only the 3rd step for SEs
ACE_doublelasso_boot <- function(y_form, data, indices){
  data_o = data[indices,]
  fit.out.reg <- glm(as.formula(y_form), data=data_o, family=gaussian())
  Y1 = predict(fit.out.reg, newdata = data.frame(mytx = 1, data_o), type="response")
  Y0 = predict(fit.out.reg, newdata = data.frame(mytx = 0, data_o), type="response")
  ACE=mean(Y1)-mean(Y0)
 return(ACE=ACE) 
}

####

ACE.boot <- boot(y_form=y_form, data=XAY,ACE_doublelasso_boot, R = 1999)
boot.ci.ace<-boot.ci(ACE.boot, conf = 0.95, type = "perc")
ACE.LCI<-boot.ci.ace$percent[4]
ACE.UCI<-boot.ci.ace$percent[5]

###


## so we do lasso for A and Y to select variables
## take the union and 
## re-run models for Y and A using this control set BUT 
## not lasso, just regular regression
#### double lasso with INTERACTIONS####

#we can extend X by adding interactions and other non-linear terms and try double-lasso
##Extending with interactions and non-linear terms ##



mydata <-bdata[,!colnames(bdata) %in% c(potential_confounders_cols)]

mydata_inter <- t(apply(mydata, 1, combn, 2, prod))
colnames(mydata_inter) <- paste("Inter.V", combn(dim(mydata)[2], 2, paste, collapse="V"), sep="")

mydata_main_inter<-as.data.frame(cbind(mydata,mydata_inter))

mymatrix_inter <-as.matrix(mydata_main_inter)
#higher order terms and
#interactions of the control variables defined above.

set.seed (1)
train=sample (1: nrow(mymatrix_inter), nrow(mymatrix_inter)/2)
test=(- train )
y.test=myoutcome[test]

# STEP 1: select variables that predict outcomes using lasso
cv.out =cv.glmnet(mymatrix_inter[train ,],myoutcome[train],family="gaussian",alpha =1)
#select best value for the penalty lambda using CV
bestlam =cv.out$lambda.min
lasso.out=glmnet(mymatrix_inter,myoutcome,alpha =1, family="gaussian",lambda =bestlam)

lasso.pred.out=predict(lasso.out ,s=bestlam ,newx=mymatrix_inter[test ,])
mean(( lasso.pred.out -y.test)^2) #test MSE

#extracting vars chosen for the outcome model

selected.vars.out<-t(data.frame(name = lasso.coef.out@Dimnames[[1]][lasso.coef.out@i + 1], coefficient = lasso.coef.out@x))[1,]

# STEP 2: select variables that predict treatment using lasso
tx.test=mytx[test]


set.seed (101)
cv.tx =cv.glmnet(mymatrix_inter[train ,],mytx[train],family="gaussian",alpha =1)
#select best value for the penalty lambda using CV
bestlam.tx =cv.tx$lambda.min
lasso.tx =glmnet (mymatrix_inter[train ,],mytx[train], family="gaussian", alpha =1, lambda =bestlam.tx)
                    
lasso.pred.tx=predict(lasso.tx ,s=bestlam.tx ,newx=mymatrix_inter[test ,])
mean((lasso.pred.tx-tx.test)^2) #test MSE

out.tx=glmnet(mymatrix_inter,mytx,family="gaussian",alpha =1, lambda =bestlam.tx)
lasso.coef.tx=predict(out.tx ,type ="coefficients",s=bestlam.tx )

selected.vars.tx<-t(data.frame(name = lasso.coef.tx@Dimnames[[1]][lasso.coef.tx@i + 1], coefficient = lasso.coef.tx@x))[1,]

##which variables are not selected

# STEP 3: linear regression with both sets of variables, no mediator
use.inter = union(selected.vars.out,selected.vars.tx)

#bind treatment with covarites
X.inter <-  select(as.data.frame(mymatrix_inter), one_of(use.inter[-1]))
#names(X)
XAY.inter<- as.data.frame(cbind(X.inter,mytx,myoutcome))
#regresors<- paste(names(X),collapse=" + ")
y_form_inter  = paste("myoutcome", "~","mytx","+", paste(names(X.inter),collapse=" + "), collapse = "") 
a_form_inter  = paste("mytx", "~", paste(names(X.inter),collapse=" + "), collapse = "") 

fit.out.inter <- glm(as.formula(y_form_inter), data=XAY.inter, family=gaussian())
summary(fit.out.inter) # show results
ATE.outreg.inter<-fit.out.inter$coefficients["mytx"]
SE.outreg.inter<-sqrt(vcov(fit.out.inter)[2,2])



###bootstrapping only the 3rd step for SEs
ACE_doublelasso_boot <- function(formula, data, indices){
  data_o = data[indices,]
  fit.out.reg <- glm(as.formula(formula), data=data_o, family=gaussian())
  Y1 = predict(fit.out.reg, newdata = data.frame(mytx = 1, data_o), type="response")
  Y0 = predict(fit.out.reg, newdata = data.frame(mytx = 0, data_o), type="response")
  ACE=mean(Y1)-mean(Y0)
  return(ACE=ACE) 
}

####

set.seed(101)
ACE.boot.inter <- boot(formula=y_form_inter, data=XAY.inter,ACE_doublelasso_boot, R = 1999)
boot.ci.ace.inter<-boot.ci(ACE.boot.inter, conf = 0.95, type = "perc")

ACE.LCI.inter<-boot.ci.ace.inter$percent[4]
ACE.UCI.inter<-boot.ci.ace.inter$percent[5]

