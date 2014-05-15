#-----------------------------------------------------------------------------------------------------#
# Author: Machteld Varewyck
# Creation date : MV 10-06-2013
# R version: 2.14.1
#	Subject: Accompanying R-code to paper
#             "On shrinkage and model extrapolation in statistical methods
#              for assessing clinical center performance"
#              Varewyck, M., Goetghebeur, E., Eriksson, M. and Vansteelandt, S.
#
#
# For mixed-effects method you need additional Bug-files: mixed-classical.bug and mixed-clustered.bug
#     located at your working directory
#-----------------------------------------------------------------------------------------------------#


#Install packages
install.packages(c("rjags","brglm","nnet"))
#Set your working directory
setwd("C:/Users/mmvarewy/Documents/Doctoraat/2012/Code/Package")

#Load packages
library(rjags)    #01 Bayesian in JAGS
library(brglm)		#02/03 Firth correction
library(nnet)     #03 Multinomial regression model

#Define functions
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))
dx.expit <- function(x) expit(x)*(1-expit(x))



### INPUT ###
### ----- ###

# You can assess the center performance based on simulated data or your own data 

## START - Simulated Data ##

set.seed(130513)

n <- 10000    #total nr of patients e.g. 5000
m <- 100      #initial nr of centers (min. 10) e.g. 70
l <- 5        #number of patient-specific covariates in outcome model
l.PS <- l     #number of patient-specific covariates in PS model


#Patient specific characteristics L
L <- matrix(c(rnorm(n,0,10),rbinom(n,1,0.6),rmultinom(n,1,c(0.6,0.3,0.1))), nrow=n)
#Binary center indicators per individual C
C <- t(rmultinom(n, 1, rnorm(m,1/m,1/(4*m))))
#Binary outcome Y: depends on center effect and covariate-specific effect
psi.true <- 5*rbeta(m,0.5,0.5)   #center effect
beta.true <- c(-1.5,-0.08,-0.02,0.5,2)   #covariate-specific effects
Y <- rbinom(n,1,expit(-15+cbind(L,C)%*%c(beta.true,psi.true)))

## END - Simulations ##


## START - Own Data ##

stroke.data <- read.table("C:/Users/mmvarewy/Documents/Doctoraat/2012/Riks-Stroke/cleaned_stroke.txt",
                          header=T, sep = "\t")
ids <- sample(nrow(stroke.data), 10000)
data <- stroke.data[ids,]

#Center per individual
center <- as.factor(data$sjukhuskod)
#Binary center indicators per individual
C <- diag(nlevels(center))[center,]
#Patient specific characteristics
L <- cbind(data$alder, data$sex, data$smoker)
#Binary outcome
Y <- ifelse(data$timetoevent<31 & data$censor==0, 1, 0)

n <- dim(L)[1]
l <- dim(L)[2]
m <- nlevels(center)

## END - Own Data ##



# Clinical tolerance level, lambda
lambda <- 0.2
# Statistical tolerance level, k
k <- 0.75



### INSTRUCTIONS ###
### ------------ ###

#First, run the function 'model.fitting()' below, 
# next, choose the statistical method to analyze the data
# then you get numerical output (in Console), 
# output data in file "Output-methodXX.txt" (in working directory)
# and the following figures (in working directory):
#   - Plot E{Y(c)} vs center size + 50% CIs "FullPopulationRisk_CI-methodXX.png"
#   - Plot observed vs expected outcome per center "Observed_Expected-methodXX.png"


#Statistical method
#  method: 1 = classical Bayesian normal mixed effects method
#          2 = clustered (N=3) Bayesian normal mixed effects method  !Attention: takes very long time to fit!
#          3 = fixed effects method
#          4 = doubly robust propensity scores method
# -- Additionally for methods 3 and 4 --
#  Firth: F = Do not apply the Firth correction
#         T = Apply the Firth correction on outcome model (for method 3),
#             Apply the Firth correction on outcome model and logistic PS model (for method 4)
# -- Additionally for method 4 --
#  multinomial: F = Fit a logistic model for the propensity scores 
#               T = Fit a multinomial model for the propensity scores


results.method1 <- model.fitting(method=1)
results.method2 <- model.fitting(method=2)
results.method3 <- model.fitting(method=3, Firth=T)
results.method4 <- model.fitting(method=4, Firth=T, multinomial=T)



## START - Function 'model.fitting' ##

model.fitting <- function(method, Firth, multinomial){

  
### CALCULATIONS ###
### ------------ ###
  
## General calculations ##

#Model matrix for PS
X.PS <- cbind(1,L)
#Model matrix for outcome, center 1 as reference
X <- cbind(1,L,C[,-1])
##alternative: X <- model.matrix(~ L + C[,-1])

center <- apply(t(C)*c(1:m), 2, sum)
center.size <- table(center)

C.pop <- diag(m)[rep(1:m,each=n),]
X.pop <- cbind(1,L[rep(1:n,times=m),],C.pop[,-1])


time1 <- Sys.time()


if(method == 1){
  
  # Method.1: Classical Bayesian normal mixed-effects 
  # -------------------------------------------------
  
  input <- list("n"=n, "m"=m, "l"=l, "y"=Y, "center"=factor(center), "L"=L, "lambda"=lambda)
  #  n = the number of patients
  #  m = the number of centers
  #  l = the number of patient-specific covariates
  #  y = observed binary outcome for each patient
  #  center = observed center number (1..m) for each patient
  #  L = observed patient characteristics
  #  lambda = clinical tolerance level
  
  #parameters that need to be saved
  parameters <- c("center.effect", "beta", "mean.pop.out", "mean.out.lambdaM", "mean.out.lambdaP")
  
  #JAGS model fitting
  model.1 <- jags.model("mixed-classical.bug", data = input, n.chains = 1, n.adapt = 3000, quiet=F)
  mixed.sim <- jags.samples(model.1, parameters, n.iter=5000)
  
  # Estimate outcomes
  mean.pop.out <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.pop.out)), 2, mean)
  
  # Estimate variance
  #saves Var[(mean.pop.out)]
  var.pop.out <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.pop.out)), 2, var)
  
  # Calculate 95% credible interval limits for mean.pop.out
  lower.pop.CI <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.pop.out)), 2, quantile, 0.025)
  upper.pop.CI <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.pop.out)), 2, quantile, 0.975)
  
  #P(E[Y(c)] > (1+lambda)*E[Y])
  probM <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.out.lambdaM))<0, 2, mean)
  #P(E[Y(c)] < (1-lambda)*E[Y])
  probP <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.out.lambdaP))>0, 2, mean)
  

  
}else if(method == 2){
  
  # Method.2: Clustered Bayesian normal mixed-effects 
  # -------------------------------------------------
  
  input <- list("n"=n, "m"=m, "l"=l, "y"=Y, "center"=factor(center), "L"=L, "lambda"=lambda, "N"=3)
  
  #  n = the number of patients
  #  m = the number of centers
  #  l = the number of patient-specific covariates
  #  y = observed binary outcome for each patient
  #  center = observed center number (1..m) for each patient
  #  L = observed patient characteristics
  #  lambda = clinical tolerance level
  #  N = the number of center clusters
  
  #parameters that need to be saved
  parameters <- c("center.effect", "beta", "mean.pop.out", "mean.out.lambdaM", "mean.out.lambdaP")
  
  #JAGS model fitting
  model.2 <- jags.model("mixed-clustered.bug", data = input, n.chains = 1, n.adapt = 3000, quiet=F)
  mixed.sim <- jags.samples(model.2, parameters, n.iter=5000)
  
  # Estimate outcomes
  mean.pop.out <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.pop.out)), 2, mean)
  
  # Estimate variance
  #saves Var[(mean.pop.out)]
  var.pop.out <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.pop.out)), 2, var)
  
  # Calculate 95% credible interval limits for mean.pop.out
  lower.pop.CI <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.pop.out)), 2, quantile, 0.025)
  upper.pop.CI <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.pop.out)), 2, quantile, 0.975)
  
  #P(E[Y(c)] > (1+lambda)*E[Y])
  probM <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.out.lambdaM))<0, 2, mean)
  #P(E[Y(c)] < (1-lambda)*E[Y])
  probP <- apply(as.matrix(as.mcmc.list(mixed.sim$mean.out.lambdaP))>0, 2, mean)
  
  
  
}else if(method == 3){
  
  
  # Method.3: Fixed-effects without/with Firth correction
  # -----------------------------------------------------
  
  # Model fitting  
  if(Firth==F){
    model.3 <- glm(Y ~ X[,-1], family = binomial)
    b <- as.matrix(model.3$coef)
  }else if(Firth==T){
    model.3 <- brglm(Y ~ X[,-1], family = binomial(logit), method = "brglm.fit", pl=TRUE)
    b <- as.matrix(model.3$coef)
  }
  
  Lb <- cbind(1,L)%*%b[1:(l+1)]    # n x 1
  psi <- c(0,b[-(1:(l+1))])    # length = m
  
  
  # Estimate outcomes
  p.expit <- expit(X%*%b)  ##U(\theta,X)
  p.pop.expit <- matrix(expit(X.pop%*%b),nrow=n)  ##U_c(\theta,X)
  
  mean.pop.out <- apply(p.pop.expit, 2, mean)   #E[Y(c)]
  #check: mean(p.expit); mean(fitted.values(model.1)); mean(Y)
  
  # Estimate variance
  R <- X*as.vector(Y - p.expit)
  deriv <- dx.expit(as.vector(X%*%b))*X
  A <- 1/n*t(X)%*%deriv
  Ainv <- solve(A)   #square matrix with dim l+m
  #check: summary(round(vcov(model.0)) == round(Ainv/n))
  #symmetric matrix: summary(t(round(Ainv,5)) == round(Ainv,5))
  
  help <- dx.expit(as.vector(X.pop%*%b))*X.pop
  deriv.pop <- matrix(apply(matrix(help, nrow=n), 2, mean), nrow=(l+m), byrow=T)
  #returns l+m x m matrix, with each column partial derivatives for one center
  #as.vector plaatst de kolommen (nt rijen) van matrix achter elkaar -> OK
  
  #saves Var[logit(mean.pop.out)]
  var.pop.out <- 1/(n*(mean.pop.out*(1-mean.pop.out))^2)*apply(p.pop.expit + R%*%Ainv%*%deriv.pop,2,var)
  
  #saves Var[logit(Y_center) - logit(1-lambda)*(Y_overall)]
  var.lambdaM.pop <- 1/n*apply(t(1/(mean.pop.out*(1-mean.pop.out))*t(p.pop.expit + R%*%Ainv%*%deriv.pop)
                                 - matrix(rep(Y/(mean(Y)*(1-(1-lambda)*mean(Y))), m), nrow=m, byrow=T)),2,var)
  #saves Var[logit(Y_center) - logit(1+lambda)*(Y_overall)]
  var.lambdaP.pop <- 1/n*apply(t(1/(mean.pop.out*(1-mean.pop.out))*t(p.pop.expit + R%*%Ainv%*%deriv.pop)
                                 - matrix(rep(Y/(mean(Y)*(1-(1+lambda)*mean(Y))), m), nrow=m, byrow=T)),2,var)
  

  
}else if(method == 4){
  
  
  # Method 4: Double robust PS method
  # ---------------------------------
  
  # Model fitting for PS
  
  ## Multinomial regression for PS ##
  if(multinomial==T){
  model.PS <- multinom(center ~ X.PS[,-1], trace=F)
  PS.all <- model.PS$fitted.values
  }
  
  ## Logistic regression for PS ##
  if(multinomial==F){
  PS.all.log <- matrix(NA, nrow=n, ncol=m)
  b.PS <- matrix(NA, nrow=l.PS+1, ncol=m)
  
  if(Firth==F){
    for(i in 1:m){  
      model.PS <- glm((center==i)*1 ~ X.PS - 1, family=binomial(link="logit"))
      PS.all.log[,i] <- predict(model.PS, newdata=as.data.frame(X.PS), type="response")
      b.PS[,i] <- as.vector(model.PS$coef)
    }
  }else if(Firth==T){
    for(i in 1:m){  
      model.PS <- brglm((center==i)*1 ~ X.PS - 1, family=binomial(link="logit"), method = "brglm.fit", pl=TRUE)
      PS.all.log[,i] <- predict(model.PS, newdata=as.data.frame(X.PS), type="response")
      b.PS[,i] <- as.vector(model.PS$coef)
    }
  }
  PS.all <- (PS.all.log/apply(PS.all.log, 1, sum))
  }
  
  #stabilized weights
  w.all <- t(as.vector(center.size)/(n*t(PS.all)))*C
  w <- apply(w.all, 1, sum)    #vector of individual weights (1/PS)
  

  # Model fitting for outcome
  if(Firth==F){
    model.4 <- glm(Y ~ X-1, family = binomial("logit"), weights=w)
  }else if(Firth==T){
    model.4 <- brglm(Y ~ X-1, family = binomial("logit"), weights=w)
  }
  b <- as.matrix(model.4$coef)
  
  Lb <- cbind(1,L)%*%b[1:(l+1)]    # n x 1
  psi <- c(0,b[-(1:(l+1))])    # length = m
  
  
  # Estimate outcomes
  p.expit <- expit(X%*%b)  ##U(\theta,X)
  p.pop.expit <- matrix(expit(X.pop%*%b),nrow=n)  ##U_c(\theta,X)
  
  mean.pop.out <- apply(p.pop.expit, 2, mean)   #E[Y(c)]
  #check: mean(p.expit); mean(fitted.values(model.3)); mean(Y)
  
  
  # Estimate variance

  # Part1: From the estimating equations for PS model #
  
  L1 <- cbind(1,L)
  m1 <- diag(m-1)
  
  ##MULTINOMIAL##
  if(multinomial==T){
  resid.PS <- C[,2:m] - PS.all[,2:m]
  n1 <- L1[,rep(1:(l+1),times=m-1)]
  t4 <- function(ind){
    t1 <- as.vector(matrix(L1[ind,],ncol=1) %*% matrix(PS.all[ind,2:m], nrow=1))
    t2 <- m1 - PS.all[ind,2:m]
    t3 <- t(L1[ind,]*t2[rep(1:(m-1), each=l+1),])
    t1 * t3[rep(1:(m-1),each=l+1),]
  }
  A1 <- matrix(0, nrow=(m-1)*(l+1), ncol=(m-1)*(l+1))
  for(j in 1:n){A1 <- A1 + t4(j)}   #600 ind/minute (10 000 is 1h30)
  A1inv <- solve(A1/n) #singular? library(MASS); A1inv <- ginv(A1)
  #check: summary(round(vcov(model.PS)) == round(A1/n))
  R1 <- n1*resid.PS[,rep(1:(m-1),each=l+1)]
  }
  
  ##LOGISTIC##
  if(multinomial==F){
  resid.PS <- (C-PS.all.log)
  A1 <- matrix(0, nrow=m*(l.PS+1), ncol=m*(l.PS+1))
  for(i in 1:m){
    deriv <- dx.expit(as.vector(L1%*%b.PS[,i]))*L1
    A.log <- 1/n*t(L1)%*%deriv
    A1[(((i-1)*(l.PS+1)+1):(i*(l.PS+1))),(((i-1)*(l.PS+1)+1):(i*(l.PS+1)))] <- A.log
  }
  A1inv <- solve(A1)
  #check for last center only: summary(round(vcov(model.PS)) == round(A.log/n))
  R1 <- L1[,rep(1:(l+1),times=m)]*resid.PS[,rep(1:m,each=(l+1))]
  #R1.sum <- t(L1)%*%residx
  }
  remove(A1); gc()
  rho.PS <- 1/sqrt(n)*apply(R1,2,sum)%*%A1inv  #vector: center2, intercept; center2, Li1; ...; center3, intercept; ...
  
  # Part2: From the estimating equations for outcome model #
  
  resid.out <- as.vector(Y - p.expit)
  R2 <- w*resid.out*X
  
  ##MULTINOMIAL##
  if(multinomial==T){
  t1 <- w*resid.PS
  t2 <- matrix(t1[,rep(1:(m-1),each=l+1)], nrow=n)
  t3 <- ((n1*t2)*resid.out)%*%t(rho.PS)
  R3 <- R2 - 1/sqrt(n)*X*as.vector(t3)
  }
  
  ##LOGISTIC##
  if(multinomial==F){
  t1 <- c()
  for(i in 1:n){
    t0 <- matrix(L1[i,])%*%t(matrix((as.vector(center.size[center[i]])/n-w.all[i,])
                                    *dx.expit(as.vector(L1[i,]%*%b.PS))/PS.all.log[i,center[i]]))
    t1[i] <- matrix(t0,nrow=1)%*%t(rho.PS)*resid.out[i]
  }
  R3 <- R2 + 1/sqrt(n)*X*t1
  #R3.sum <- apply(R3, 2, sum)
  }
  
  deriv <- w*dx.expit(as.vector(X%*%b))*X
  A2 <- 1/n*t(X)%*%deriv  
  A2inv <- solve(A2)
  #check impossible (geen MLE)
  theta.out <- 1/sqrt(n)*apply(R3,2,sum)%*%A2inv
  
  
  # Part 3: Calculate sample variance
  
  ##MULTINOMIAL##
  if(multinomial==T){
  t1 <- -apply(as.vector(center.size/n)*t(C),2,sum)/PS.all[,2:m]*resid.out*resid.PS
  t2 <- n1*t1[,rep(1:(m-1),each=l+1)]
  E.dV.rho <- matrix(NA, nrow=m, ncol=ncol(t2))
  for(i in 1:m){
    E.dV.rho[i,] <- apply(t2*C[,i], 2, mean)
  }
  #t3 <- t2[,rep(1:ncol(t2),times=m)]*C[,rep(1:m,each=ncol(t2))]
  #E.dV.rho <- matrix(apply(t3, 2, mean),nrow=m,byrow=T)
  
  t3 <- dx.expit(as.vector(X.pop%*%b))*X.pop  
  E.dV.theta <- matrix(apply(matrix(t3*(1-w*as.vector(t(C))), nrow=n),2,mean),nrow=m)
  remove(t1,t2,t3); gc()
  }
  

  ##LOGISTIC##
  if(multinomial==F){
  t1 <- dx.expit(X.PS%*%b.PS)
  E.dV.rho <- matrix(NA, nrow=(l.PS+1)*m, ncol=m)
  E.dV.theta <- matrix(NA, nrow=(l+m),ncol=m)
  for(i in 1:m){
    E.dV.rho[,i] <- 1/n*as.vector(t(X.PS[which(center==i),])%*%((as.vector(center.size[center])/n-w.all)
                                *t1*resid.out/apply(PS.all.log*C,1,sum))[center==i,])
    E.dV.theta[,i] <- 1/n*t(X.pop[((i-1)*n+1):(i*n),])%*%((1-w.all[,i])*dx.expit(as.vector(Lb + psi[i])))
  }
  E.dV.rho <- t(E.dV.rho)
  E.dV.theta <- t(E.dV.theta)
  }
  
  remove(X.pop)
  #saves Var[logit(mean.pop.out)]
  var.pop.out <- 1/(n*(mean.pop.out*(1-mean.pop.out))^2)*apply(p.pop.expit + w*C*resid.out + (R1%*%A1inv)%*%t(E.dV.rho) + (R3%*%A2inv)%*%t(E.dV.theta), 2, var)
  #saves Var[logit(Y_center) - logit(1-lambda)*(Y_overall)]
  denom1 <- matrix(rep(mean.pop.out*(1-mean.pop.out), n), nrow=n, byrow=T)
  var.lambdaM.pop <- 1/n*apply(1/denom1*(p.pop.expit + w*C*resid.out + (R1%*%A1inv)%*%t(E.dV.rho) + (R3%*%A2inv)%*%t(E.dV.theta))
                               - matrix(rep(Y/(mean(Y)*(1-(1-lambda)*mean(Y))), m), nrow=n),2,var)
  #saves Var[logit(Y_center) - logit(1+lambda)*(Y_overall)]
  var.lambdaP.pop <- 1/n*apply(1/denom1*(p.pop.expit + w*C*resid.out + (R1%*%A1inv)%*%t(E.dV.rho) + (R3%*%A2inv)%*%t(E.dV.theta))
                               - matrix(rep(Y/(mean(Y)*(1-(1+lambda)*mean(Y))), m), nrow=n),2,var)
  
  
  
}



### OUTPUT ###
### ------ ###

cat("Time needed to perform model fitting: \n")
print(Sys.time()-time1)

cat("Mean incidence rate \n")
print(round(mean(Y), 3))
cat("Center size \n")
print(round(center.size))


#Check VAR
if(method %in% c(3,4)){
  
  # Calculate 95% confidence interval limits for mean.pop.out
  lower.pop.CI <- expit(logit(mean.pop.out) - qnorm(0.975)*sqrt(var.pop.out))
  upper.pop.CI <- expit(logit(mean.pop.out) + qnorm(0.975)*sqrt(var.pop.out))
  
  #P(E[Y(c)] > (1+lambda)*E[Y])
  probP <- 1-(pnorm(-logit(mean.pop.out)+logit(mean(Y)*(1+lambda)), 0, sqrt(var.lambdaP.pop)))
  #P(E[Y(c)] < (1-lambda)*E[Y])
  probM <- pnorm(-logit(mean.pop.out)+logit(mean(Y)*(1-lambda)), 0, sqrt(var.lambdaM.pop))
    
}

ratio <- tapply(Y, center, sum)/as.vector(center.size)
cat("Center-specific risk E(Y|C=c) \n")
print(round(ratio,3))
cat("Potential full population risk E{Y(c)} \n")
print(round(mean.pop.out, 3))
cat("Estimated 95% lower bound for the potential full population risk \n")
print(round(lower.pop.CI, 3))
cat("Estimated 95% upper bound for the potential full population risk \n")
print(round(upper.pop.CI, 3))


#Each center gets an indication (L=low mortality risk, A=accepted, H=high mortality risk)
result <- rep('A', m)
result[which(probM>k)] <- 'L'
result[which(probP>k)] <- 'H'
result[which(is.na(probP)==T | is.na(probM)==T)] <- NA

cat("Center classification labels \n")
print(result)

#Save output in file (in working directory)
write.table(data.frame(center.size=as.vector(center.size), center.factor=levels(as.factor(center)), result=result,
         lower.pop.CI=lower.pop.CI, upper.pop.CI=upper.pop.CI, mean.pop.out=mean.pop.out,
         ratio=ratio, probM=probM, probP=probP),                 
    file=paste("Output-method",method,".txt",sep=""))

#This data-file can be read as e.g.
#data <- read.table("Output-method1.txt", header=T)



par()$mar #5.1 4.1 4.1 2.1


##Plot E{Y(c)} vs center size + 50% CIs

png(filename=paste("FullPopulationRisk_CI-method",method,".png",sep=""),
    width = 1000, height = 800)
par(mar=c(5.1,5.1,2.1,2.1))

center.low <- center.size[result=='L']
center.acc <- center.size[result=='A']
center.high <- center.size[result=='H']

plot(mean.pop.out[result=="L"], center.low, pch=18, col="blue", log="y",
     ylim=c(min(center.size),max(center.size)),xlim=c(min(lower.pop.CI),max(upper.pop.CI)),
     cex=1.5,cex.axis=1.5,cex.lab=1.5, ylab="Center size",xlab="Ê{Y(c)} and 50% CI")
points(mean.pop.out[result=="A"], center.acc,pch=19,col="black")  
points(mean.pop.out[result=="H"], center.high,pch=18,cex=1.5,col="red")
#50% confidence limits
lower.pop.CI50 <- expit(logit(mean.pop.out) - qnorm(0.75)*sqrt(var.pop.out))
upper.pop.CI50 <- expit(logit(mean.pop.out) + qnorm(0.75)*sqrt(var.pop.out))
for(i in 1:m){
  x.index <- as.vector(center.size)[i]
  lines(c(lower.pop.CI50[i],upper.pop.CI50[i]),c(x.index,x.index))
}
abline(v=(1-lambda)*mean(Y), col="darkgray", lty=3,lwd=2)
abline(v=(1+lambda)*mean(Y), col="darkgray", lty=3,lwd=2)

dev.off()


##Plot observed vs expected outcome per center

png(filename=paste("Observed_Expected-method",method,".png",sep=""),
    width = 1000, height = 800)
par(mar=c(5.1,5.1,2.1,2.1))

ratio.low <- as.vector(ratio)[result=="L"]
ratio.acc <- as.vector(ratio)[result=="A"]
ratio.high <- as.vector(ratio)[result=="H"]
plot(mean.pop.out[result=="L"], ratio.low, pch=18,col='blue',
     xlim=c(min(mean.pop.out),max(mean.pop.out)), ylim=c(min(ratio),max(ratio)),
     cex=1.5,cex.axis=1.5,cex.lab=1.5, ylab="Ê(Y|C=c)",xlab="Ê{Y(c)}")
points(mean.pop.out[result=="A"], ratio.acc, pch=19,col="black",cex.axis=1.5)  
points(mean.pop.out[result=="H"], ratio.high, pch=18,col="red",cex=1.5,cex.axis=1.5)
lines(loess.smooth(mean.pop.out, ratio),lwd=2)
points(mean(Y),mean(Y),pch='M',col="darkgray",cex=3)
abline(a=0,b=1,lwd=2,lty=2)

dev.off()

}

## END - Function 'model.fitting' ##