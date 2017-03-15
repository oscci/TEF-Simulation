################################################################################
#
# Simulation of the Mock TEF results - Benchmarking
#
# by Paul Thompson
#
################################################################################

# latest version 12-07-2016

#install.packages(c("extraDistr","copula","plyr"))

#timing the run time for code!
ptm <- proc.time()

require(extraDistr)
require(copula)
require(plyr)
#library(foreach)
#library(doParallel)

#registerDoParallel(cores=4)
#cl <- makeCluster(4)
#registerDoParallel(cl)

#set seed to that simulated values are reproducible 
set.seed(1981)

#set number of UK students
n <- 300000 
nuniv=120
gamma_val=0.5
# 7 variables for benchmarking:

# years (Categorical)
# disabled (Binomial)
# ethnicity (Binomial)
# level (Binomial)
# age (Binomial)
# domicile (Categorical)
# NSS (Gaussian)

#Simulate a population of students at the 120 institutions, so that the demographic is accurate. The correlations between variables can be altered to change the demographic. If specific types of institutions are required, then this must be altered  manually. For instance, an institution that has a higher than average number of international students, or a postgraduate only institution.


## Build a normal copula for the underlying distibution:
# `param` is the lower triangle of a correlation matrix
underlying <- normalCopula(dim = 7, 
                           param = c(-0.4,0,-0.2,-0.5,-0.2,0.2,   # years
                                     0.01,0.05,0.1,0.1,0.1,       # disabled
                                     0.4,0.6,0.3,0.2,             # ethnicity
                                     0.3,0.1,0.3,                 # level
                                     0.5,0.4,                     # age
                                     0.6),                        #domicile
                           dispstr = "un") # The Copula

## build the marginal distibutions from the underlying copula
marginals <- mvdc(copula = underlying, 
                  margins = c("norm","binom", "binom", "binom", "binom", "norm","beta"), 
                  paramMargins = list(list(mean = 0, sd = 1), # parameters for the marginal distributions...
                                      list(size = 1, prob = 0.9),
                                      list(size = 1, prob = 0.2),
                                      list(size = 1, prob = 0.05),
                                      list(size = 1, prob = 0.15),
                                      list(mean = 0, sd = 1),
                                      list(shape1 = 10, shape2 = 3))) #edit these parameters for NSS distribution.

## generate the correlated random variables from the marginal distributions
zdat <- as.data.frame(rMvdc(n, marginals)) 

## rename the variables
names(zdat) <- c( "years", "disabled", "ethnicity", "level", "age", "domicile","NSS")


############# CATEGORIES procedure #########################################

###cite: http://finzi.psych.upenn.edu/library/extraDistr/html/Categorical.html ## for categorical quantile function (qcat).

##  years variable ##########


u<-pnorm(zdat[,1])

years = qcat(u,prob=c(0.4,0.3,0.3))

round(cor(zdat),2)

cor(years,zdat[,2:7],method="spearman")

zdat$years<-years

## Domicile variable ##########

u<-pnorm(zdat[,6])

domicile = qcat(u,prob=c(0.85,0.1,0.05))

round(cor(zdat),2)

cor(domicile,zdat[,c(1:5,7)],method="spearman")

zdat$domicile<-domicile

################################################################################


# convert other variables to factors and clean up
zdat$years <- ordered(zdat$years)
zdat$domicile <- factor(zdat$domicile)
levels(zdat$domicile) <- c("UK", "Other EU","Non EU")
zdat$disabled <- factor(zdat$disabled)
levels(zdat$disabled) <- c("Yes", "No")
zdat$ethnicity <- factor(zdat$ethnicity)
levels(zdat$ethnicity) <- c("White", "BME")
zdat$age <- factor(zdat$age)
levels(zdat$age) <- c("Young", "Mature")
zdat$level <- factor(zdat$level)
levels(zdat$level) <- c("first degree", "Other UG")
zdat <- zdat[complete.cases(zdat),]

################################################################################

#Add university indicator

zdat$university<-rcat(n,rdirichlet(1, sort(rpois(120,4)+1,decreasing =T)))
zdat$NSS<-ifelse(zdat$NSS>0.75,1,0)

################################################################################


#Check that the simulated data has realistic and representative frequencies within each population category.

countt<-count_(zdat,c("years", "disabled", "ethnicity"))
countt

########################################################################################################

#### BENCHMARKING: Following Draper & Gittoes (2004) and BIS: TEF technical consultation for year two 

########################################################################################################

#Basic TEF grids


#Function to calculate PCF cell means
mean_PT <- function(data,para){
  ifelse(dim(filter(data,data$years==para$years, data$disabled==para$disabled, data$ethnicity==para$ethnicity))[1]==0,NA,mean(filter(data,data$years==para$years, data$disabled==para$disabled, data$ethnicity==para$ethnicity)$NSS,na.rm=T))
}

library(dplyr)
#Function to calculate PCF cell counts
counter_PT <- function(data, para){
  ifelse(dim(filter(data,data$years==para$years, data$disabled==para$disabled, data$ethnicity==para$ethnicity))[1]==0,0,dim(filter(data,data$years==para$years, data$disabled==para$disabled, data$ethnicity==para$ethnicity))[1])
}


#create a grid of all combinations of parameters
# parameters <- expand.grid(years=c("1","2","3"),disabled=c("Yes","No"),ethnicity=c("White","BME"))
# parameters$years<-ordered(parameters$years)

parameters<-countt[,1:3]
dim.para<-dim(parameters)

zdat_r<-zdat[,-c(4:6)]

Count_mat<-PCF_mat<-matrix(NA,nrow=120,ncol=dim.para[1])

#Calculate of PCF means and counts for each institution
# ptm<-proc.time()
for(i in 1:120){
  for(j in 1:dim.para[1]){
    PCF_mat[i,j] <- mean_PT(data=zdat_r[zdat_r$university==i,],para=parameters[j,])
    Count_mat[i,j] <- counter_PT(data=zdat_r[zdat_r$university==i,],para = parameters[j,])
  }
}
# ptm<-proc.time()

#calculate the column means and Column sums

col_sums<-apply(Count_mat,2,sum,na.rm=T)
row_sums<-apply(Count_mat,1,sum,na.rm=T)


weight_col_means<-sapply(1:dim.para[1],function(j){(1/col_sums[j])*sum(Count_mat[,j]*PCF_mat[,j],na.rm=T)})


########### Observed and Expected NSS rates (new)

O_hat<-E_hat<-D_hat<-SE_D<-z<-vector(mode="numeric",length=120)


lambda <- array(NA,c(120,120,dim.para[1]))

for(i in 1:120){
  for(k in 1:120){
    for(j in 1:dim.para[1]){
      lambda[i,k,j]<-ifelse(i==k,(Count_mat[i,j]/row_sums[i])*(1-(Count_mat[i,j]/col_sums[j])),
                            -(Count_mat[i,j]*Count_mat[k,j])/(row_sums[i]*col_sums[j]))
    }
  }
}



################################################################################

################################################################################

#tuning constant determined by simulation and comparison to real data. We do not have this so I have simulated dummy tuning constants. This part is key in determining things.......
#gamma<-matrix(rdirichlet(1, rep(1,120*144)),120,144)
#gamma<-matrix(gamma_val,120,dim.para[1])

#calculate P_hat_star

P_hat_star <- (gamma_val*mean(weight_col_means))+(gamma_val*PCF_mat)

#calculate the variance of P_hat

var_P_hat<-(P_hat_star*(1-P_hat_star))/Count_mat

#calculate the variance of D_hat

Var_D_hat<-D_hat<-vector(mode="numeric",length=120)

for (i in 1:120){
  D_hat[i]<-sum(lambda[i,,]*(PCF_mat), na.rm =TRUE)
  Var_D_hat[i]<-sum((lambda[i,,]^2)*var_P_hat, na.rm =TRUE)
}


SE_D <- sqrt(Var_D_hat)

z<-D_hat/SE_D

myresult=cbind(row_sums,D_hat,z)
myresult=myresult[order(row_sums),]
colnames(myresult)=c('N','Diff','z')

plot(myresult[,1],myresult[,3],pch=20,main="TEF Benchmark (each point is a different institution)",xlab="Number of students at institution (N)",ylab="Z scored benchmarked difference")
abline(h=0)
# Stop the clock
proc.time() - ptm

#stopCluster(cl) 

