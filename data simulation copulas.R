

require(copula)

set.seed(1981)

n <- 120 

# 6 variables:

# years (Gaussian)
# disabled (Binomial)
# ethnicity (Binomial)
# level (Binomial)
# age (Binomial)
# domicile (Gaussian)


## build a normal copula for the underlying distibution:
# `param` is the lower triangle of a correlation matrix
underlying <- normalCopula(dim = 6, 
                           param = c(-0.4,0,-0.2,-0.5,-0.2,   # years
                                     0.3,0.05,0.5,0.7,        # disabled
                                     0.4,0.6,0.3,             # ethnicity
                                     0.3,0.1,                 # level
                                     0.5),                    # age
                           dispstr = "un") # The Copula

## build the marginal distibutions from the underlying copula
marginals <- mvdc(copula = underlying, 
                  margins = c("norm","binom", "binom", "binom", "binom", "norm"), 
                  paramMargins = list(list(mean = 0, sd = 1), # parameters for the marginal distributions...
                                      list(size = 1, prob = 0.5),
                                      list(size = 1, prob = 0.5),
                                      list(size = 1, prob = 0.5),
                                      list(size = 1, prob = 0.5),
                                      list(mean = 0, sd = 1))) 

## generate the correlated random variables from the marginal distributions
dat <- as.data.frame(rMvdc(n, marginals)) 

## rename the variables
names(dat) <- c( "years", "disabled", "ethnicity", "level", "age", "domicile")

## Stratify the years variable
dat$years <- cut(round(dat$years, 5),
                      seq(from=min(dat$years), to = max(dat$years),length.out = 4))

levels(dat$years) <- 1:3

## Stratify the domicile variable
dat$domicile <- cut(round(dat$domicile, 5),
                    seq(from=min(dat$domicile), to = max(dat$domicile),length.out = 4))

levels(dat$domicile) <- 1:3

# convert others to factors and clean up
dat$years <- ordered(dat$years)
dat$domicile <- ordered(dat$domicile)
levels(dat$domicile) <- c("UK", "Other EU","Non EU")
dat$disabled <- factor(dat$disabled)
levels(dat$disabled) <- c("Yes", "No")
dat$ethnicity <- factor(dat$ethnicity)
levels(dat$ethnicity) <- c("White", "BME")
dat$age <- factor(dat$age)
levels(dat$ethnicity) <- c("young", "Mature")
dat$level <- factor(dat$level)
levels(dat$level) <- c("first degree", "Other UG")
dat <- dat[complete.cases(dat),]
