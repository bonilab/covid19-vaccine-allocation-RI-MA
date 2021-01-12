#!/usr/bin/env Rscript

### MA-v5.R
### last edited: 01 Nov 2020
### (edited from Dr. Hank's example script)

# folder to store output Rdata
args = commandArgs(trailingOnly=TRUE) 
if (length(args)==0) {
  output_dir = "output"
} else {
  output_dir = args[1]
}

# Number of mcmc chains to run in parallel (Note: MUST BE LESS THAN 20!)
n.chs = 5

###############################################################################
### 1. Preliminaries
###############################################################################

### necessary packages

## library for vectorized multinomial
library(mc2d)
## library for multivariate normal distribution
library(mvtnorm)
## library to read in .xlsx files (not needed to read in .csv files)
library(readxl)
## library for spline-based expansions
library(fda)
## packages for parallelization:
library(snow); library(doParallel); library(foreach)
## package for multivariate normal
library(msm)

### compile the odesim model
odepath="../../../../cpp-v6-test-vaccination/"
# odepath <- "../../../../cpp-v5-discharges-nonhospdeaths/"

## uncomment to recompile odesim
# system(paste("rm ",odepath,"odesim",sep=""))
# system(paste("make -C", odepath))

### read in data
# ma.data <- read.csv("../../../../data/Massachusetts/Past_data/MA_20200906_data_7daysmoothed.csv")
ma.data <- read.csv("../../../../data/Massachusetts/MA_20201206_7daysmoothed.csv")

source("../../../code/data.process.R")
dp <- data.process(ma.data, loc="MA")

## removing days before 61
idx.remove <- which(ma.data$daynum < 61)
if(length(idx.remove) > 0){
  ma.data <- ma.data[-idx.remove,]
}

n.days <- nrow(ma.data)
days <- ma.data$daynum

### mcmc initialization (don't change)
end.day <- max(days, na.rm = TRUE) + 1
num.days <- end.day - 60 ## number of days with data

### Create spline expansions

# Cubic spline with one basis function every 7 days
bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days/7))
# zero-order spline with one rate
bspl.rr <- create.bspline.basis(c(61,end.day), nbasis = 1, norder = 1)

# create basis evaluation function matrices
Z <- eval.basis(bspl, 61:end.day)
Z.rr <- eval.basis(bspl.rr, 61:end.day)


### Creating testing delay probabilities

# list of 4 vectors, each corresponding to a reporting period:
#   p1: Beginning - March 13,
#   p2: March 14 - March 23,
#   p3: March 24 - March 27,
#   p4: March 28 - present
# the 1st term = prob of 0 day delay, 2nd = prob. of 1 day delay, etc.
delay.probs <- list(p1 = c(1,rep(0,7)),
                    p2 = c(1,rep(0,7)),
                    p3 = c(1,rep(0,7)),
                    p4 = c(1,rep(0,7)))


###############################################################################
### 2. Starting values for parameters
###############################################################################

# starting values for beta (length should == ncol(Z))
beta.strt <- c(1.96598288105542, 1.01826909513005, 0.919589388966286, 0.519400789304884, 0.155984509495529, 
               0.26363296888402, 0.06180450536809, 0.08448471975314, 0.076203176346524, 0.054456240146284, 
               0.140994710559201, 0.003733955067284, 0.127610988369687, 0.086681814165828, 0.127108018854174, 
               0.069998388194707, 0.165285424498244, 0.128451599421262, 0.151316269260835, 0.120111187813653, 
               0.192210592571, 0.140019758898831, 0.192710792859994, 0.081526797652464, 0.165577516101401, 
               0.160103839640496, 0.209316701755828, 0.157663005363419, 0.151558461223054, 0.214848166058495,
               0.122263232752122, 0.201642710421095, 0.252261825964921, 0.131641732458921, 0.160019439075722, 
               0.495662603263847, 0.209316701755828, 0.357663005363419, 0.351558461223054, 0.214848166058495)


## starting values for odesim parameters and other parameters
prms <- c(2.9846698, 1.6669738, 2.7871050, 2.9889799, 2.9384883, 
          2.5081480, 2.1894337, 3.5101630, 7.0658870, 
          2.9975107, 4.2186056, 2.1276788, 2.0301146, 1.2005998, 
          1.1127662, 1.6358084, 1.2231795, 4.1241666, 
          135.12482220369,
          7.19868338, 0.11060653, 0.34329851, 
          3.97273602, 1.11190707, 0.94737756, 0.62379872, 166.58897982, 
          0.78035244, 0.84262255,
          0.03417957, 0.03211873, 0.05652323, 0.08024622, 0.11088909, 0.20420583, 0.29560123, 0.21313143)

names(prms) <- c("contact-coeff-00","contact-coeff-10","contact-coeff-20","contact-coeff-30","contact-coeff-40",
                 "contact-coeff-50","contact-coeff-60","contact-coeff-70","contact-coeff-80",
                 "contact-coeff-postld-00","contact-coeff-postld-10","contact-coeff-postld-20","contact-coeff-postld-30","contact-coeff-postld-40",
                 "contact-coeff-postld-50","contact-coeff-postld-60","contact-coeff-postld-70","contact-coeff-postld-80",
                 "firstlockdown-endday",
                 "mean-time-vent", "death-prob-home-70", "death-prob-home-80", 
                 "time-symp-to-hosp", "dev-len-hospstay", "dev-icu-frac", "dev-icu-frac-phase2", "dev-icu-frac-phase2beginday",
                 "prob-icu-vent", "dev-ventdeath-mid", 
                 "hosp-frac-10", "hosp-frac-20", "hosp-frac-30", "hosp-frac-40", "hosp-frac-50", "hosp-frac-60", "hosp-frac-70", "hosp-frac-80")


## prior limits for other parameters (all uniform priors)
params.prior.min=c(.1, .1, .1, .1, .1, .1, .1, .1, .1,
                   .1, .1, .1, .1, .1, .1, .1, .1, .1,
                   125,
                   5, 0.01, 0.1, 
                   2, 0.4, 0.1, 0.1, 100, 
                   0.5, 0.6, 
                   0.01, 0.01, 0.01, 0.01, 0.01, 0.15, 0.15, 0.15)

params.prior.max=c(10, 10, 10, 10, 10, 10, 10, 10, 20,
                   10, 10, 10, 10, 10, 10, 10, 10, 20,
                   200,
                   15, 0.4, 0.5, 
                   10, 1.5, 1.5, 1.5, 200, 
                   0.9, 1.5, 
                   0.10, 0.10, 0.10, 0.15, 0.20, 0.30, 0.40, 0.40)

## rough calculation for variance
range=params.prior.max-params.prior.min
ode.params.var.start=range^2/(mean(range^2))
names(ode.params.var.start)=names(prms)
## changing this for "time-symp-to-hosp"
#ode.params.var.start["time-symp-to-hosp"]=1

## starting values for non ode and constant parameters
non.odesim.params = NULL
## medians from MA-1110-nov06S run
const.prms.start=c(2, "")
names(const.prms.start) = c("steps-per-day", "symp-frac-davies")
# const.prms.start=c(2, "", 1.5, 2.4, 1.97, 1.2, 1.2, 0.605, 0.638, 1.1)
# names(const.prms.start) = c("steps-per-day", "symp-frac-davies", "contact-rate-postld-10","contact-rate-postld-20","contact-rate-postld-30","contact-rate-postld-40", "contact-rate-postld-50","contact-rate-postld-60","contact-rate-postld-70","contact-rate-postld-80")

## starting values for NB dispersion parameters
# nb.disp.start=c(.3,.3,.3) # first = new case counts, second = new hosp. cases, third = total deaths (if lik.tot.deaths=TRUE)
# nb.disp.start=c(4.844052, 9.192843, 19.860591)

############################################################
##
## Empirical start values for NB dispersion parameters
## (added 8/21/20)
############################################################


# par(mfrow=c(2,2))
## sympt
sympt.loess = loess(dp$tot.sympt.new~dp$days)
sympt.hat=predict(sympt.loess)
resids.sympt=dp$tot.sympt.new[!is.na(dp$tot.sympt.new)]-sympt.hat
days.sympt=dp$days[!is.na(dp$tot.sympt.new)]
## estimating size parameter
sympt.nb.size.hat=1/median((resids.sympt^2-sympt.hat)/sympt.hat^2)
## check with plots and simulations
mean.sympt=sympt.hat
mean.sympt[mean.sympt<=0]=.0001
# plot(days.sympt,rnbinom(length(days.sympt),mu=mean.sympt,size=sympt.nb.size.hat),type="p",main="sympt")
# points(days.sympt,sympt.hat,type="l")
# points(dp$days,dp$tot.sympt.new,type="b",col="red")

## hosp
hosp.loess = loess(dp$tot.hosp.new~dp$days)
hosp.hat=predict(hosp.loess)
resids.hosp=dp$tot.hosp.new[!is.na(dp$tot.hosp.new)]-hosp.hat
days.hosp=dp$days[!is.na(dp$tot.hosp.new)]
## estimating size parameter
hosp.nb.size.hat=1/median((resids.hosp^2-hosp.hat)/hosp.hat^2)
## check with plots and simulations
mean.hosp=hosp.hat
mean.hosp[mean.hosp<=0]=.0001
# plot(days.hosp,rnbinom(length(days.hosp),mu=mean.hosp,size=hosp.nb.size.hat),type="p",main="hosp")
# points(days.hosp,hosp.hat,type="l")
# points(dp$days,dp$tot.hosp.new,type="b",col="red")

## deaths
deaths.loess = loess(dp$tot.deaths.new~dp$days)
deaths.hat=predict(deaths.loess)
resids.deaths=dp$tot.deaths.new[!is.na(dp$tot.deaths.new)]-deaths.hat
days.deaths=dp$days[!is.na(dp$tot.deaths.new)]
## estimating size parameter
deaths.nb.size.hat=1/median((resids.deaths^2-deaths.hat)/deaths.hat^2)
## check with plots and simulations
mean.deaths=deaths.hat
mean.deaths[mean.deaths<=0]=.0001
# plot(days.deaths,rnbinom(length(days.deaths),mu=mean.deaths,size=deaths.nb.size.hat),type="p",main="deaths")
# points(days.deaths,deaths.hat,type="l")
# points(dp$days,dp$tot.deaths.new,type="b",col="red")




## starting values for NB dispersion parameters
## first = new case counts
## second = new hospitalization counts
## third = total deaths (if lik.tot.deaths==TRUE)
## third = hosp deaths (if lik.hosp.deaths==TRUE and lik.home.deaths==TRUE)
## fourth = home deaths (if lik.hosp.deaths==TRUE and lik.home.deaths==TRUE)
## fifth = discharges
## ...so nb.disp.start with have different lengths depending on what data are used
nb.disp.start=c(sympt.nb.size.hat,hosp.nb.size.hat,deaths.nb.size.hat)
nb.disp.start

## starting values for reporting rates
rr.start <-  0.734666488082196
rr.daily <- Z.rr %*% rr.start

# ## checking visual fit
# source("../../traj.process.R")
# source("../../traj.from.params.R")
# source("../../plots.odesim.R")
# tp=traj.process(traj.from.params(beta=beta.strt,params=prms,const.params=const.prms.start,non.odesim.params=non.odesim.params,
#                                  loc="MA",odepath=odepath,tf=end.day),loc="MA",odesim.version="v5")
# par(mfrow=c(2,3))
# plots.odesim(dp,tp,rr=rr.daily[-1])
# plot(tp$tot.hosp.new.odesim)

# create a "list" of starting values for all chains
betas.0=list()
ode.params.0=list()
rr.params.0=list()
lik.params.0=list()
s2.hosp.0=list()
s2.vent.0=list()
s2.icu.0=list()
for(k in 1:n.chs){
  betas.0[[k]]=beta.strt
  ode.params.0[[k]]=prms
  rr.params.0[[k]]= rr.start # rep(.7,ncol(Z.rr))
  lik.params.0[[k]]=nb.disp.start
  s2.hosp.0[[k]]=1
  s2.vent.0[[k]]=1
  s2.icu.0[[k]]=1
}


##
## proposal variance object
##
## order= 
var.tune.0=c(.0001,.001,.007,25,.01)
Sigma.tune.0=list(Sigma.tune.beta=diag(length(betas.0[[1]])),
                  Sigma.tune.ode.params=diag(ode.params.var.start),
                  Sigma.tune.rr.params=diag(length(rr.params.0[[1]])),
                  Sigma.tune.s2=diag(3),
                  Sigma.tune.ll.params=diag(length(lik.params.0[[1]]))
)


source("../../traj.from.params.R")
source("../../loglik.odesim.4.0.R", chdir=TRUE)
source("../../mcmc.odesim.2.0.R")
source("../../traj.process.R")
source("../../data.process.R")
source("../../plots.odesim.R")

# fit <- mcmc.odesim(n.mcmc=300,
#                    beta.start=betas.0[[1]],
#                    report.rate.params.start=rr.params.0[[1]],
#                    ode.params.start=ode.params.0[[1]],
#                    s2.hosp.start=s2.hosp.0[[1]],
#                    s2.vent.start=s2.vent.0[[1]],
#                    s2.icu.start=s2.icu.0[[1]],
#                    lik.params.start=lik.params.0[[1]],
#                    Sigma.tune=Sigma.tune.0,
#                    var.tune=var.tune.0,
#                    ## mcmc.odesim options below
#                    df=ma.data, ## data
#                    loc="MA",
#                    lik.tot = FALSE, ## only use total new case data, not age-structured
#                    lik.age = TRUE, ## use age-structured new case data
#                    lik.hosp.new = FALSE, ## use daily new hospitalization data
#                    lik.hosp.curr = TRUE, ## use current hospitalization data
#                    lik.icu.curr = TRUE, ## use current ICU data
#                    lik.vent.curr = TRUE, ## use current on-vent data
#                    lik.tot.deaths = FALSE, ## use total deaths
#                    lik.home.deaths = FALSE,
#                    lik.hosp.deaths = FALSE,
#                    lik.age.deaths = TRUE,
#                    total.size.constraint = TRUE,
#                    fixed.nb.disp = TRUE,
#                    spline.beta=bspl, ## fda spline object for time-varying betas
#                    spline.rr=bspl.rr, ## fda spline object for time-varying betas
#                    ode.params.prior.min=params.prior.min, ## prior minima  odesim inputs
#                    ode.params.prior.max=params.prior.max, ## prior maxima  odesim inputs
#                    start.day=61,
#                    end.day=end.day,
#                    introday=55, ## timing info
#                    odepath=odepath, ## path to odesim
#                    odesim.ver="v6",
#                    p.vecs=delay.probs, ## testing delay probabilities
#                    thin=10, ## thinning rate for saving posterior samples
#                    plot.rate=100,
#                    print.iter=TRUE,
#                    c0=1,
#                    c1=.8,
#                    adapt.type="ShabyWells",
#                    adapt.iter=100,
#                    # hosp.report.rate=0.4,
#                    const.params=const.prms.start,
#                    non.odesim.params=non.odesim.params,
#                    sf.choice=FALSE
# )


###############################################################################
### 3. MCMC chains in parallel
###############################################################################

n.iter = 500
n.mcmc.per.iter = 1000

# library(doMC)
# registerDoMC(cores=n.chs)

source("../../multichain.mcmc.odesim.2.0.R")

### Run MCMC
fit <- multichain.mcmc.odesim(parallel.type = "psock", ## "psock" for aci or "doMC"
                              n.chains = n.chs, ## number of independent MCMC chains
                              n.cores = n.chs, ## number of cores to use (keep equal to n.chains)
                              n.iter = n.iter, ## number of times to repeat n.mcmc iterations
                              n.mcmc.per.iter = n.mcmc.per.iter, # number of mcmc iterations to run before saving
                              resample = FALSE, ## ignore for now
                              df = ma.data, ## data
                              loc = "MA",
                              save.file.name = paste(output_dir,"out",sep="/"),  ## output file directory
                              inf.file.dir = "../../",
                              odesim.ver = "v6", ## either "v4" or "v5"
                              introday = 55, ## day of start of epidemic
                              betas.lst = betas.0, ## starting values
                              ode.params.lst = ode.params.0, ## starting values
                              lik.params.lst = lik.params.0, ## starting values
                              rr.params.lst = rr.params.0, ## starting values
                              s2.hosp.lst = s2.hosp.0, ## starting values
                              s2.vent.lst = s2.vent.0, ## starting values
                              s2.icu.lst = s2.icu.0, ## starting values
                              Sigma.tune = Sigma.tune.0,
                              var.tune = var.tune.0, ## starting value for proposal variance.  Can adjust
                              lik.tot = FALSE, ## only use total new case data, not age-structured
                              lik.age = TRUE, ## use age-structured new case data
                              lik.hosp.new = FALSE, ## use daily new hospitalization data
                              lik.hosp.curr = TRUE, ## use current hospitalization data
                              lik.icu.curr = TRUE, ## use current ICU data
                              lik.vent.curr = TRUE, ## use current on-vent data
                              lik.tot.deaths = FALSE, ## use total deaths
                              lik.home.deaths = FALSE,
                              lik.hosp.deaths = FALSE,
                              lik.age.deaths = TRUE,
                              lik.hosp.discharges=FALSE,
                              total.size.constraint = TRUE,
                              fixed.nb.disp = TRUE,
                              spline.beta = bspl, ## fda spline object for time-varying betas
                              spline.rr = bspl.rr, ## fda spline object for time-varying reporting rate
                              ode.params.prior.min = params.prior.min, ## prior minima  odesim inputs
                              ode.params.prior.max = params.prior.max, ## prior maxima  odesim inputs
                              start.day = 61, ## don't change this!
                              end.day = max(ma.data$daynum,na.rm=TRUE)+1, ## timing info
                              odepath = odepath, ## path to odesim
                              p.vecs = delay.probs, ## testing delay probabilities
                              adapt = "ShabyWells", ## adaptive tuning type.  "None" means no tuning.
                              adapt.iter = 100, ## number of iterations between tuning (100 or 1000 is good)
                              thin.rt = 10, ## only save every "thin.rt"-th iteration
                              plot.save = FALSE, ## write out plots occasionally,
                              const.params = const.prms.start,
                              non.odesim.params = non.odesim.params,
                              sf.choice = FALSE
)
