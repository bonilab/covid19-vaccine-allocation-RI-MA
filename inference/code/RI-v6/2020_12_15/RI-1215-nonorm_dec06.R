#!/usr/bin/env Rscript

### RI-v5.R
### last edited: 31 July 2020
### Ephraim Hanks, Nathan Wikle, Emily Strong

# folder to store output Rdata
args = commandArgs(trailingOnly=TRUE) 
if (length(args)==0) {
  output_dir = "output"
} else {
  output_dir = args[1]
}

### Number of mcmc chains to run in parallel (Note: MUST BE LESS THAN 20!)
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
## package for iSplines:
library(splines2)

### compile the odesim model
odepath="../../../../cpp-v6-test-vaccination/"
# odepath="../../../../cpp-v5-discharges-nonhospdeaths/"

# ## uncomment to recompile odesim
# system(paste("rm ",odepath,"odesim",sep=""))
# system(paste("make -C", odepath))

### read in data
ri.data <- read.csv("../../../../data/Rhode_Island/RI_formatted_20201206.csv")
ri.data <- data.frame(ri.data)

## removing days before 61
idx.remove=which(ri.data$daynum<61)
if(length(idx.remove)>0){
  ri.data=ri.data[-idx.remove,]
}

n.days <- nrow(ri.data)
days <- ri.data$daynum

### mcmc initialization (don't change)
end.day <- max(days, na.rm = TRUE) + 1
num.days <- end.day - 60 ## number of days with data

source("../../data.process.R")
dp=data.process(ri.data)
str(dp)

### Create spline expansions

# Cubic spline with one basis function every 7 days
bspl <- create.bspline.basis(c(61, end.day), nbasis=round(num.days/7))


# # cubic spline with multiple rates
# bspl.rr <- create.bspline.basis(c(61,end.day),
#                                 #nbasis=,
#                                 breaks=c(61, 84, 88, 92, 96, 100, 130, 160, 190, end.day),
#                                 norder=4)

### create I-spline basis functions, evaluated every day from 61:end.day
### create I-spline basis functions, evaluated every day from 61:end.day
knots <- c(84, 92, 100, 108, 130, 160, 190)
Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
Z.rr <- Z.rr[,c(1,4:7)]

# create basis evaluation function matrices
Z <- eval.basis(bspl,61:end.day)
# Z.rr <- eval.basis(bspl.rr,61:end.day)

### Creating testing delay probabilities

# list of 4 vectors, each corresponding to a reporting period:
#   p1: Beginning - March 13,
#   p2: March 14 - March 23,
#   p3: March 24 - March 27,
#   p4: March 28 - present
# the 1st term = prob of 0 day delay, 2nd = prob. of 1 day delay, etc.
delay.probs <- list(p1 = c(1, rep(0,7)),
                    p2 = c(1, rep(0,7)),
                    p3 = c(1, rep(0,7)),
                    p4 = c(1, rep(0,7)))


###############################################################################
### 2. Starting values for parameters
###############################################################################

# starting values for reporting rate
rr.strt <- c(0.174444783140871, 0.299576103938855, 0.396107034129121, 0.076464365746264, 0.005996618391096) 
rr.daily = Z.rr %*% rr.strt

# starting values for beta (length should == ncol(Z))
## starting values for beta
# beta.strt <- c(1.9, 1.15, 2.2, 0.45, .45,
#                .35, .16, .17, .14, .05,
#                .13, .08, .05, .12, .07,
#                .08, .12, .2, .2, 0.27, .3, .5)

# beta.strt <- c(0.8, 0.18, 0.09, 0.09, 0.05, 0.06, 0.03, 0.04, 0.05, 0.09,
#                0.1, 0.1, 0.1, 0.1, 0.04, 0.06, 0.14, 0.12, 0.12, 0.19,
#                0.38, 0.5)

beta.strt <- c(1.93390475379047, 0.295310009797186, 0.444736967546712, 0.134063765725156, 0.184334822781263, 
               0.154480036651505, 0.106316468470633, 0.099758437348227, 0.071756257568282, 0.062483738375899, 
               0.050688211668911, 0.07919226761951, 0.060334667262311, 0.352427214003411, 0.130800373081394, 
               0.28967899723136, 0.146940113034839, 0.212396370024068, 0.510531848806016, 0.297908475747728, 
               0.293169714913363, 0.262321690675385, 0.240944797537493, 0.259266714118101, 0.248944344905161, 
               0.32969245770994, 0.286059435558952, 0.430894637998807, 0.260165006242562, 0.434615304451599, 
               0.326060956145053, 0.325714966935963, 0.319484970750116, 0.70015463266141, 0.577914164960418, 
               0.426735088964163, 0.286059435558952, 0.230894637998807, 0.260165006242562, 0.234615304451599)

## starting values for odesim parameters and other parameters
prms.start <- c(2.4511093, 2.5629740, 3.5772877, 3.6135832, 2.2581054, 
                2.8738528, 2.1503628, 3.1708887, 5.3727806, 
                1.4665137, 2.7094735, 1.4205091, 1.2113444, 0.7627325,
                0.5832549, 0.4391605, 0.5036369, 2.2530055,
                150.109718)


names(prms.start) <- c("contact-coeff-00","contact-coeff-10","contact-coeff-20","contact-coeff-30","contact-coeff-40",
                       "contact-coeff-50","contact-coeff-60","contact-coeff-70","contact-coeff-80",
                       "contact-coeff-postld-00","contact-coeff-postld-10","contact-coeff-postld-20","contact-coeff-postld-30","contact-coeff-postld-40",
                       "contact-coeff-postld-50","contact-coeff-postld-60","contact-coeff-postld-70","contact-coeff-postld-80",
                       "firstlockdown-endday")


## ## prior limits for other parameters (all uniform priors)
params.prior.min <- c(.1, .1, .1, .1, .1, .1, .1, .1, .1,
                      .1, .1, .1, .1, .1, .1, .1, .1, .1,
                      140)

params.prior.max <- c(10, 10, 10, 10, 10, 10, 10, 10, 20,
                      10, 10, 10, 10, 10, 10, 10, 10, 20,
                      200)


## names of parameters that are NOT odesim parameters, but are still to be estimated
non.odesim.params <- NULL
## inputs to odesim that are to be FIXED, and not estimated


const.prms.start <- c("", 2,
                      8.9172987, 0.0298327, 0.1432668, 0.294211, 
                      2.3946916, 0.7002739, 0.6303916, 0.4973039, 153.9066363, 
                      0.6585476, 0.7388647, 
                      0.0302321, 0.0349996, 0.0559191, 0.0840281, 0.1220664, 0.2030505, 0.3012512, 0.1981929)
names(const.prms.start) <- c("symp-frac-davies", "steps-per-day", 
                             "mean-time-vent", "death-prob-home-60", "death-prob-home-70","death-prob-home-80",
                             "time-symp-to-hosp","dev-len-hospstay","dev-icu-frac", "dev-icu-frac-phase2", "dev-icu-frac-phase2beginday",
                             "prob-icu-vent", "dev-ventdeath-mid", 
                             "hosp-frac-10","hosp-frac-20","hosp-frac-30","hosp-frac-40","hosp-frac-50", "hosp-frac-60","hosp-frac-70", "hosp-frac-80")


## check to make sure priors and starting values line up
rbind(params.prior.min, prms.start, params.prior.max)

##
##
## checking visual fit (good for choosing prms and beta.strt above)
##
##
source("../../traj.process.R")
source("../../data.process.R")
source("../../traj.from.params.R")
source("../../plots.odesim.R")


## rough calculation for tuning variance - make proposal variance smaller for parameters with smaller prior range
range <- params.prior.max-params.prior.min
ode.params.var.start <- range^2/(mean(range^2))
names(ode.params.var.start) <- names(prms.start)

## starting values for NB dispersion parameters
## first = new case counts
## second = new hospitalization counts
## third = total deaths (if lik.tot.deaths==TRUE)
## third = hosp deaths (if lik.hosp.deaths==TRUE and lik.home.deaths==TRUE)
## fourth = home deaths (if lik.hosp.deaths==TRUE and lik.home.deaths==TRUE)
## ...so nb.disp.start with have different lengths depending on what data are used
nb.disp.start <- c(.3,.3,.3,.3,.3) 

# create a "list" of starting values for all chains
betas.0 <- list()
ode.params.0 <- list()
const.params.0 <- list()
const.params.0 <- const.prms.start
rr.params.0 <- list()
lik.params.0 <- list()
s2.hosp.0 <- list()
s2.vent.0 <- list()
s2.icu.0 <- list()
for(k in 1:n.chs){
  betas.0[[k]] <- beta.strt
  ode.params.0[[k]] <- prms.start
  rr.params.0[[k]] <- rr.strt
  lik.params.0[[k]] <- nb.disp.start
  s2.hosp.0[[k]] <- 1
  s2.vent.0[[k]] <- 1
  s2.icu.0[[k]] <- 1
}


# initialize proposal densities

### initialize proposal variance
var.tune.0=c(.0001,.001,.007,25,.01)
Sigma.tune.0=list(Sigma.tune.beta=diag(length(betas.0[[1]])),
                  Sigma.tune.ode.params=diag(ode.params.var.start),
                  Sigma.tune.rr.params=diag(length(rr.params.0[[1]])),
                  Sigma.tune.s2=diag(3),
                  Sigma.tune.ll.params=diag(length(lik.params.0[[1]]))
)

# var.tune.0 <- c(0.0009628478, 0.0024778840, 0.0210342994, 0.1223803880, 0.0627313899)
# 
# Sigma.tune.beta.0 <- diag(length(beta.strt))
# 
# Sigma.tune.ode.params.0 <- diag(c(1.153230e+01, 5.125468e-02, 5.125468e-02, 2.883076e+00, 1.281367e+00,
#                                   3.203417e-01, 1.153230e-01, 6.278698e-01, 4.625735e-04,
#                                   1.948959e-03, 2.050187e-03, 3.203417e-03, 3.203417e-03, 4.612921e-03,
#                                   3.203417e-03, 3.203417e-03, 
#                                   .01, .01, .01, .01, .01, .01, .01, .01,
#                                   .01, .01, .01, .01, .01, .01, .01, .01))
# 
# Sigma.tune.rr.0 <- diag(ncol(Z.rr))
# 
# Sigma.tune.s2.0 <- diag(3) # number of "current" data streams used (hosp, ICU, vent)
# 
# Sigma.tune.ll.params.0 <- diag(length(nb.disp.start))
# 
# Sigma.tune.0=list(Sigma.tune.beta=Sigma.tune.beta.0,
#                   Sigma.tune.ode.params=Sigma.tune.ode.params.0,
#                   Sigma.tune.rr.params=Sigma.tune.rr.0,
#                   Sigma.tune.s2=Sigma.tune.s2.0,
#                   Sigma.tune.ll.params=Sigma.tune.ll.params.0)



source("../../traj.from.params.R")
source("../../loglik.odesim.4.0.R", chdir=TRUE)
source("../../mcmc.odesim.2.0.R")
source("../../traj.process.R")
source("../../data.process.R")
source("../../plots.odesim.R")

# fit <- mcmc.odesim(n.mcmc=500,
#                    beta.start=betas.0[[1]],
#                    report.rate.params.start=rr.params.0[[1]],
#                    ode.params.start=ode.params.0[[1]],
#                    const.params=const.params.0,
#                    non.odesim.params=non.odesim.params,
#                    s2.hosp.start=s2.hosp.0[[1]],
#                    s2.vent.start=s2.vent.0[[1]],
#                    s2.icu.start=s2.icu.0[[1]],
#                    lik.params.start=lik.params.0[[1]],
#                    Sigma.tune=Sigma.tune.0,
#                    var.tune=var.tune.0,
#                    ## mcmc.odesim options below
#                    df=ri.data, ## data
#                    loc="RI",
#                    lik.tot=FALSE, ## only use total new case data, not age-structured
#                    lik.age=TRUE,
#                    lik.hosp.new=TRUE, ## use daily new hospitalization data
#                    lik.hosp.curr=TRUE,
#                    lik.icu.curr=TRUE,
#                    lik.vent.curr=TRUE,
#                    lik.tot.deaths=FALSE,
#                    lik.age.deaths=TRUE,
#                    lik.home.deaths = TRUE,
#                    lik.hosp.deaths = TRUE,
#                    lik.hosp.discharges=TRUE,
#                    total.size.constraint = TRUE, ## forces sympt to be close to data
#                    active.surv = FALSE,
#                    p.asympt = 0.4, ## change if active.serv if TRUE...
#                    spline.beta=bspl, ## fda spline object for time-varying betas
#                    spline.rr=Z.rr, ## fda spline object for time-varying betas
#                    ode.params.prior.min=params.prior.min, ## prior minima  odesim inputs
#                    ode.params.prior.max=params.prior.max, ## prior maxima  odesim inputs
#                    start.day=61,
#                    end.day=max(ri.data$daynum,na.rm=TRUE)+1,
#                    introday=55, ## timing info
#                    odepath=odepath, ## path to odesim
#                    odesim.ver="v6",
#                    p.vecs=delay.probs, ## testing delay probabilities
#                    thin=1, ## thining rate for saving posterior samples
#                    c0=1,
#                    c1=.8,
#                    adapt.type="ShabyWells",
#                    adapt.iter=100,
#                    s2.beta.start=.01,
#                    s2.rr.start=.01,
#                    t.adapt.start=0,
#                    prop.type="tnorm",
#                    plot.save=TRUE,
#                    plot.rate=200,
#                    print.iter=TRUE,
#                    plot.name="trace.plots.pdf",
#                    sf.choice = FALSE
# )
# str(fit)
# save(fit,file="fit.Rdata")



###############################################################################
### 3. MCMC chains in parallel
###############################################################################

source("../../multichain.mcmc.odesim.2.0.R")

# library(doMC)
# registerDoMC(cores=n.chs)

### Run MCMC
fit <- multichain.mcmc.odesim(parallel.type = "psock", ## "psock" for aci or "doMC"
                              n.chains = n.chs, ## number of independent MCMC chains
                              n.cores = n.chs, ## number of cores to use (keep equal to n.chains)
                              n.iter = 500, ## number of times to repeat n.mcmc iterations
                              n.mcmc.per.iter = 1000, ## number of mcmc iterations to run before saving
                              adapt.iter = 100, ## number of iterations between tuning (100 or 1000 is good)
                              thin.rt = 10, ## only save every "thin.rt"-th iteration
                              resample = FALSE, ## ignore for now
                              df = ri.data, ## data
                              loc = "RI",
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
                              lik.hosp.new = TRUE, ## use daily new hospitalization data
                              lik.hosp.curr = TRUE, ## use current hospitalization data
                              lik.icu.curr = TRUE, ## no current ICU data
                              lik.vent.curr = TRUE, ## No current on-vent data
                              lik.tot.deaths = FALSE, ## total deaths
                              lik.age.deaths = TRUE, ## age-structured deaths
                              lik.home.deaths = TRUE,
                              lik.hosp.deaths = TRUE,
                              lik.hosp.discharges = TRUE, ## use new hosp discharges
                              total.size.constraint = TRUE, ## forces sympt to be close to data
                              active.surv = FALSE,
                              p.asympt = 0.4, ## change if active.serv if TRUE...
                              spline.beta = bspl, ## fda spline object for time-varying betas
                              spline.rr = Z.rr, ## fda spline object for time-varying reporting rate
                              ode.params.prior.min = params.prior.min, ## prior minima  odesim inputs
                              ode.params.prior.max = params.prior.max, ## prior maxima  odesim inputs
                              const.params = const.params.0,
                              non.odesim.params = non.odesim.params,
                              start.day = 61, ## don't change this!
                              end.day = max(ri.data$daynum,na.rm=TRUE) + 1, ## timing info
                              odepath = odepath, ## path to odesim
                              p.vecs = delay.probs, ## testing delay probabilities
                              adapt = "ShabyWells", ## adaptive tuning type.  "None" means no tuning.
                              plot.save = FALSE, ## write out plots occasionally
                              sf.choice = FALSE
)


