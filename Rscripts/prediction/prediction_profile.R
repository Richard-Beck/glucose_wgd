#
setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/00_modelling_v2/")



wrap_constrained_fit <- function(zi,p0,time,gluconc,fit_dat,ploidy){
  constrained_fit <- function(par,p0,zi,time,gluconc,fit_dat,ploidy){
    err <- tryCatch({
    lli <- fit_full_model(par,dat=fit_dat$bx,parNames = names(p0),
                          ploidy=ploidy,glu_dat=fit_dat$gx)
    z <- run_full_model(par,names(p0),gluconc = gluconc)
    z <- z$N[z$time==time]
    err <- 10^6*(log(z)-log(zi))^2+lli
    if(length(err)!=1) return(10^12)
    if(is.na(err)) return(10^12)
    err
    },error = function(e) return(10^9))
    #print(err)
    return(err)
  }
  opt <- optimr(p0,constrained_fit,"grcentral",method="BFGS",control=list(trace=0),
         p0=p0,zi=zi,time=time,gluconc=gluconc,fit_dat=fit_dat,ploidy=ploidy)
  
  par <- opt$par
  z <- run_full_model(par,names(p0),gluconc = gluconc)
  z <- z$N[z$time==time]
  lli <- fit_full_model(par,dat=fit_dat$bx,parNames = names(p0),
                        ploidy=ploidy,glu_dat=fit_dat$gx)
  data.frame(z_target=zi,z_actual=z,ll=lli)
}

model <- "modelB"
ploidy <- "2N"

source(paste0("Rscripts/models/",model,".R"))
source("Rscripts/optimisation/fitting_funcs.R")
source("Rscripts/utils.R")
library(lhs)
library(parallel)

fit_dat <- load_dat()
p0 <- readRDS(paste0("best_pars/",model,"/",ploidy,".Rds"))
ll0 <- p0$value
p0 <- p0$par

time <- 24
gluconc <- 0.05
z <- run_full_model(p0,names(p0),gluconc = gluconc)
z <- z$N[z$time==time]
z_range <- seq(10,1400,10)






cl <- makeCluster(getOption("cl.cores", 70))
clusterCall(cl, function(model) {
  library(deSolve)
  library(optimx)
  source(paste0("Rscripts/models/",model,".R"))
  source("Rscripts/optimisation/fitting_funcs.R")
},model=model)


prof <- parLapplyLB(cl=cl,X=z_range,fun = wrap_constrained_fit,p0=p0,time=time,gluconc=gluconc,fit_dat=fit_dat,ploidy=ploidy)
prof <- do.call(rbind,prof)
saveRDS(prof,paste0("predictions/test1_",model,"_",ploidy,".Rds"))