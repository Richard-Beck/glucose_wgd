#
setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/00_modelling_v2/")


wrap_mix_err <- function(zi,p0i,parnames,time,G,fit_dat){
  mixmod_err <- function(par,zi,parnames,time,G,fit_dat){
    err <- tryCatch({
      p02 <- par[1:(length(par)/2)]
      p04 <- par[(length(par)/2+1):length(par)]
      names(p02) <- parnames
      names(p04) <- parnames
      
      p0 <- list(p02,p04)
      y0 <- list(c(N=500,D=0),c(N=500,D=0))
      
      err2n <- fit_full_model(pars=p02,dat=fit_dat$bx,parNames=parnames,ploidy="2N",glu_dat = fit_dat$gx)
      err4n <- fit_full_model(pars=p04,dat=fit_dat$bx,parNames=parnames,ploidy="4N",glu_dat = fit_dat$gx)
      
      z <- run_mixmod(p0,times=1:time,y0=y0,G=G)
      z <- as.numeric(z[nrow(z),"N1"]-z[nrow(z),"N2"])
      errz <- 10^6*((z-zi)/zi)^2
      
      err <- err2n+err4n+errz
      if(length(err)!=1) stop()
      if(!is.finite(err)) stop()
      err
    },error = function(e) return(10^12))
    return(err)
  }
  optimr(p0i,mixmod_err,"grcentral",method="BFGS",control=list(trace=0),
         zi=zi,parnames=parnames,time=time,G=G,fit_dat=fit_dat)
}


model <- "modelB"
source(paste0("Rscripts/models/",model,".R"))
source("Rscripts/optimisation/fitting_funcs.R")
source("Rscripts/utils.R")
library(lhs)
library(parallel)
library(optimx)
p2N <- readRDS("best_pars/modelB/2N.Rds")$par
p4N <- readRDS("best_pars/modelB/4N.Rds")$par
fit_dat <- load_dat()

p0 <- list(p2N,p4N)
y0 <- list(c(N=500,D=0),c(N=500,D=0))
time <- 72
G <- 1

z <- run_mixmod(p0,times=1:time,y0=y0,G=G)
z <- as.numeric(z[nrow(z),"N1"]-z[nrow(z),"N2"])

z_range <- seq(0,500,5)

z1 <- z-z_range
z2 <- z+z_range

opt_par <- unlist(p0)

opt <- list()
opt_smm <- data.frame(z_target=numeric(0),z_actual=numeric(0),ll=numeric(0))
for(zi in z1){
  opti <- wrap_mix_err(zi,opt_par,parnames=names(p2N),time=time,G=G,fit_dat=fit_dat)
  opt <- c(opt,list(opti))
  opt_par <- opti$par
  p02 <- opt_par[1:(length(opt_par)/2)]
  p04 <- opt_par[(length(opt_par)/2+1):length(opt_par)]
  names(p02) <- names(p2N)
  names(p04) <- names(p4N)
  p0 <- list(p02,p04)
  y0 <- list(c(N=500,D=0),c(N=500,D=0))
  z <- run_mixmod(p0,times=1:time,y0=y0,G=G)
  z <- as.numeric(z[nrow(z),"N1"]-z[nrow(z),"N2"])
  errz <- 10^6*((z-zi)/zi)^2
  ll <- opti$value-errz
  dfi <- data.frame(z_target=zi,z_actual=z,ll=ll)
  opt_smm <- rbind(opt_smm,dfi)
  print(opt_smm)
  result <- list(opt=opt,opt_smm=opt_smm)
  saveRDS(result,"predictions/72h_delta_type_modelB_run2.Rds")
}

p0 <- list(p2N,p4N)
y0 <- list(c(N=500,D=0),c(N=500,D=0))
opt_par <- unlist(p0)

for(zi in z2){
  opti <- wrap_mix_err(zi,opt_par,parnames=names(p2N),time=time,G=G,fit_dat=fit_dat)
  opt <- c(opt,list(opti))
  opt_par <- opti$par
  p02 <- opt_par[1:(length(opt_par)/2)]
  p04 <- opt_par[(length(opt_par)/2+1):length(opt_par)]
  names(p02) <- names(p2N)
  names(p04) <- names(p4N)
  p0 <- list(p02,p04)
  y0 <- list(c(N=500,D=0),c(N=500,D=0))
  z <- run_mixmod(p0,times=1:time,y0=y0,G=G)
  z <- as.numeric(z[nrow(z),"N1"]-z[nrow(z),"N2"])
  errz <- 10^6*((z-zi)/zi)^2
  ll <- opti$value-errz
  dfi <- data.frame(z_target=zi,z_actual=z,ll=ll)
  opt_smm <- rbind(opt_smm,dfi)
  print(opt_smm)
  result <- list(opt=opt,opt_smm=opt_smm)
  saveRDS(result,"predictions/72h_delta_type_modelB_run2.Rds")
}

result <- list(opt=opt,opt_smm=opt_smm)
saveRDS(result,"predictions/72h_delta_type_modelB_run2.Rds")





