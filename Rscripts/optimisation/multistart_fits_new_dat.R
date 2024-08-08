setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/glucose_wgd/")

wrap_opt <- function(ploidy,model,Ntrial=100){
  library(lhs)
  source(paste0("Rscripts/models/model",model,".R"))
  source("Rscripts/utils.R")
  dat <- readRDS("data/fitting/fitting_data.Rds")
  dat <- melt_data_for_plotting(dat)
  gpar <- readRDS("data/fitted_parameters/glu_pars.Rds")
  info <- model_info()
    
  parmax <- log(c(info$upper,0.3))
  parmin <- log(c(info$lower,0.0000001))
  parNames <- c(info$parnames,"delta_G")
  par_range <- parmax-parmin
  
  init_pars <- randomLHS(Ntrial,length(parmax))
  colnames(init_pars) <- parNames
  init_pars <- t(apply(init_pars,1, function(ppi){
    ppi*par_range + parmin
  }))
  
  x0 <- apply(init_pars,1,function(p0){
    fit_full_model(p0,dat=dat$bx,parNames=parNames,ploidy=ploidy,glu_dat=dat$gx, theta_hat=NULL,gpar=gpar)  
  })
  
  best <- as.numeric(init_pars[x0==min(x0),])
  opt <- optim(best,fn = fit_full_model,dat=dat$bx,parNames=parNames,ploidy=ploidy,glu_dat=dat$gx,gpar=gpar)
  
  res <- c(opt$par,opt$value)
  names(res) <- c(parNames,"value")
  return(res)
}
library(parallel)
ncores <- 3
Ntrial <- 100
Nstarts <- 3
cl <- makeCluster(getOption("cl.cores", ncores))

opt_2n_b <- do.call(rbind,parLapplyLB(cl,X=rep("2N",Nstarts),
                                      fun = wrap_opt, model="B",Ntrial=Ntrial))
saveRDS(opt_2n_b,file="data/fitted_parameters/new_data/opt_2N_B.Rds")

opt_4n_b <- do.call(rbind,parLapplyLB(cl,X=rep("4N",Nstarts),
                                      fun = wrap_opt, model="B",Ntrial=Ntrial))
saveRDS(opt_4n_b,file="data/fitted_parameters/new_data/opt_4N_B.Rds")

opt_2n_a3 <- do.call(rbind,parLapplyLB(cl,X=rep("2N",Nstarts),
                                      fun = wrap_opt, model="A3",Ntrial=Ntrial))
saveRDS(opt_2n_a3,file="data/fitted_parameters/new_data/opt_2N_A3.Rds")

opt_4n_a3 <- do.call(rbind,parLapplyLB(cl,X=rep("4N",Nstarts),
                                       fun = wrap_opt, model="A3",Ntrial=Ntrial))
saveRDS(opt_4n_a3,file="data/fitted_parameters/new_data/opt_4N_A3.Rds")







