setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/glucose_wgd/")



wrap_opt <- function(ploidy,model,Ntrial=100,LR=FALSE){
  library(lhs)
  source(paste0("Rscripts/models/model",model,".R"))
  source("Rscripts/utils.R")
  dat <- readRDS("data/fitting/fit_dat.Rds")
  gs <- readRDS("data/fitted_parameters/glucose_struct.Rds")
  info <- model_info()
  icpar <- get_ics(dat,ploidy,LR=LR)
    
  parmax <- log(c(info$upper,icpar$upper))
  parmin <- log(c(info$lower,icpar$lower))
  parNames <- c(info$parnames,icpar$parnames)
  par_range <- parmax-parmin
  
  init_pars <- randomLHS(Ntrial,length(parmax))
  colnames(init_pars) <- parNames
  init_pars <- t(apply(init_pars,1, function(ppi){
    ppi*par_range + parmin
  }))
  
  x0 <- apply(init_pars,1,function(p0){
    masterfit(pars=p0,dat=dat,parNames=parNames,gs = gs,ploidy = ploidy)  
  })
  
  best <- as.numeric(init_pars[x0==min(x0),])
  opt <- optim(best,fn = masterfit,dat=dat,parNames=parNames,gs = gs,ploidy = ploidy)
  
  res <- c(opt$par,opt$value)
  names(res) <- c(parNames,"value")
  return(res)
}
library(parallel)
ncores <- 70
Ntrial <- 100
Nstarts <- 20000
cl <- makeCluster(getOption("cl.cores", ncores))
outdir <- "data/fitted_parameters/"
dir.create(outdir)

opt_2n_0LR <- do.call(rbind,parLapplyLB(cl,X=rep("2N",Nstarts),
                                      fun = wrap_opt, model="0LR",Ntrial=Ntrial,LR=TRUE))
saveRDS(opt_2n_0LR,file=paste0(outdir,"/opt_2N_0LR.Rds"))

opt_4n_0LR <- do.call(rbind,parLapplyLB(cl,X=rep("4N",Nstarts),
                                      fun = wrap_opt, model="0LR",Ntrial=Ntrial,LR=TRUE))
saveRDS(opt_4n_0LR,file=paste0(outdir,"/opt_4N_0LR.Rds"))
stop()

opt_2n_0 <- do.call(rbind,parLapplyLB(cl,X=rep("2N",Nstarts),
                                      fun = wrap_opt, model="0",Ntrial=Ntrial))
saveRDS(opt_2n_0,file=paste0(outdir,"/opt_2N_0.Rds"))

opt_4n_0 <- do.call(rbind,parLapplyLB(cl,X=rep("4N",Nstarts),
                                      fun = wrap_opt, model="0",Ntrial=Ntrial))
saveRDS(opt_4n_0,file=paste0(outdir,"/opt_4N_0.Rds"))
stop()
opt_2n_b <- do.call(rbind,parLapplyLB(cl,X=rep("2N",Nstarts),
                                      fun = wrap_opt, model="B",Ntrial=Ntrial))
saveRDS(opt_2n_b,file=paste0(outdir,"/opt_2N_B.Rds"))

opt_4n_b <- do.call(rbind,parLapplyLB(cl,X=rep("4N",Nstarts),
                                      fun = wrap_opt, model="B",Ntrial=Ntrial))
saveRDS(opt_4n_b,file=paste0(outdir,"/opt_4N_B.Rds"))

opt_2n_a3 <- do.call(rbind,parLapplyLB(cl,X=rep("2N",Nstarts),
                                      fun = wrap_opt, model="A",Ntrial=Ntrial))
saveRDS(opt_2n_a3,file=paste0(outdir,"/opt_2N_A.Rds"))

opt_4n_a3 <- do.call(rbind,parLapplyLB(cl,X=rep("4N",Nstarts),
                                       fun = wrap_opt, model="A",Ntrial=Ntrial))
saveRDS(opt_4n_a3,file=paste0(outdir,"/opt_4N_A.Rds"))







