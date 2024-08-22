setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/glucose_wgd/")

##handles both ploidies at the same time.
## assumes seperate bystander killings per experiment,
## identical glucose and R per ploidy (but not experiment)

init_file <- function(writecon,model,LR){
  source(paste0("Rscripts/models/model",model,".R"))
  source("Rscripts/utils.R")
  dat <- readRDS("data/fitting/fit_dat.Rds")
  info <- model_info()
  icpar <- get_ics_all(dat,LR=LR)
  
  ibys <- which(info$parnames=="kbys")
  maxbys <- info$upper[ibys]
  minbys <- info$lower[ibys]
  
  info$parnames <- c(info$parnames[-ibys],"es_exp1_kbys","es_exp2_kbys")
  info$upper <- c(info$upper[-ibys],maxbys,maxbys)
  info$lower <- c(info$lower[-ibys],minbys,minbys)
  
  iN <- grepl("ic_N",icpar$parnames)
  
  info$parnames <- c(info$parnames,icpar$parnames[iN])
  icpar$parnames <- icpar$parnames[!iN]
  
  parNames <- c(paste0("p2N_",info$parnames),paste0("p4N_",info$parnames),icpar$parnames,"value")
  write(paste(parNames,collapse=","),writecon,append=F)
  return(0)
}

wrap_opt <- function(model,writecon,Ntrial=100,LR=FALSE){
  library(lhs)
  source(paste0("Rscripts/models/model",model,".R"))
  source("Rscripts/utils.R")
  dat <- readRDS("data/fitting/fit_dat.Rds")
  gs <- readRDS("data/fitted_parameters/glucose_struct.Rds")
  info <- model_info()
  icpar <- get_ics_all(dat,LR=LR)
  
  ibys <- which(info$parnames=="kbys")
  maxbys <- info$upper[ibys]
  minbys <- info$lower[ibys]
  
  info$parnames <- c(info$parnames[-ibys],"es_exp1_kbys","es_exp2_kbys")
  info$upper <- c(info$upper[-ibys],maxbys,maxbys)
  info$lower <- c(info$lower[-ibys],minbys,minbys)
  
  iN <- grepl("ic_N",icpar$parnames)
  
  info$parnames <- c(info$parnames,icpar$parnames[iN])
  icpar$parnames <- icpar$parnames[!iN]
  
  info$upper <- c(info$upper,icpar$upper[iN])
  icpar$upper <- icpar$upper[!iN]
  
  info$lower <- c(info$lower,icpar$lower[iN])
  icpar$lower <- icpar$lower[!iN]
    
  parmax <- log(c(info$upper,info$upper,icpar$upper))
  parmin <- log(c(info$lower,info$upper,icpar$lower))
  parNames <- c(paste0("p2N_",info$parnames),paste0("p4N_",info$parnames),icpar$parnames)
  par_range <- parmax-parmin
  
  init_pars <- randomLHS(Ntrial,length(parmax))
  colnames(init_pars) <- parNames
  init_pars <- t(apply(init_pars,1, function(ppi){
    ppi*par_range + parmin
  }))
  
  x0 <- apply(init_pars,1,function(p0){
    masterfit_all(pars=p0,dat=dat,parNames=parNames,gs = gs)  
  })
  
  best <- as.numeric(init_pars[x0==min(x0),])
  opt <- optim(best,fn = masterfit_all,dat=dat,parNames=parNames,gs = gs)
  
  res <- c(opt$par,opt$value)
  names(res) <- c(parNames,"value")
  write(paste(res,collapse=","),writecon,append=T)
  return(res)
}
library(parallel)
ncores <- 70
Ntrial <- 200
Nstarts <- 10000
cl <- makeCluster(getOption("cl.cores", ncores))
outdir <- "data/fitted_parameters/"
dir.create(outdir)

writecon <- paste0(outdir,"opt_1LR.csv")
init_file(writecon,model="1LR",LR=T)
opt_1LR <- do.call(rbind,parLapplyLB(cl,X=rep("1LR",Nstarts),
                                      fun = wrap_opt, writecon=writecon,Ntrial=Ntrial,LR=TRUE))
saveRDS(opt_1LR,file="opt_1LR.Rds")

writecon <- paste0(outdir,"opt_0LR.csv")
init_file(writecon,model="0LR",LR=T)
opt_0LR <- do.call(rbind,parLapplyLB(cl,X=rep("0LR",Nstarts),
                                     fun = wrap_opt, writecon=writecon,Ntrial=Ntrial,LR=TRUE))
saveRDS(opt_0LR,file="opt_1LR.Rds")
