
setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/00_modelling_v2/")

run_model <- "modelA3.R"
source_model <- "modelA"

print("enter ploidy, 1=2N; 2=4N: ")
ploidy <- c("2N","4N")[as.numeric(readLines(con = "stdin", n = 1))]

#ploidy <- "2N"
comments <- ploidy

print("enter number of starts: ")
nstarts <- as.numeric(readLines(con = "stdin", n = 1))
#nstarts <- 3

print("enter number of cores: ")
ncores <- as.numeric(readLines(con = "stdin", n = 1))

#ncores <- 3

source(paste0("Rscripts/models/",run_model))
source("Rscripts/optimisation/fitting_funcs.R")
source("Rscripts/utils.R")
library(lhs)
library(parallel)

fit_dat <- load_dat()
info <- model_info()

p0 <- data.frame(readRDS(paste0("best_pars/",source_model,"/",ploidy,".Rds"))$res)
p0 <- head(p0[order(p0$value),],nstarts)

p0 <- p0[,!colnames(p0)%in%c("value","fevals","gevals","convergence")]
p0$delta_G <- -3
p0 <- p0[,info$parnames]

p0 <- lapply(1:nrow(p0), function(i) unlist(p0[i,]))

#fit_full_model(p0[[1]],dat=fit_dat$bx,parNames=info$parnames,ploidy=ploidy,glu_dat=fit_dat$gx, theta_hat=NULL)


wrap_opt <- function(p0,fit_dat,info,ploidy){
  optimr(p0,fit_full_model,"grcentral",lower=log(info$lower),upper=log(info$upper),method="L-BFGS-B",control=list(trace=0),dat=fit_dat$bx,parNames=info$parnames,ploidy=ploidy,glu_dat=fit_dat$gx, theta_hat=NULL)
}
cl <- makeCluster(getOption("cl.cores", ncores))
clusterCall(cl, function(model) {
  library(deSolve)
  library(optimx)
  source(paste0("Rscripts/models/",model))
  source("Rscripts/optimisation/fitting_funcs.R")
},model=run_model)
 
#v0 <- apply(p0,1,function(p0i) fit_full_model(p0i,dat=fit_dat$bx,parNames=names(info$parnames),ploidy=ploidy,glu_dat=fit_dat$gx))









eval_p0 <- function(p0,fit_dat,info,ploidy){
  fit_full_model(p0,dat=fit_dat$bx,parNames=info$parnames,ploidy=ploidy,glu_dat=fit_dat$gx, theta_hat=NULL)
}



opt <- parLapplyLB(cl=cl,X=p0,fun = wrap_opt,fit_dat=fit_dat,info=info,ploidy=ploidy)
  
res <- do.call(rbind,lapply(opt, function(oi) c(oi$par,oi$value,oi$counts,oi$convergence)))
  
colnames <- sapply(1:length(info$parnames), function(i) paste0("p",i))
colnames <- c(colnames,"value","fevals","gevals","convergence")
  
colnames(res) <- colnames


colnames(res)[1:length(info$parnames)] <- info$parnames

out <- list(res=res,comments=comments)
model <- unlist(strsplit(run_model,"[.]"))[1]

saveid <- length(list.files(paste0("opt_out/",model)))
library(stringr)
saveid <- str_pad(saveid,width=3,pad=0)

fnameout <- paste0("opt_out/",model,"/",saveid,".Rds")
saveRDS(out,fnameout)


