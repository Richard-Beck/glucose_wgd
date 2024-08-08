#

setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/glucose_wgd/")
available_models <- list.files("Rscripts/models/")
model_df <- cbind(available_models,1:length(available_models))
colnames(model_df) <- c("model_name","code")
print(model_df)
print("Enter model code: ")
model <- available_models[as.numeric(readLines(con = "stdin", n = 1))]

print(paste("using",model))

print("enter any comments about this fit that may be helpful later: ")
comments <- readLines(con = "stdin", n = 1)


source(paste0("Rscripts/models/",model))
source("Rscripts/optimisation/fitting_funcs.R")
source("Rscripts/utils.R")
library(lhs)
library(parallel)

fit_dat <- load_dat()
info <- model_info()

print("enter number of starts: ")
nstarts <- as.numeric(readLines(con = "stdin", n = 1))



#v0 <- apply(p0,1,function(p0i) fit_full_model(p0i,dat=fit_dat$bx,parNames=names(info$parnames),ploidy=ploidy,glu_dat=fit_dat$gx))

print("enter ploidy, 1=2N; 2=4N: ")
ploidy <- c("2N","4N")[as.numeric(readLines(con = "stdin", n = 1))]


print("enter number of cores: ")
ncores <- as.numeric(readLines(con = "stdin", n = 1))




eval_p0 <- function(p0,fit_dat,info,ploidy){
  fit_full_model(p0,dat=fit_dat$bx,parNames=info$parnames,ploidy=ploidy,glu_dat=fit_dat$gx, theta_hat=NULL)
}

wrap_opt <- function(p0,fit_dat,info,ploidy){
    optimr(p0,fit_full_model,"grcentral",lower=log(info$lower),upper=log(info$upper),method="L-BFGS-B",control=list(trace=0),dat=fit_dat$bx,parNames=info$parnames,ploidy=ploidy,glu_dat=fit_dat$gx, theta_hat=NULL)
}
cl <- makeCluster(getOption("cl.cores", ncores))
clusterCall(cl, function(model) {
  library(deSolve)
  library(optimx)
  source(paste0("Rscripts/models/",model))
  source("Rscripts/optimisation/fitting_funcs.R")
},model=model)


p0 <- lhs_guess(nstarts,log(info$upper),log(info$lower))
p0 <- lapply(1:nrow(p0), function(i) p0[i,])

ll0 <- unlist(parLapplyLB(cl=cl,X=p0,fun = eval_p0,fit_dat=fit_dat,info=info,ploidy=ploidy))

p0 <- do.call(rbind,p0)

p0 <- p0[order(ll0),]
p0 <- p0[1:nstarts,]

p0 <- lapply(1:nrow(p0), function(i) p0[i,])
opt <- parLapplyLB(cl=cl,X=p0,fun = wrap_opt,fit_dat=fit_dat,info=info,ploidy=ploidy)
  
res <- do.call(rbind,lapply(opt, function(oi) c(oi$par,oi$value,oi$counts,oi$convergence)))
  
colnames <- sapply(1:length(info$parnames), function(i) paste0("p",i))
colnames <- c(colnames,"value","fevals","gevals","convergence")
  
colnames(res) <- colnames


colnames(res)[1:length(info$parnames)] <- info$parnames

out <- list(res=res,comments=comments)
model <- unlist(strsplit(model,"[.]"))[1]

saveid <- length(list.files(paste0("opt_out/",model)))
library(stringr)
saveid <- str_pad(saveid,width=3,pad=0)

dir.create(paste0("opt_out/",model),recursive=T)
fnameout <- paste0("opt_out/",model,"/",saveid,".Rds")
saveRDS(out,fnameout)


