mod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    #if (Time > delay) {
    # G_delayed <- lagvalue(Time - delay, 3)
    #} else {
    # G_delayed <- 25  # If the time is less than the delay, use the initial value
    #}
    dN <- ka*N/(1+(g50a/G)^na)-kd*N*(1-1/(1+(g50d/G)^nd))
    dD <- kd*N*(1-1/(1+(g50d/G)^nd))
    dG <- -N*(kgr*(1/(1+(g50a/G)^na))+kgd*(1/(1+g50d/G)^nd))
    return(list(c(dN,dD,dG)))
  })
}

run_mod <- function(pars,y0){
  times <- seq(0,7,1/30)
  ode(y0,times,mod,pars)
  #dede(y0,times,mod,pars)
}

mod_err <- function(pars,varpars,yi){
  pars <- exp(pars)
  y0 <- c(N=varpars[1],D=0,G=varpars[2])
  
  out <- run_mod(pars,y0)
  
  N <- out[out[,"time"]%in%yi$bx$day,"N"]
  D <- out[out[,"time"]%in%yi$bx$day,"D"]
  
  obs <- c(yi$bx$ncells[yi$bx$metric=="alive"],
           yi$bx$ncells[yi$bx$metric=="dead"])
  sd_obs <- c(yi$bx$sd[yi$bx$metric=="alive"],
              yi$bx$sd[yi$bx$metric=="dead"])
  pred <- c(N,D)
  if(!is.null(yi$gx)){
    G <- out[out[,"time"]%in%yi$gx$day,"G"]
    obs <- c(obs,yi$gx$G)
    sd_obs <- c(sd_obs,yi$gx$sd)
    pred <- c(pred,G)
  }
  pred[is.na(pred)] <- 10^5
  sum((pred-obs)^2/sd_obs^2)
}

wrap_fit_varpars <- function(varpars,pars,yi){
  mod_err(pars,varpars,yi)
}

wrap_fit_mod <- function(pars,yxn){
  N0 <- exp(pars[grepl("inicell",names(pars))])
  G0 <- exp(pars[grepl("iniglu",names(pars))])
  
  pars <- pars[!grepl("inicell",names(pars))]
  pars <- pars[!grepl("iniglu",names(pars))]
  
  errs <- sapply(1:length(yxn),function(i){
    yi <- yxn[[i]]
    #N0 <- yi$bx$ncells[yi$bx$day==min(yi$bx$day)&yi$bx$metric=="alive"]
    varpars <- as.numeric(c(N0[i],G0[i]))
    tryCatch({
      mod_err(pars,varpars,yi)},
      error=function(e) return(10^10))
  })
  print(errs)
  print(sum(errs))
  return(sum(errs))
}

gen_pars <- function(yxn,Ntrial){
  g0 <- do.call(cbind,lapply(yxn,function(yi){
    g0 <- yi$gx[yi$gx$day==0,c("G","sd")]
    log(pmax(0.0001,rnorm(Ntrial,mean=as.numeric(g0["G"]),
                          sd=as.numeric(g0["sd"]))))  
  }))
  colnames(g0) <- paste0("iniglu_",colnames(g0))
  
  n0 <- do.call(cbind,lapply(yxn,function(yi){
    n0 <- yi$bx[yi$bx$day==min(yi$bx$day)&yi$bx$metric=="alive",c("ncells","sd")]
    log(pmax(1,rnorm(Ntrial,mean=as.numeric(n0["ncells"]),
                     sd=as.numeric(n0["sd"]))))  
  }))
  colnames(n0) <- paste0("inicell_",colnames(n0))
  
  
  parmax <- c(kgr=1/100,kgd=1/1000,kd=100,g50d=0.5,nd=20,ka=1,g50a=1,na=20)
  parmin <- c(kgr=1/10000,kgd=1/10000000,kd=0.1,g50d=0.001,nd=1,ka=0.4,g50a=0.001,na=1)
  parmax <- log(parmax)
  parmin <- log(parmin)
  par_range <- parmax-parmin
  
  init_pars <- randomLHS(Ntrial,length(parmax))
  init_pars <- t(apply(init_pars,1, function(ppi){
    ppi*par_range + parmin
  }))
  
  init_pars <- cbind(g0,n0,init_pars)
}

trial_pars <- function(init_pars,yxn){
  errs <- apply(init_pars,1,wrap_fit_mod,yxn=yxn)
}
 
wrap_full_fit <- function(Ntrial,yxn){
  init_pars <- gen_pars(yxn,Ntrial)
  errs <- trial_pars(init_pars,yxn)
  
  pars <- init_pars[errs==min(errs),]
  opt <- optim(pars,fn=wrap_fit_mod,yxn=yxn)
  res <- list(pars=opt$par,
              init=pars,
              value = opt$value,
              init_value=min(errs))
  return(res)
}

library(parallel)
ncores <- 60
Ntrial <- 100
Nstarts <- 200
setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/glucose_wgd/")
y <- readRDS("data/fitting/fitting_data.Rds")

cl <- makeCluster(getOption("cl.cores", ncores))
clusterCall(cl, function() {
  library(deSolve)
  library(lhs)
})
clusterExport(cl,varlist = c("trial_pars","gen_pars",
                             "wrap_fit_mod","mod_err",
                             "run_mod","mod"))



for(ploidy in c("2N","4N")){
  yxn <- y[grepl(ploidy,names(y))]
  yxn <- yxn[sapply(yxn,function(yi) !is.null(yi$gx))]
  res <- parLapplyLB(cl,X=rep(Ntrial,Nstarts),fun=wrap_full_fit,yxn=yxn)
  saveRDS(res,paste0("data/fitted_parameters/multistart_modB_",ploidy,".Rds"))
}










