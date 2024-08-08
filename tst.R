get_best <- function(parlist){
  errs <- sapply(parlist,function(pp) {
    pp$value
  })
  parlist[[which(errs==min(errs))]]$par
}
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

ll_wrapper <- function(iniglu_2N.0, iniglu_2N.0.1, iniglu_2N.0.25, iniglu_2N.0.5, iniglu_2N.1,
                       inicell_2N.0, inicell_2N.0.1, inicell_2N.0.25, inicell_2N.0.5, inicell_2N.1,
                       kgr, kgd, kd, g50d, nd, ka, g50a, na) {
  pars <- c(iniglu_2N.0 = iniglu_2N.0, iniglu_2N.0.1 = iniglu_2N.0.1, iniglu_2N.0.25 = iniglu_2N.0.25,
            iniglu_2N.0.5 = iniglu_2N.0.5, iniglu_2N.1 = iniglu_2N.1, inicell_2N.0 = inicell_2N.0,
            inicell_2N.0.1 = inicell_2N.0.1, inicell_2N.0.25 = inicell_2N.0.25, inicell_2N.0.5 = inicell_2N.0.5,
            inicell_2N.1 = inicell_2N.1, kgr = kgr, kgd = kgd, kd = kd, g50d = g50d, nd = nd,
            ka = ka, g50a = g50a, na = na)
  res <- -wrap_fit_mod(pars,yxn)
  return(res)
}

y <- readRDS("data/fitting/fitting_data.Rds")
pars<-get_best(readRDS(paste0("data/fitted_parameters/multistart_modB_2N.Rds")))
yxn <- y[grepl("2N",names(y))]
yxn <- yxn[sapply(yxn,function(yi) !is.null(yi$gx))]

opt=mle2(ll_wrapper,start = as.list(pars))
as.list(pars)
