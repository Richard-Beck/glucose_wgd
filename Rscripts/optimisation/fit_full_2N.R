mod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    #if (Time > delay) {
    # G_delayed <- lagvalue(Time - delay, 3)
    #} else {
    # G_delayed <- 25  # If the time is less than the delay, use the initial value
    #}
    dN <- ka*N*(1-N/theta)/(1+(g50a/G)^na)-kd*N*(1-1/(1+(g50d/G)^nd))
    dD <- kd*N*(1-1/(1+(g50d/G)^nd))
    dG <- -kg*N*G^ng/(G^ng+g50c^ng)
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
  N0 <- exp(pars[grepl("init_cells",names(pars))])
  G0 <- exp(pars[grepl("init_G",names(pars))])
  
  pars <- pars[!grepl("init_cells",names(pars))]
  pars <- pars[!grepl("init_G",names(pars))]
  
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

library(deSolve)
setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/glucose_wgd/")
y <- readRDS("data/fitting/fitting_data.Rds")
y2n <- y[grepl("2N",names(y))]
y2n <- y2n[sapply(y2n,function(yi) !is.null(yi$gx))]
pars <- readRDS("data/fitted_parameters/init_2N.Rds")
pars["theta"] <- 30000

N0 <- sapply(y2n,function(yi) yi$bx$ncells[yi$bx$metric=="alive"&yi$bx$day==min(yi$bx$day)])
G0 <- sapply(y2n,function(yi) yi$gx$G[1])
names(N0) <- paste0("init_cells",1:length(N0))
names(G0) <- paste0("init_G",1:length(G0))
pars <- c(N0,G0,pars)

#pars["delay"] <- 1
#pars["N0"] <- 600
pars <- log(abs(pars))
opt2N <- optim(pars,fn=wrap_fit_mod,yxn=y2n)
saveRDS(opt2N,"data/fitted_parameters/final_2N_delay.Rds")
