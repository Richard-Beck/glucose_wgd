lhs_guess <- function(n,upper,lower){
  ## generate a random parameter set according to latin hypercube sampling
  scale <- upper-lower
  p0 <- randomLHS(n,length(upper))
  p0 <- do.call(cbind,lapply(1:ncol(p0), function(i) p0[,i]*scale[i]+lower[i]))
  p0
}

fit_full_model <- function(pars,dat,parNames,fixpars=NULL,ploidy="2N",glu_dat=NULL,
                           theta_hat=NULL){
  
  Vpen <- 0
  if(!is.null(theta_hat)){
    Vpen <- (sqrt(sum((pars-theta_hat)^2))-1)^2
  }
  
  dat <- dat[dat$ploidy==ploidy,]
  pars <- exp(pars)
  names(pars) <- parNames
  

  if(!is.null(fixpars)){
    pars[names(fixpars)]<-exp(fixpars)
  }
  
  delta_G <- 0
  if("delta_G"%in%names(pars)){
    delta_G <- pars["delta_G"]
    pars <- pars[!names(pars)=="delta_G"]
  }
  
  y0 <- pars[c("N","D")]
  pars <- pars[!names(pars)%in%c("N","D")]

  dat <- split(dat, f = dat$glucose)
  
  errs <- tryCatch(sapply(dat,function(di){
    gi <- di$glucose[1]
    y0["G"] <- gi+delta_G
    times <- seq(0,144,1)
    x   <- data.frame(run_mod(pars, times, y0))
    colnames(x) <- c("time","N","D","G")
    ya <- x$N[x$time%in%(8+di$frame[di$metric=="alive"]*8)]
    yd <- x$D[x$time%in%(8+di$frame[di$metric=="dead"]*8)]
    
    xa <- di$ncells[di$metric=="alive"]
    xd <- di$ncells[di$metric=="dead"]
    
    sda <- di$sd[di$metric=="alive"]
    sdd <- di$sd[di$metric=="dead"]
    
    xx <- c(xa,xd)
    yy <- c(ya,yd)
    sdx <- c(sda,sdd)
    
    if(!is.null(glu_dat)){
      gxi <- glu_dat[glu_dat$glucose==gi&glu_dat$ploidy==ploidy,]
      xg <- gxi$G
      sdg <- gxi$sd
      yg <- x$G[x$time%in%(gxi$day*24)]
      xx <- c(xx,xg)
      yy <- c(yy,yg)
      sdx <- c(sdx,sdg)
    }
    
    ll <- sum((xx-yy)^2/sdx^2)
    ll[is.na(ll)] <- 10^11
    return(ll)
  }),error=function(e) return(10^11),
  warning = function(w) return(10^9))
  res <- sum(errs/2)+Vpen
  return(res)
}

get_ics <- function(bx,ploidy){
  bx <- bx[bx$ploidy==ploidy,]
  bx <- split(bx, f = bx$glucose)
  y0 <- lapply(bx,function(di){
    gi <- di$glucose[1]
    ni <- di$ncells[di$metric=="alive"&di$frame==0]
    Di <- di$ncells[di$metric=="dead"&di$frame==0]
    c(N=ni,D=Di,G=gi)
  })
  return(y0)
}

run_full_model <- function(pars,parNames,gluconc=c(0.1,0.5,1,5,25),times = seq(0,144,1)){
  pars <- exp(pars)
  names(pars) <- parNames
  
  delta_G <- 0
  if("delta_G"%in%names(pars)){
    delta_G <- pars["delta_G"]
    pars <- pars[!names(pars)=="delta_G"]
  }
  
  y0 <- pars[names(pars)%in%c("N","D")]
  pars <- pars[!names(pars)%in%c("N","D")]
  di <- do.call(rbind,lapply(gluconc,function(gi){
    y0["G"] <- gi+delta_G
    x   <- data.frame(run_mod(pars, times, y0))
    colnames(x) <- c("time","N","D","G")
    x$glucose <- gi
    x
    
  }))
  return(di)}

wrap_2N4N_fit <- function(pars,parNames,dat,switchpar=NULL){
  bx <- dat$bx
  gx <- dat$gx
  pars2N <- pars[1:length(parNames)]
  pars4N <- pars2N
  pars4N[parNames%in%switchpar] <- tail(pars,(length(pars)-length(pars2N)))
  
  ll2N <- fit_full_model(pars2N,bx,parNames,fixpars=NULL,ploidy="2N",glu_dat=gx)
  ll4N <- fit_full_model(pars4N,bx,parNames,fixpars=NULL,ploidy="4N",glu_dat=gx)
  return(ll2N+ll4N)
}

wrap_bootstrap_2N4N <- function(switchpar,pars,parnames,fit_dat){
  p0 <- c(pars,pars[parnames%in%switchpar])
  opt <- optim(p0,fn=wrap_2N4N_fit,parNames=parnames,dat=fit_dat,switchpar=switchpar)
  print(c(switchpar,opt$value))
  list(pars=unpack_2N4N_pars(opt$par,parnames,switchpar),
       switchpar=switchpar,
       value=opt$value)
}

unpack_2N4N_pars <- function(pars,parNames,switchpar=NULL){

  pars2N <- pars[1:length(parNames)]
  names(pars2N) <- parNames
  pars4N <- pars2N
  pars4N[parNames%in%switchpar] <- tail(pars,(length(pars)-length(pars2N)))
  
  list(pars2N=pars2N,pars4N=pars4N)
}


