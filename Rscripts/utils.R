data_subset <- function(X,metrics=c("dead","nuclei"),ploidy=2,glucose=5){
  X <- X[X$metric%in%metrics&X$ploidy%in%ploidy&X$glucose%in%glucose,]
  X
}
convert_glucose <- function(G,censor=TRUE) {
  G <- 12.33+G*14.9
  if(censor) G[G<20] <- 0
  G
}
media_refresh <- function(pars,parNames,g0,times=seq(0,144,1),refresh_times){
  names(pars) <- parNames
  pars <- exp(pars)
  y0 <- pars[names(pars)%in%c("N","D")]
  pars <- pars[!parNames%in%c("N","D")]
  y0['G']<- g0
  
  t0 <- times[times<=refresh_times[1]]
  x0   <- data.frame(run_mod(pars, t0, y0))
  
  for(i in 2:length(refresh_times)){
    ti <- times[times<refresh_times[i]]
    ti <- c(refresh_times[i-1],ti[!ti%in%t0])
    y0['N'] <- tail(x0$N,1)
    y0['D'] <- tail(x0$D,1)
    xi   <- data.frame(run_mod(pars, ti, y0))
    xi <- xi[-1,]
    x0 <- rbind(x0,xi)
    t0 <- c(t0,ti)
  }
  
  return(x0)}


load_dat <- function(data_dir = "~/projects/006_ploidyModeling/06_ploidy_glucose/glucose_wgd/data/processed/fitting_datasets/"){
  bx <- readRDS(paste0(data_dir,"cellcount_data.Rds"))
  gx <- readRDS(paste0(data_dir,"glucose_data.Rds"))
  list(bx=bx,gx=gx)
  
}

load_dat_old <- function(data_dir = "~/projects/006_ploidyModeling/06_ploidy_glucose/glucose_wgd/data/prior_dat/"){
  bx <- readRDS(paste0(data_dir,"01_aggregated_counts.Rds"))
  glu_dat <- readRDS(paste0(data_dir,"04_fitting_glucose.Rds"))
  fit_sd <- lm(sd~G,data=glu_dat)
  pred_sd <- predict(fit_sd,glu_dat)
  glu_dat$sd[glu_dat$sd==0] <- as.numeric(pred_sd)[glu_dat$sd==0]
  
  gx2 <- readRDS(paste0(data_dir,"05_updated_glucose.Rds"))
  
  gx2i <- aggregate(list(G=gx2$glucose),by=list(day=gx2$day,ploidy=gx2$ploidy,glucose=gx2$G),mean)
  gx2sd <- aggregate(list(sd=gx2$glucose),by=list(day=gx2$day,ploidy=gx2$ploidy,glucose=gx2$G),sd)
  gx2i$sd <- gx2sd$sd
  
  gx <- rbind(gx2i[gx2i$glucose%in%c(0.1,0.5,1,5),],
              glu_dat[glu_dat$glucose==25,])
  
  list(bx=bx,gx=gx,gx1=glu_dat,gx2=gx2i)
  
}

load_opt <- function(pars,parNames,ploidy,vars=c("alive","dead")){
  
  x <- run_full_model(pars,parNames=names(pars))
  colnames(x)[2:4] <- c("alive","dead","G")
  x <- x[,c("time",vars,"glucose")]
  x <- reshape2::melt(x,id.vars=c("time","glucose"))
  x$ploidy <- ploidy
  colnames(x)[3:4] <- c("metric","value")
  x
}

melt_data_for_plotting <- function(x){
  bx <- do.call(rbind,lapply(x,function(xi){
    id <- "newexp"
    if(is.null(xi$gx)){
      id <- "oldexp" 
    }
    bx <- xi$bx
    bx$id <- id
    bx
  }))
  gx <- do.call(rbind,lapply(x,function(xi){
    id <- "newexp"
    if(is.null(xi$gx)){
      id <- "oldexp" 
    }
    if(id=="oldexp") return(NULL)
    gx <- xi$gx
    gx$id <- id
    gx
  })) 
  list(bx=bx,gx=gx)
}

transform_fluor_for_plotting <- function(xg,gpar){
  qvals <- data.frame(do.call(rbind,lapply(xg$fluor,function(fi){
    (exp(log(fi)-qnorm(p=c(0.95,0.5,0.05),sd=gpar["sd"]))-gpar["b"])/gpar["a"]
    
  })))
  qvals <- qvals*xg$dilution
  
  colnames(qvals)<-c("Gmin","G","Gmax")
  xg <- cbind(xg,qvals)
  return(xg)
}

run_old_model <- function(pars,parNames,gluconc=c(0.1,0.5,1,5,25),times = seq(0,144,1),ploidy="2N"){
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
  di$time <- di$time/24
  colnames(di)[1:3] <- c("day","alive","dead")
  di$ploidy <- ploidy
  gx <- di[,c("day","ploidy","glucose","G")]
  bx <- di[,c("day","ploidy","glucose","alive","dead")]
  bx <- reshape2::melt(bx,id.vars=c("day","ploidy","glucose"))
  colnames(bx)[4:5] <- c("metric","ncells")
  list(bx=bx,gx=gx)
  }

fit_full_model <- function(pars,dat,parNames,fixpars=NULL,ploidy="2N",glu_dat=NULL,
                           theta_hat=NULL,gpar=NULL){
  print("FUNCTION DEPRECATED! HOPE THERES A GOOD REASON FOR RUNNING THIS...")
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
    ya <- x$N[x$time%in%(di$day[di$metric=="alive"]*24)]
    yd <- x$D[x$time%in%(di$day[di$metric=="dead"]*24)]
    
    xa <- di$ncells[di$metric=="alive"]
    xd <- di$ncells[di$metric=="dead"]
    
    sda <- di$sd[di$metric=="alive"]
    sdd <- di$sd[di$metric=="dead"]
    
    xx <- c(xa,xd)
    yy <- c(ya,yd)
    sdx <- c(sda,sdd)
    
    negll_errs <- c()
    if(!is.null(glu_dat)){
      if(is.null(gpar)) stop("glucose calibration parameters (gpar) required!")
      gxi <- glu_dat[glu_dat$glucose==gi&glu_dat$ploidy==ploidy,]
      logf <- log(gxi$fluor)
      yg <- x$G[x$time%in%(gxi$day*24)]
      logpred <- log(yg*gpar["a"]/gxi$dilution+gpar["b"])##
      errs <- logpred-logf
      negll_errs <- -dnorm(errs,mean=0,sd=gpar["sd"],log=T)
    }
    
    ll <- c((xx-yy)^2/sdx^2,negll_errs)
    ll[is.na(ll)] <- 10^11
    ll <- sum(ll)
    return(ll)
  }),error=function(e) return(10^11),
  warning = function(w) return(10^9))
  res <- sum(errs/2)+Vpen
  print(res)
  return(res)
}

## this gets names and ranges for initial value parameters.
## Current assumption is every glucose batch has a unique G0 value,
## and every experiment/ploidy has a different initial ncells0
get_ics <- function(dat,ploidy,LR=FALSE){
  ncells1 <- mean(dat$bx1$ncells[dat$bx1$ploidy==ploidy&dat$bx1$day==min(dat$bx1$day)&dat$bx1$label=="alive"])
  ncells2 <- mean(dat$bx2$ncells[dat$bx2$ploidy==ploidy&dat$bx2$day==min(dat$bx2$day)&dat$bx2$label=="alive"])
  gconc1 <- as.numeric(unique(dat$gx1$glucose))
  gconc2 <- as.numeric(unique(dat$gx2$glucose))
  
  parnames <- c(paste0("ic_",c("N1","N2")),
                paste0("ic_",c(paste0("G1_",gsub("[.]","p",gconc1)),
                               paste0("G2_",gsub("[.]","p",gconc2)))))
  
  parvals <- c(ncells1,ncells2,gconc1,gconc2)
  upper <- parvals*1.3
  upper[upper==0] <- 0.01
  lower <- parvals*0.7
  lower[lower==0] <- 0.000001
  if(LR){
    parnames <- c("ic_R1","ic_R2",parnames)
    lower <- c(0.000001,0.000001,lower)
    upper <- c(5,5,upper)
  }
  list(parnames=parnames,upper=upper,lower=lower)
}

unpack_pars <- function(pars){
  icpars <- pars[grepl("ic_",names(pars))]
  pars <- pars[!grepl("ic_",names(pars))]

  exp1 <- data.frame(cbind(N=icpars[grepl("_N1",names(icpars))],G=icpars[grepl("G1_",names(icpars))]))
  exp1$glucose <- sapply(rownames(exp1),function(ri){
    gsub("p",".",tail(unlist(strsplit(ri,split="_")),1))
  })
  if(sum(grepl("_R1",names(icpars)))){
    exp1$R <- icpars[grepl("_R1",names(icpars))]
  }
  exp2 <- data.frame(cbind(N=icpars[grepl("_N2",names(icpars))],G=icpars[grepl("G2_",names(icpars))]))
  exp2$glucose <- sapply(rownames(exp2),function(ri){
    gsub("p",".",tail(unlist(strsplit(ri,split="_")),1))
  })
  if(sum(grepl("_R2",names(icpars)))){
    exp2$R <- icpars[grepl("_R2",names(icpars))]
  }
  list(pars=pars,exp1=exp1,exp2=exp2)
}

err_experiment <- function(par,gs,ics,bx,gx,verbose=FALSE){
  
  errs <- tryCatch(sapply(1:nrow(ics),function(i){
    gi <- ics$glucose[i]
    bxi <- bx[bx$glucose==gi,]
    gxi <- gx[gx$glucose==gi,]
    y0 <- c(N=ics$N[i],D=0,G=ics$G[i])
    if("R"%in%colnames(ics)) y0 <- c(y0,R=ics$R[i])
    times <- seq(0,144,1)
  
    x   <- data.frame(run_mod(par, times, y0))
    colnames(x)[1:4] <- c("time","N","D","G")
    ya <- x$N[x$time%in%(bxi$day[bxi$label=="alive"]*24)]
    yd <- x$D[x$time%in%(bxi$day[bxi$label=="dead"]*24)]
    
    xa <- bxi$ncells[bxi$label=="alive"]
    xd <- bxi$ncells[bxi$label=="dead"]
    
    sda <- bxi$sd[bxi$label=="alive"]
    sdd <- bxi$sd[bxi$label=="dead"]
    
    xx <- c(xa,xd)
    yy <- c(ya,yd)
    sdx <- c(sda,sdd)
    
    bx_errs <- (xx-yy)^2/sdx^2
    gx_errs <- 0
    
    if("G"%in%colnames(gxi)){ ## fitting experiment 1
      yg <- x$G[x$time%in%(gxi$day*24)]
      gx_errs <- (yg-gxi$G)^2/gxi$sd^2
    }
    if("lum"%in%colnames(gxi)){##fitting to experiment 2
      logf <- log(gxi$lum)
      yg <- x$G[x$time%in%(gxi$day*24)]
      gx_errs <- -gs$pg(glu = yg/gxi$dilution,lum = gxi$lum,gpar = gs$gpar,aslog=T)
      
    }
    
    ll <- c(bx_errs,gx_errs)
    ll[is.na(ll)] <- 10^11
    ll <- sum(ll)
    if(verbose){
      print(paste(ics[i,],collapse=" "))
      print(paste("cell error:",sum(bx_errs)))
      print(paste("glucose error:",sum(gx_errs)))
    }
    return(ll)
  }),error=function(e) return(10^11),
  warning = function(w) return(10^9))
  
  res <- sum(errs)
  return(res)
}

run_experiment <- function(par,gs,ics,bx,gx){
  
  x <- tryCatch(lapply(1:nrow(ics),function(i){
    gi <- ics$glucose[i]
    bxi <- bx[bx$glucose==gi,]
    gxi <- gx[gx$glucose==gi,]
    y0 <- c(N=ics$N[i],D=0,G=ics$G[i])
    if("R"%in%colnames(ics)) y0 <- c(y0,R=ics$R[i])
    times <- seq(0,144,1)
    
    x   <- data.frame(run_mod(par, times, y0))
    colnames(x)[1:4] <- c("day","N","D","G")
    x$day <- x$day/24
    
    bx_alive <- x[,c("day","N")]
    bx_dead <- x[,c("day","D")]
    bx_alive$label <- "alive"
    bx_dead$label <- "dead"
    colnames(bx_alive)[2] <- "ncells"
    colnames(bx_dead)[2] <- "ncells"
    bxi <- rbind(bx_alive,bx_dead)
    bxi$glucose <- gi
    
    gxi <- x[,c("day","G")]
    gxi$glucose <- gi
    
    return(list(gx=gxi,bx=bxi))
    
    
  }),error=function(e) return(NULL))
  
  bx <- do.call(rbind,lapply(x,function(xi) xi$bx))
  gx <- do.call(rbind,lapply(x,function(xi) xi$gx))
  
  return(list(bx=bx,gx=gx))
}

masterfit <- function(pars,dat,parNames,gs,ploidy="2N",verbose=F){
  
  pars <- exp(pars)
  names(pars) <- parNames
  pars <- unpack_pars(pars)
  err1 <- err_experiment(par = pars$pars,gs,ics=pars$exp1,
                         bx=dat$bx1[dat$bx1$ploidy==ploidy,],
                         gx=dat$gx1[dat$gx1$ploidy==ploidy,],verbose=verbose)
  err2 <- err_experiment(par=pars$pars,gs,ics=pars$exp2,
                         bx=dat$bx2[dat$bx2$ploidy==ploidy,],
                         gx=dat$gx2[dat$gx2$ploidy==ploidy,],verbose=verbose)
  
  print(paste0("exp1 err: ",err1))
  print(paste0("exp2 err: ",err2))
  return(err1+err2)
}

masterrun <- function(pars,dat,parNames,gs,ploidy="2N"){
  
  pars <- exp(pars)
  names(pars) <- parNames
  pars <- unpack_pars(pars)
  x1 <- run_experiment(par = pars$pars,gs,ics=pars$exp1,
                         bx=dat$bx1[dat$bx1$ploidy==ploidy,],
                         gx=dat$gx1[dat$gx1$ploidy==ploidy,])
  x1$bx$id <- "exp 1"
  x1$gx$id <- "exp 1"
  x2 <- run_experiment(par=pars$pars,gs,ics=pars$exp2,
                         bx=dat$bx2[dat$bx2$ploidy==ploidy,],
                         gx=dat$gx2[dat$gx2$ploidy==ploidy,])
  x2$bx$id <- "exp 2"
  x2$gx$id <- "exp 2"
  bx<-rbind(x1$bx,x2$bx)
  gx<-rbind(x1$gx,x2$gx)
  bx$ploidy <- ploidy
  gx$ploidy <- ploidy
  cx <- plot_curves(log(pars$pars),ploidy = ploidy)
  x <- list(bx=bx,gx=gx,cx=cx)
  return(x)
}