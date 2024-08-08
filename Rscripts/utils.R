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