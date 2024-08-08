run_fit <- function(theta,ploidy,gluconc=c(0.1,0.5,1,5,25)){

  run_full_model(as.numeric(theta),names(theta),gluconc = gluconc,times=seq(0,144,1))
}

viz_fits <- function(N,ploidy,gluconc=c(0.1,0.5,1,5,25),model="modelB"){
  source(paste0("Rscripts/models/",model,".R"))
  chain <- readRDS(paste0("best_pars/",model,"/",ploidy,"_chain.Rds"))[[1]]
  ll <- chain$ll
  par <- chain$par
  par <- par[ll<(min(ll)+50),]
  par <- tail(par,round(nrow(par)/2))
  par <- par[sample(1:nrow(par),N),]
  res <- do.call(rbind,lapply(1:N, function(i){
    x <- run_fit(par[i,],ploidy,gluconc)
    x$ploidy <- ploidy
    x <- reshape2::melt(x,id.vars=c("time","glucose","ploidy"))
    colnames(x)[4] <- "metric"
    x$day <- x$time/24
    x$metric <- c(N="alive",D="dead",G="G")[x$metric]
    x$rep <- i
    x
  }))
}

delta_fits <- function(N,ploidy,gluconc=c(0.1,0.5,1,5,25),model="modelB",rev="rev0"){
  source(paste0("Rscripts/models/",model,".R"))
  chain <- readRDS(paste0("best_pars/",model,"/",rev,"/",ploidy,"_chain.Rds"))[[1]]
  ll <- chain$ll
  par <- chain$par
  par <- par[ll<(min(ll)+50),]
  par <- tail(par,round(nrow(par)/2))
  par <- par[sample(1:nrow(par),N),]
  res <- lapply(1:N, function(i){
    r <- run_fit(par[i,],ploidy,gluconc)

  })
  

  N <- do.call(cbind,lapply(res, function(ri) ri$N))
  N <- data.frame(t(apply(N,1, quantile,probs=c(0.05,0.5,0.95))))
  colnames(N) <- c("lower","middle","upper")
  N$metric <- "alive"
  N$time <- res[[1]][,"time"]
  N$glucose <- res[[1]][,"glucose"]
  
  D <- do.call(cbind,lapply(res, function(ri) ri$D))
  D <- data.frame(t(apply(D,1, quantile,probs=c(0.05,0.5,0.95))))
  colnames(D) <- c("lower","middle","upper")
  D$metric <- "dead"
  D$time <- res[[1]][,"time"]
  D$glucose <- res[[1]][,"glucose"]
  
  G <- do.call(cbind,lapply(res, function(ri) ri$G))
  G <- data.frame(t(apply(G,1, quantile,probs=c(0.05,0.5,0.95))))
  colnames(G) <- c("lower","middle","upper")
  G$metric <- "G"
  G$time <- res[[1]][,"time"]
  G$glucose <- res[[1]][,"glucose"]
  
  x <- rbind(N,D,G)
  x$ploidy <- ploidy
  x$model <- model
  x$day <- x$time/24
  x <- x[!colnames(x)=="time"]
  x
  
}



