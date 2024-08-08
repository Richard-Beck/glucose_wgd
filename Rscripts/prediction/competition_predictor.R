mixexp <- function(p0,G=1,times=1:144){
  #p00 <- p0
  y0 <- list(c(N=250,D=0),c(N=250,D=0))
  
  run_mixmod(p0,times=times,y0=y0,G=G)
  
}

viz_mixexp <- function(G,N=100,times=1:144,model="modelB",rev="rev0"){
  source(paste0("Rscripts/models/",model,".R"))
  c2<- readRDS(paste0("best_pars/",model,"/",rev,"/2N_chain.Rds"))[[1]]
  c4<- readRDS(paste0("best_pars/",model,"/",rev,"/4N_chain.Rds"))[[1]]
  
  c2$par <- c2$par[c2$ll<(min(c2$ll)+50),]
  c4$par <- c4$par[c4$ll<(min(c4$ll)+50),]
  
  c4 <- tail(c4$par,round(nrow(c4$par)/2))
  c2 <- tail(c2$par,round(nrow(c2$par)/2))
  
  c2 <- c2[sample(1:nrow(c2),N),]
  c4 <- c4[sample(1:nrow(c4),N),]
  

  
  z <- do.call(rbind,lapply(1:N, function(i){
    p0 <- list(c2[i,],
               c4[i,])
    z <- mixexp(p0,G,times)
    df2 <- data.frame(z[,c("time","N1","D1")])
    df4 <- data.frame(z[,c("time","N2","D2")])
    colnames(df2) <- c("hours","alive","dead")
    colnames(df4) <- c("hours","alive","dead")
    df2$ploidy <- "2N"
    df4$ploidy <- "4N"
    df <- rbind(df2,df4)
    df <- reshape2::melt(df,id.vars=c("hours","ploidy"))
    df$rep <- i
    df
  }))
  z$G <- G
  return(z)
}

delta_mixexp <- function(G,N=100,times=1:144, model="modelB",rev="rev0"){
  source(paste0("Rscripts/models/",model,".R"))
  c2<- readRDS(paste0("best_pars/",model,"/",rev,"/2N_chain.Rds"))[[1]]
  c4<- readRDS(paste0("best_pars/",model,"/",rev,"/4N_chain.Rds"))[[1]]
  
  c2$par <- c2$par[c2$ll<(min(c2$ll)+50),]
  c4$par <- c4$par[c4$ll<(min(c4$ll)+50),]
  
  c4 <- tail(c4$par,round(nrow(c4$par)/2))
  c2 <- tail(c2$par,round(nrow(c2$par)/2))
  
  c2 <- c2[sample(1:nrow(c2),N),]
  c4 <- c4[sample(1:nrow(c4),N),]
  
  z <- do.call(cbind,lapply(1:N, function(i){
    p0 <- list(c2[i,],
               c4[i,])
    z <- mixexp(p0,G,times)
    as.numeric(z[,"N1"])-as.numeric(z[,"N2"])
  }))
  z <- data.frame(t(apply(z,1,quantile,probs=c(0.05,0.5,0.95),na.rm=T)))
  colnames(z) <- c("lower","middle","upper")
  z$time <- times
  z$G <- G
  return(z)
}


