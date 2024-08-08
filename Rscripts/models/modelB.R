#constrained glucose consumption and confluency cell death
library(deSolve)
library(sundialr)

mod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dN <- kp*N*(1-N/theta)*1/(1+(g50a/G)^na)-kd*N*(1-1/(1+(g50d/G)^nd))-kd2*N^2/theta
    dD <- kd*N*(1-1/(1+(g50d/G)^nd))+kd2*N^2/theta
    dG <- -N*(v1*(1/(1+g50a^na/G^na))+v2*(1/(1+g50d^nd/G^nd)))/2
    return(list(c(dN,dD,dG)))
  })
}

mod <- function(t, x, theta){
  dx <- vector(mode = "numeric", length = length(x))
  x[x<0] <- 0
  dx[1] <-  theta[1]*x[1]*(1-x[1]/theta[2])*1/(1+(theta[3]/x[3])^theta[4])-theta[5]*x[1]*(1-1/(1+(theta[6]/x[3])^theta[7]))-theta[8]*x[1]^2/theta[2]
  dx[2] <- theta[5]*x[1]*(1-1/(1+(theta[6]/x[3])^theta[7]))+theta[8]*x[1]^2/theta[2]
  dx[3] <- -x[1]*(theta[9]/(1+(theta[3]/x[3])^theta[4])+theta[10]/(1+(theta[6]/x[3])^theta[7]))/2
  
  dx
}

mixmod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    ##joint carrying capacity
    jcc <- 0.5*(N1/theta1+N2/theta2)
    
    dN1 <- kp1*N1*(1-jcc)*1/(1+(g50a1/G)^na1)-kd1*N1*(1-1/(1+(g50d1/G)^nd1))-kd21*N1*jcc
    dD1 <- kd1*N1*(1-1/(1+(g50d1/G)^nd1))+kd21*N1*jcc
    dG1 <- -N1*(v11*(1/(1+g50a1^na1/G^na1))+v21*(1/(1+g50d1^nd1/G^nd1)))/2
    
    dN2 <- kp2*N2*(1-jcc)*1/(1+(g50a2/G)^na2)-kd2*N2*(1-1/(1+(g50d2/G)^nd2))-kd22*N2*jcc
    dD2 <- kd2*N2*(1-1/(1+(g50d2/G)^nd2))+kd22*N2*jcc
    dG2 <- -N2*(v12*(1/(1+g50a2^na2/G^na2))+v22*(1/(1+g50d2^nd2/G^nd2)))/2
    
    dG <- dG1+dG2
    return(list(c(dN1,dD1,dN2,dD2,dG)))
  })
}

mod_trans <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    ## define x5 <- G^na and x6 <- G^nd
    dN <- kp*N*(1-N/theta)*x4/(x4+1)-kd*N*(1-x5/(x5+1))-kd2*N^2/theta
    dD <- kd*N*(1-x5/(x5+1))+kd2*N^2/theta
    dG <- -N*(v1*x4/(x4+1)+v2*x5/(x5+1))/2
    dx4 <- -N*(v1*x4/(x4+1)+v2*x5/(x5+1))/2*na*x4/G
    dx5 <- -N*(v1*x4/(x4+1)+v2*x5/(x5+1))/2*nd*x5/G

    
    return(list(c(dN,dD,dG,dx4,dx5)))
  })
}

## wrapper to run ODE model
run_mod <- function(p0,times,y0){
  base_y0 <- c(N=1000,D=0,G=5)
  base_pars<-c(kp=0.1,kd=0.1,kd2=0.01,v1=0.00001,v2=0.00001,theta=10000,g50a=0.2,na=1,g50d=0.01,nd=1)
  base_pars[names(p0)] <- p0
  base_y0[names(y0)] <- y0
  ode(base_y0, times, mod, base_pars)
}

## wrapper to run ODE model
run_mod <- function(p0,times,y0){
  theta <- c(kp=0.1,kd=0.1,kd2=0.01,v1=0.00001,v2=0.00001,theta=10000,g50a=0.2,na=1,g50d=0.01,nd=1)
  theta <- theta[c("kp","theta","g50a","na","kd","g50d","nd","kd2","v1","v2")]
  theta[names(p0)] <- p0
  #base_y0[names(y0)] <- y0
  reltol <- 1e-04
  abstol <- c(1e-8,1e-8,1e-6)
  
  ## Solving the ODEs using cvode function
  cvode(times, y0, mod , theta, reltol, abstol) 
  #ode(y0, times, mod, theta)
}

run_mixmod <- function(p0,times,y0,G){
  
  p02 <- unlist(p0[[1]])
  p04 <- unlist(p0[[2]])
  p02 <- p02[!names(p02)%in%c("N","D")]
  p04 <- p04[!names(p04)%in%c("N","D")]
  names(p02) <- paste0(names(p02),"1")
  names(p04) <- paste0(names(p04),"2")
  p0 <- exp(c(p02,p04))
  
  y02 <- y0[[1]]
  y04 <- y0[[2]]
  y02 <- y02[c("N","D")]
  y04 <- y04[c("N","D")]
  names(y02) <- paste0(names(y02),"1")
  names(y04) <- paste0(names(y04),"2")
  y0 <- c(y02,y04)
  
  base_y0 <- c(N1=1000,D1=0,N2=1000,D2=0,G=G)
  base_pars<-c(kp1=0.1,kd1=0.1,kd21=0.01,v11=0.00001,v21=0.00001,theta1=10000,g50a1=0.2,na1=1,g50d1=0.01,nd1=1,
               kp2=0.1,kd2=0.1,kd22=0.01,v12=0.00001,v22=0.00001,theta2=10000,g50a2=0.2,na2=1,g50d2=0.01,nd2=1)
  base_pars[names(p0)] <- p0
  base_y0[names(y0)] <- y0
  ode(base_y0, times, mixmod, base_pars)
}

## wrapper to run transformed ODE model
run_mod_trans <- function(p0,times,y0){
  base_y0 <- c(N=1000,D=0,G=5,x4=5,x5=5)
  base_pars<-c(kp=0.1,kd=0.1,kd2=0.01,v1=0.00001,v2=0.00001,theta=10000,g50a=0.2,na=1,g50d=0.01,nd=1)
  base_pars[names(p0)] <- p0
  base_y0[names(y0)] <- y0
  
  base_y0["x4"] <- (base_y0["G"]/base_pars["g50a"])^base_pars["na"]
  base_y0["x5"] <- (base_y0["G"]/base_pars["g50d"])^base_pars["nd"]
  
  ode(base_y0, times, mod_trans, base_pars)
}

##info required for running or fitting model

plot_curves <- function(pars,ploidy){
  pars <- exp(pars)
  with(as.list(pars),{
    G <- seq(0.001,1,0.001)
    df <- data.frame(G,ploidy,
                     div=kp/(1+(g50a/G)^na),
                     death=kd*(1-1/(1+(g50d/G)^nd)),
                     cons=v1*(1/(1+(g50a/G)^na))+v2*(1/(1+g50d/G)^nd))
    df$ploidy <- ploidy
    df <- reshape2::melt(df,id.vars=c("G","ploidy"))
  })  
}

model_info <- function(){
  parnames <- c("theta","kp","kd","kd2","g50a","g50d","na","nd","v1"  ,"v2" ,"N" ,"D")
  lower <-    c(6000   ,0.01,0.01,0.001, 0.001, 0.001, 0.1,0.1 , 1e-8 ,1e-8 ,100 ,1)
  upper <-    c(50000  ,1   ,1   ,    1, 5    , 1    , 25 ,25  , 0.001,0.001,1000,500)
  list(parnames=parnames,upper=upper,lower=lower)
}


