library(deSolve)
library(sundialr)



mod <- function(t, x, theta){
  dx <- vector(mode = "numeric", length = length(x))
   
  dx[1] <- theta[1]*x[1]*(1-x[1]/theta[2])*x[3]^theta[3]/(x[3]^theta[3]+theta[4]^theta[3])-theta[5]*x[1]*(theta[12]/theta[5]*x[4]/(x[4]+theta[11])+1-x[3]^theta[6]/(x[3]^theta[6]+theta[7]^theta[6]))
  dx[2] <- theta[5]*x[1]*(theta[12]/theta[5]*x[4]/(x[4]+theta[11])+1-x[3]^theta[6]/(x[3]^theta[6]+theta[7]^theta[6]))
  dx[3] <- -theta[8]*x[1]*x[3]^theta[9]/(x[3]^theta[9]+theta[10]^theta[9])
  dx[4] <- x[1]
  dx
}

mixmod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    ## hypothetical toxic compound that is being generated
    dY <- N1+N2
    ##joint carrying capacity
    jcc <- 0.5*(N1/theta1+N2/theta2)
    
    dN1 <- kp1*N1*(1-jcc)*G^na1/(G^na1+g50a1^na1)-kd1*N1*(ky1/kd1*Y/(Y+g50y1)+1-G^nd1/(G^nd1+g50d1^nd1))
    dD1 <- kd1*N1*(ky1/kd1*Y/(Y+g50y1)+1-G^nd1/(G^nd1+g50d1^nd1))
    
    dN2 <- kp2*N2*(1-jcc)*G^na2/(G^na2+g50a2^na2)-kd2*N2*(ky2/kd2*Y/(Y+g50y2)+1-G^nd2/(G^nd2+g50d2^nd2))
    dD2 <- kd2*N2*(ky2/kd2*Y/(Y+g50y2)+1-G^nd2/(G^nd2+g50d2^nd2))
    
    dG <- -v1*N1*G^m1/(G^m1+g50c1^m1)-v2*N2*G^m2/(G^m2+g50c2^m2)
    
    return(list(c(dN1,dD1,dN2,dD2,dG,dY)))
  })
}


##reparameterise with:
## x4=(G/g50a)^n
## x5=(G/g50d)^n
## x6=(G/g50c)^n

mod_trans <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    ## hypothetical toxic compound that is being generated
    dY <- N
    dNY <- ky*Y/(Y+g50Y)
    
    dN <- kp*N*(1-N/theta)*x4/(x4+1)-kd*N*(dNY/kd+1-x5/(x5+1))
    dD <- kd*N*(dNY/kd+1-x5/(x5+1))
    dG <- -v*N*x6/(x6+1)
    
    dx4 <- -v*N*x4/(x4+1)*n*x4/G
    dx5 <- -v*N*x5/(x5+1)*n*x5/G
    dx6 <- -v*N*x6/(x6+1)*m*x6/G
    
    return(list(c(dN,dD,dG,dY,dx4,dx5,x6)))
  })
}
## wrapper to run ODE model
run_mod <- function(p0,times,y0){
  theta <- c(kp=0.1,theta=10000,na=1,g50a=0.1,kd=0.1,nd=1,g50d=0.1,v=1e-6,m=1,g50c=0.1,g50y=8000,ky=0.1)
  theta[names(p0)] <- p0
  reltol <- 1e-8
  abstol <- c(1e-8,1e-8,1e-8,1e-8)
  if(!"Y"%in%names(y0)) y0 <- c(y0,c(Y=0))
  y0["G"] <- y0["G"]
    
  ## Solving the ODEs using cvode function
  cvode(times, y0, mod , theta, reltol, abstol) 
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
  
  base_y0 <- c(N1=1000,D1=0,N2=1000,D2=0,G=G,Y=0)
  base_pars_02<-c(kp=0.1,theta=10000,na=1,g50a=0.1,kd=0.1,nd=1,g50d=0.1,v=1e-6,m=1,g50c=0.1,g50y=8000,ky=0.1)
  names(base_pars_02) <- paste0(names(base_pars_02),"1")
  base_pars_04<-c(kp=0.1,theta=10000,na=1,g50a=0.1,kd=0.1,nd=1,g50d=0.1,v=1e-6,m=1,g50c=0.1,g50y=8000,ky=0.1)
  names(base_pars_04) <- paste0(names(base_pars_04),"2")
  base_pars <- c(base_pars_02,base_pars_04)
  base_pars[names(p0)] <- p0
  base_y0[names(y0)] <- y0
  ode(base_y0, times, mixmod, base_pars)
}

##info required for running or fitting model
model_info <- function(){
  parnames <- c("kp", "theta", "na", "g50a", "kd", "nd", "g50d", "v"  , "m", "g50c", "g50y", "ky" )
  lower <-    c(0.01, 6000   , 0.1 , 0.001 , 0.01, 0.1 , 0.001 , 1e-8 , 0.1, 0.001 , 5000  , 0.001)
  upper <-    c(1   , 50000  , 25  , 5     , 1   , 25  , 5     , 0.001, 25 , 5     , 500000, 1)
  list(parnames=parnames,upper=upper,lower=lower)
}

