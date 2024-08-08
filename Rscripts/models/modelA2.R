library(deSolve)
library(sundialr)

mod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    ## hypothetical toxic compound that is being generated
    dY <- N
    dGi <- v*G^m/(G^m+g50c^m)-Gi*(kgi+kp*(1-N/theta)*(Gi)^na/((Gi)^na+g50a^na))
    dN <- kp*N*(1-N/theta)*(Gi)^na/((Gi)^na+g50a^na)-kd*N*(ky/kd*Y/(Y+g50y)+1-(Gi)^nd/((Gi)^nd+g50d^nd))
    dD <- kd*N*(ky/kd*Y/(Y+g50y)+1-(Gi)^nd/((Gi)^nd+g50d^nd))
    dG <- -v*N*G^m/(G^m+g50c^m)
    return(list(c(dN,dD,dG,dY,dGi)))
  })
}
##steady state internal glucose concentration
get_Gi <- function(G=5,pars){
  with(as.list(pars),{
    return(v*G^m/(G^m+g50c^m)/(kgi+kp))
  })
}
## wrapper to run ODE model
run_mod <- function(p0=NULL,times=0:144,y0=NULL){
  base_y0 <- c(N=1000,D=0,G=1,Y=0,Gi=0)
  base_pars<-c(kp=0.05,theta=15000,na=5,g50a=0.1,kd=0.1,nd=5,g50d=0.05,v=1e-4,m=1,g50c=0.1,g50y=8000,ky=0,kgi=0.01)
  base_pars[names(p0)] <- p0
  base_pars["g50a"] <- get_Gi(base_pars["g50a"],base_pars)
  base_pars["g50d"] <- get_Gi(base_pars["g50d"],base_pars)
  base_y0[names(y0)] <- y0
  base_y0["Gi"] <- get_Gi(5,base_pars)
  ode(base_y0, times, mod, base_pars)
}




##reparameterise with:
## x4=(Gi/g50a)^na
## x5=(Gi/g50d)^nd
## x6=(G/g50c)^n

mod_trans <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    ## hypothetical toxic compound that is being generated
    dY <- N
    dGi <- v*G^m/(G^m+g50c^m)-Gi*(kgi+kp*(1-N/theta)*(Gi)^na/((Gi)^na+g50a^na))
    dN <- kp*N*(1-N/theta)*Gi^na/(Gi^na+g50a^na)-kd*N*(ky/kd*Y/(Y+g50y)+1-Gi^nd/(Gi^nd+g50d^nd))
    dD <- kd*N*(ky/kd*Y/(Y+g50y)+1-Gi^nd/(Gi^nd+g50d^nd))
    dG <- -v*N*G^m/(G^m+g50c^m)
    return(list(c(dN,dD,dG,dY,dGi)))
  })
}

mod <- function(t, x, theta){
  dx <- vector(mode = "numeric", length = length(x))
  
  dx[1] <- theta[1]*x[1]*(1-x[1]/theta[2])*x[5]^theta[3]/(x[5]^theta[3]+theta[4]^theta[3])-theta[5]*x[1]*(theta[12]/theta[5]*x[4]/(x[4]+theta[11])+1-x[5]^theta[6]/(x[5]^theta[6]+theta[7]^theta[6]))
  dx[2] <- theta[5]*x[1]*(theta[12]/theta[5]*x[4]/(x[4]+theta[11])+1-x[5]^theta[6]/(x[5]^theta[6]+theta[7]^theta[6]))
  dx[3] <- -theta[8]*x[1]*x[3]^theta[9]/(x[3]^theta[9]+theta[10]^theta[9])
  dx[4] <- x[1]
  dx[5] <- theta[8]*x[3]^theta[9]/(x[3]^theta[9]+theta[10]^theta[9])-x[5]*(theta[13]+theta[1]*(1-x[1]/theta[2])*x[5]^theta[3]/(x[5]^theta[3]+theta[4]^theta[3]))
  dx
}

## wrapper to run ODE model
run_mod <- function(p0,times=0:144,y0){
  base_y0 <- c(N=1000,D=0,G=1,Y=0,Gi=0)
  base_y0[names(y0)] <- y0
  theta <- c(kp=0.05,theta=15000,na=5,g50a=0.1,kd=0.1,nd=5,g50d=0.05,v=1e-4,m=1,g50c=0.1,g50y=8000,ky=0,kgi=0.01)
  theta[names(p0)] <- p0
  theta["g50a"] <- get_Gi(theta["g50a"],theta)
  theta["g50d"] <- get_Gi(theta["g50d"],theta)
  reltol <- 1e-8
  abstol <- c(1e-8,1e-8,1e-8,1e-8,1e-8)
  base_y0["Gi"] <- get_Gi(5,theta)
  
  ## Solving the ODEs using cvode function
  cvode(times, base_y0, mod , theta, reltol, abstol) 
}







##info required for running or fitting model
model_info <- function(){
  parnames <- c("kp", "theta", "na", "g50a", "kd", "nd", "g50d", "v"  , "m", "g50c", "g50y", "ky" , "kgi","N" , "D")
  lower <-    c(0.01, 6000   , 0.1 , 0.001 , 0.01, 0.1 , 0.001 , 1e-8 , 0.1, 0.001 , 5000  , 0.001, 0.01 ,100 , 1)
  upper <-    c(1   , 50000  , 150  , 5     , 1   , 150  , 5     , 0.001, 25 , 5     , 500000, 1    , 10   ,1000, 500)
  list(parnames=parnames,upper=upper,lower=lower)
}





