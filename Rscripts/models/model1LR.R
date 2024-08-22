#as model 1 but now decouple birth and death.
library(deSolve)

mod <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    Sd <- (1-(R+G)/((R+G)+Gmin_d))*tanh(Time*kp)
    Sp <- (R+G)/((R+G)+Gmin_p)*tanh(Time*kp)
    
    dG <- -v*N*G/(G+Gstar)
    dL <- -2*dG ## lactate
    dR <- -vr*N*R/(R+Rstar) ## glucose
    
    dN <- kp*N*(1-N/theta)*Sp-kd*N*Sd-kbys*N*D/(D+N)-kdl*L/(L+L50)
    dD <- kd*N*Sd+kbys*N*D/(D+N)+kdl*L/(L+L50)
    
    return(list(c(dN,dD,dG,dR,dL)))
  })
}

## wrapper to run ODE model
run_mod <- function(p0,times,y0){
  base_y0 <- c(N=1000,D=0,G=5,R=10,L=0)
  base_pars<-c(kp=0.1,kd=0.1,kbys=0.1,theta=15000,v=0.00001,Gmin_p=0.1,gmin_d=0.01,Gstar=0.1,
               vr=0.00001,Rstar=0.1,kdl=0.1,L50=25)
  base_pars[names(p0)] <- p0
  base_y0[names(y0)] <- y0
  ode(base_y0, times, mod, base_pars)
}

##info required for running or fitting model

model_info <- function(){
  parnames <- c("theta", "kp", "kd", "kbys", "v"  , "Gmin_p","Gmin_d", "Gstar","vr" ,"Rstar","kdl","L50")
  lower <-    c(6000   , 0.01, 0.01, 0.001 , 1e-8 , 0.01    ,0.01    , 0.01   ,1e-8 ,0.01   ,0.001,5)
  upper <-    c(50000  , 1   , 5   , 1     , 0.001, 20      ,20      , 20     ,0.001,20     ,1    ,50)
  list(parnames=parnames,upper=upper,lower=lower)
}

plot_curves <- function(pars,ploidy,G=seq(0.001,1,0.001)){
  pars <- exp(pars)
  with(as.list(pars),{
    df <- data.frame(G,ploidy,
                     div=kp*G/(G+Gmin_p),
                     death=kd*(1-G/(G+Gmin_d)),
                     cons=v*G/(G+Gstar))
    df$ploidy <- ploidy
    df <- reshape2::melt(df,id.vars=c("G","ploidy"))
  })  
}

