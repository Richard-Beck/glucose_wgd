library(deSolve)
library(sundialr)

#literally identical to modelA except there is an additional parameter delta_G
## which adds a constant amount of glucose to each well. 


mod <- function(t, x, theta){
  dx <- vector(mode = "numeric", length = length(x))
  XL <- 0.5*(1+tanh(theta[12]*(x[4]-theta[13])))
  XLstar <- 1-XL
  ## the 25 in XG is chosen arbitrarily (for now). The OG authors
  ## used a switch function but I would rather use a smooth function 
  ## and just make it sufficiently steep.
  XG <- (0.5+tanh(25*(x[3]-theta[14]))/2) 
   
  dx[1] <- x[1]/theta[1]*(1-x[1]-x[2])+x[2]*XL/theta[2]-XLstar*XG*x[1]/theta[3]
  dx[2] <- x[2]/theta[4]*(1-x[1]-x[2])-x[2]*XL/theta[2]+XLstar*XG*x[1]/theta[3]
  dx[3] <- -theta[5]*x[3]*x[1]/(x[3]+theta[6]*x[4]+theta[7])-theta[8]*x[3]*x[2]/(x[3]+theta[9])
  dx[4] <- -theta[10]*x[4]/(x[3]/theta[6]+x[4]+theta[11]) +2*theta[8]*x[3]*x[2]/(x[3]+theta[9])
  dx
}

## wrapper to run ODE model
run_mod <- function(p0=NULL,times=seq(0,20,0.1),y0=c(0.01,0,1,0)){
  
  p0 <- c(to=1/log(2),tg=1/log(2),pstar=850*10^3,tog=1/24,tgo=1,gamma=100,
          lstar=10,ko=1.5,kl=5,kg=10,lambda=1,nstar=1,mstar=1,gstar=0.5,gmin=0.5)
  theta <- c(to=1,tgo=1,tog=1,tg=1, ko=1, lambda=1,nstar=1, kg = 1, gstar=1,kl=1,mstar=1,gamma=1,lstar=1,gmin=1)
  pstar <- p0["pstar"]
  theta <- p0[names(theta)]
  reltol <- 1e-8
  abstol <- c(1e-8,1e-8,1e-8,1e-8)

  ## Solving the ODEs using cvode function
  out <- data.frame(cvode(times, y0, mod , theta, reltol, abstol))
  colnames(out) <- c("time","Po","Pg","G","L")
  out$Po <- out$Po*pstar
  out$Pg <- out$Pg*pstar
  out
}

out <- run_mod()
colnames(out) <- c("time","P0","Pg","G","L")
out <- reshape2::melt(out,id.vars="time")
library(ggplot2)
p <- ggplot(out,aes(x=time,y=value))+
  facet_wrap(~variable,scales="free")+
  geom_line()
p



