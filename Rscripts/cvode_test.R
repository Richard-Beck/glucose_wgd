setwd("~/projects/006_ploidyModeling/06_ploidy_glucose/00_modelling_v2/")



modx <- function(t, x, theta){
  dx <- vector(mode = "numeric", length = length(x))
    
  dx[1] <-  theta[1]*x[1]*(1-x[1]/theta[2])*1/(1+(theta[3]/x[3])^theta[4])-theta[5]*x[1]*(1-1/(1+(theta[6]/x[3])^theta[7]))-theta[8]*x[1]^2/theta[2]
  dx[2] <- theta[5]*x[1]*(1-1/(1+(theta[6]/x[3])^theta[7]))+theta[8]*x[1]^2/theta[2]
  dx[3] <- -x[1]*(theta[9]/(1+(theta[3]/x[3])^theta[4])+theta[10]/(1+(theta[6]/x[3])^theta[7]))/2
  
  dx
}

par <- readRDS("best_pars/modelB/4N.Rds")
par <- exp(par$par[c("kp","theta","g50a","na","kd","g50d","nd","kd2","v1","v2")])
y0 <- c(N=500,D=50,G=5)

times <- seq(0,144,0.01)
reltol <- 1e-04
abstol <- c(1e-8,1e-8,1e-6)

## Solving the ODEs using cvode function
df1 <- cvode(times, y0, modx , par, reltol, abstol) 

source("Rscripts/models/modelB.R")

df2 <- ode(y0, times, mod, par)

system.time(cvode(times, y0, modx , par, reltol, abstol) )
system.time(ode(y0, times, mod, par,method="ode23"))
