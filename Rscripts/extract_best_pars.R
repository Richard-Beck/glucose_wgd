## this script doesn't really do anything. Just loads different fit outputs
## and prints the comments so they can be copy pasted into the right subfolder of 
## best_pars
setwd("C:/Users/4473331/Documents/projects/006_ploidyModeling/06_ploidy_glucose/00_modelling_v2/")

model <- "modelB"
choice <- 20

ff <- list.files(paste0("opt_out/",model))
print(ff)
print(ff[choice])
x <- readRDS(paste0("opt_out/",model,"/",ff[choice]))

print(x$comments)

