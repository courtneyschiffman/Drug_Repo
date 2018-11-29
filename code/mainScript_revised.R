rm(list=ls())

# setwd("") ## set working directory to git repo

# check installed packages
libs <- c("idr",
          "parallel",
          "snow")
notInstalled <- libs[!libs %in% rownames(installed.packages())]
if(length(notInstalled)>0) {install.packages(notInstalled)}
require(snow)

nSim <- 100 # number of simulations to perform
nPerm <- 1000 # number of permutations to perform
deg <- 49 # degrees of freedom to simulate from
nCores <- parallel::detectCores() # number of clusters to use for makeCluster()

## source script files
scriptFiles <- grep(".R", list.files(path="./code"), value=T)
scriptFiles <- grep("main", scriptFiles, invert=T, value=T)
sapply(scriptFiles, function(x) source(paste0("./code/",x)))

## load simulated data
params <- expand.grid(p1=c(0.05,0.20),
                      p2=0.05,
                      p3=c(0,0.25,0.5,0.75,1))
params_names <- apply(params, 1, function(x) paste("p1",x[1],"p2",x[2],"p3",x[3], sep="_"))
# load(file="./data/simulations.Rda")
set.seed(1)
dat <- lapply(1:nrow(params),
              function(i) {
                print(paste0("simulating ", nSim, " simulations of ", params_names[i]))
                paramSims <- lapply(1:nSim, function(x) {
                  simulateData(p1=params[i,1], p2=params[i,2], p3=params[i,3])
                })
                names(paramSims) <- paste(params_names[i], "sim", 1:nSim, sep="_")
                return(paramSims)
              })
names(dat) <- params_names

## assess significance
simFuns <- c("gsea", "IDR.func", "spearman", "trunSpearman")
simFunsDir <- c(F, T, F, F)
# stupidity check
system.time({assessSig(dat[[1]][[1]], simFuns=simFuns, simFunsDir=simFunsDir, nPerm=5)})
# initiate cluster
cl <- makeCluster(getOption("cl.cores", nCores))
clusterExport(cl, list=ls())
for(i in 1:length(params_names)) {
  print(paste("assessing for", params_names[i]))
  system.time({assign(params_names[i],
                      value=parLapply(cl, dat[[i]],
                                      function(paramSim) {
                                        assessSig(paramSim, simFuns=simFuns, simFunsDir=simFunsDir, nPerm=nPerm)
                                      }))})
  save(list=c(params_names[i]), file=paste0("./data/", params_names[i], ".Rda"))
  wflow_publish(files=paste0("./data/", params_names[i], ".Rda"))
  wflow_git_push()
}
stopCluster(cl)
