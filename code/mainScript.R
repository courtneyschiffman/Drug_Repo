rm(list=ls())

# check installed packages
libs <- c("idr",
          "parallel",
          "doParallel",
          "foreach",
          "iterators")
notInstalled <- libs[!libs %in% rownames(installed.packages())]
if(length(notInstalled)>0) {install.packages(notInstalled)}
library(doParallel)
library(foreach)
library(iterators)

nSim <- 5 # number of simulations to perform
nPerm <- 1000 # number of permutations to perform
deg <- 49 # degrees of freedom to simulate from
nCores <- parallel::detectCores() # number of clusters to use for makeCluster()

## source script files
helperFuns <- c("assessSig",
                 "GSEA",
                 "idr",
                 "simulateData",
                 "spearman")
sapply(helperFuns, function(x) source(paste0("./code/",x,".R")))

## simulate data
params <- expand.grid(p1=0.05,
                      p2=0.05,
                      p3=c(0,0.5,1),
                      p4=c(0,0.5,1))
params <- params[params$p3+params$p4<=1,]
params_names <- apply(params, 1, function(x) paste("p1",x[1],"p2",x[2],"p3",x[3], "p4", x[4], sep="_"))
set.seed(1)
dat <- lapply(1:nrow(params),
              function(i) {
                print(paste0("simulating ", nSim, " simulations of ", params_names[i]))
                paramSims <- lapply(1:nSim, function(x) {
                  simulateData(p1=params[i,1], p2=params[i,2], p3=params[i,3], p4=params[i,4])
                })
                names(paramSims) <- paste(params_names[i], "sim", 1:nSim, sep="_")
                return(paramSims)
              })
names(dat) <- params_names

## assess significance
simFuns <- c("gsea", "IDR.func", "spearman", "trunSpearman")
simFunsDir <- c(F, T, F, F)
## stupidity check
print(paste("utilizing", nCores, "cores"))
registerDoParallel(nCores)
simTime <- system.time({tmp <- assessSig(dat[[1]][[1]], simFuns=simFuns, simFunsDir=simFunsDir, nPerm=5)})
timePerSim <- round(simTime[3]*nPerm/5/60,2)
print(paste("time for 5 permutations:", round(simTime[3],2), "sec"))
print(paste("expected time per simulation:", timePerSim, "min"))
print(paste("expected time per parameter set:", round(timePerSim*nSim/60,2), "hr"))
print(paste("expected total time for", length(params_names), "parameter sets:",
            round(timePerSim*nSim/60*length(params_names),2), "hr"))
# begin permutations
for(i in 1:length(params_names)) {
  print(paste("assessing for", params_names[i]))
  tmp <- as.list(rep(NA, length(dat[[i]])))
  for(j in 1:nSim) {
    print(paste("permuting for simulation", j, "of", nSim))
    tmp[[j]] <- assessSig(dat[[i]][[j]], simFuns=simFuns, simFunsDir=simFunsDir, nPerm=nPerm)
  }
  assign(params_names[i], value=tmp)
  save(list=c(params_names[i]), file=paste0("./data/", params_names[i], ".Rda"))
}
stopImplicitCluster()
