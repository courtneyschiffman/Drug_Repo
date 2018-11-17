nSim <- 100 # number of simulations to perform
nPerm <- 1000 # number of permutations to perform

## source script files
scriptFiles <- grep(".R", list.files(path="./code"), value=T)
scriptFiles <- grep("main", scriptFiles, invert=T, value=T)
sapply(scriptFiles, function(x) source(paste0("./code/",x)))

## generate simulated data
params <- expand.grid(p1=c(0.05,0.20),
                      p2=0.05,
                      p3=c(0,0.25,0.5,0.75,1))
params_names <- apply(params, 1, function(x) paste("p1",x[1],"p2",x[2],"p3",x[3], sep="_"))
set.seed(1)
dat <- lapply(1:nrow(params),
              function(i) {
                print(paste0("simulating ",params_names[i]))
                lapply(1:nSim, function(x) {
                  simulateData(p1=params[i,1], p2=params[i,2], p3=params[i,3])
                })
              })
names(dat) <- params_names
# save(dat, file="./data/simulatedData.Rda")

## perform permutations
simFuns <- c("gsea", "spearman", "trunSpearman")
# simFuns <- c("gsea", "IDR.func", "spearman", "trunSpearman")
dat_permutations <- lapply(dat,
                           function(paramSet) {
                             lapply(paramSet,
                                    function(perm) {
                                      assessSig(perm, simFuns=simFuns, nPerm=2, deg=49)
                                    })
                           })

tmp <- lapply(dat[[1]],
              function(x) {
                assessSig(x, simFuns=simFuns, nPerm=2, deg=49)
              })
tmp <- lapply(1:length(dat[[1]]),
              function(i) {
                print(i)
                do.call("gsea", args=list(x=dat[[1]][[i]]$x,
                                          y=dat[[1]][[i]]$y,
                                          deg=49))
              })
gsea(dat[[1]][[10]]$x, dat[[1]][[10]]$y, deg=49)
tmp <- dat[[1]][[10]]
