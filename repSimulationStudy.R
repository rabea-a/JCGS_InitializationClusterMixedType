

library(clustMixType)
source("functions_initializations.R")

#########################################################################################
### SimStudy ############################################################################
#########################################################################################

simstudy_dat <- readRDS(url("https://staffcloud.hochschule-stralsund.de/index.php/s/Y8fzJfw3fSyc9ds/download/repStudy_data.rds"))
simstudy_trial_design <- readRDS("repStudy_trial_design.rds")



simstudy_results <- vector("list", nrow(simstudy_trial_design))
for(i in 1:nrow(simstudy_trial_design)){
  k <- simstudy_trial_design[i, "nC"]   # how many protos should be determined?
  simstudy_results[[i]] <- vector("list", simstudy_trial_design[i,"N"])
  
  for(m in 1:simstudy_trial_design[i,"N"]){
    
    ### object to save all results
    simstudy_results[[i]][[m]] <- vector("list", 10)
    #  1 = "nstart = 1" 
    #  2 = "nstart = 3" 
    #  3 = "nstart = 10"
    #  4 = "nstart based on expected cluster"
    #  5 = "selection by max. distance" (Kaufman, Rousseeuw; 1990)
    #  6 = "selection by centrality" (Katsavounidis et al.; 1994)
    #  7 = "probability of centrality" (Arthur, Vassilvitskii; 2007)
    #  8 = "reverse nearest neighbors" (Xu et al.; 2009)
    #  9 = "neighborhood centrality" (Ji et al.; 2015)
    # 10 = "neighborhood density" (Guo et al.; 2018)
    
    
    
    ### data to be clustered:
    dat <- simstudy_dat[[i]][[m]][,-1]
    
    # saving for each approach the following:
    # 1 = initial prototypes
    # 2 = resulting kproto-object
    # 3 = runtime
    # 4 = tot.withinss: sum of all observations' distances to their corresponding cluster prototype
    # 5 = adjusted rand index
    
    
    
    ### 1 - nstart.1, 2 - nstart.3, 3 - nstart.10:
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[1]][[2]] <- kproto(dat, k = k, nstart = 1, verbose = FALSE)
    simstudy_results[[i]][[m]][[1]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[1]][[4]] <- simstudy_results[[i]][[m]][[1]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[1]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[1]][[2]]$cluster)
    
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[2]][[2]] <- kproto(dat, k = k, nstart = 3, verbose = FALSE)
    simstudy_results[[i]][[m]][[2]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[2]][[4]] <- simstudy_results[[i]][[m]][[2]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[2]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[2]][[2]]$cluster)
    
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[3]][[2]] <- kproto(dat, k = k, nstart = 10, verbose = FALSE)
    simstudy_results[[i]][[m]][[3]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[3]][[4]] <- simstudy_results[[i]][[m]][[3]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[3]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[3]][[2]]$cluster)
    
    
    
    ### 4 - nstart.m (nstart determined by probability to include a desirable set):
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[4]][[1]] <- nstart.m(p_start = 0.9, N = simstudy_trial_design[i, "nO"],
                                                    k_true = k, n_protos = k)
    simstudy_results[[i]][[m]][[4]][[2]] <- kproto(dat, k = k, verbose = FALSE, 
                                                   nstart = simstudy_results[[i]][[m]][[4]][[1]])
    simstudy_results[[i]][[m]][[4]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[4]][[4]] <- simstudy_results[[i]][[m]][[4]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[4]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[4]][[2]]$cluster)
    
    
    
    
    
    ### 5 - sel.dist (based on Kaufman, Rousseeuw, 1990):
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[5]][[1]] <- inits_sel.dist(dat, k)
    simstudy_results[[i]][[m]][[5]][[2]] <- kproto(dat, verbose = FALSE, nstart = 1,
                                                   k = simstudy_results[[i]][[m]][[5]][[1]])
    simstudy_results[[i]][[m]][[5]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[5]][[4]] <- simstudy_results[[i]][[m]][[5]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[5]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[5]][[2]]$cluster)
    
    
    
    
    
    ### 6 - (based on Katsavounidis et al., 1994):
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[6]][[1]] <- inits_sel.cen(dat, k)
    simstudy_results[[i]][[m]][[6]][[2]] <- kproto(dat, verbose = FALSE, nstart = 1,
                                                   k = simstudy_results[[i]][[m]][[6]][[1]])
    simstudy_results[[i]][[m]][[6]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[6]][[4]] <- simstudy_results[[i]][[m]][[6]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[6]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[6]][[2]]$cluster)
    
    
    
    
    
    ### 7 - prob.cen (based on Arthur, Vassilvitskii, 2007):
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[7]][[1]] <- inits_prob.cen(dat, k)
    simstudy_results[[i]][[m]][[7]][[2]] <- kproto(dat, verbose = FALSE, nstart = 1,
                                                   k = simstudy_results[[i]][[m]][[7]][[1]])
    simstudy_results[[i]][[m]][[7]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[7]][[4]] <- simstudy_results[[i]][[m]][[7]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[7]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[7]][[2]]$cluster)
    
    
    
    
    
    ### 8 - rnn (based on Xu et al., 2009):
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[8]][[1]] <- inits_rnn(dat, k)
    simstudy_results[[i]][[m]][[8]][[2]] <- kproto(dat, verbose = FALSE, nstart = 1,
                                                   k = simstudy_results[[i]][[m]][[8]][[1]])
    simstudy_results[[i]][[m]][[8]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[8]][[4]] <- simstudy_results[[i]][[m]][[8]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[8]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[8]][[2]]$cluster)
    
    
    
    
    
    ### 9 - nbh.cen (based on Ji, Pang et al., 2015):
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[9]][[1]] <- inits_nbh.cen(dat, k)
    simstudy_results[[i]][[m]][[9]][[2]] <- kproto(dat, verbose = FALSE, nstart = 1,
                                                   k = simstudy_results[[i]][[m]][[9]][[1]])
    simstudy_results[[i]][[m]][[9]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[9]][[4]] <- simstudy_results[[i]][[m]][[9]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[9]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[9]][[2]]$cluster)
    
    
    
    
    
    ### 10 - nbh.dens (based on Guo et al., 2018): 
    time <- Sys.time()
    set.seed(221122)
    simstudy_results[[i]][[m]][[10]][[1]] <- intis_nbh_density(dat, k)
    simstudy_results[[i]][[m]][[10]][[2]] <- kproto(dat, verbose = FALSE, nstart = 1,
                                                   k = simstudy_results[[i]][[m]][[10]][[1]])
    simstudy_results[[i]][[m]][[10]][[3]] <- as.numeric(Sys.time() - time, units="secs")
    simstudy_results[[i]][[m]][[10]][[4]] <- simstudy_results[[i]][[m]][[10]][[2]]$tot.withinss
    simstudy_results[[i]][[m]][[10]][[5]] <- fossil::adj.rand.index(simstudy_dat[[i]][[m]][,1], 
                                                                   simstudy_results[[i]][[m]][[10]][[2]]$cluster)

    
    
    
    
    
    #cat(paste0("Done: ", i, " trial design out of ", nrow(simstudy_trial_design)," and N=", m, "\n"))
  }
}








