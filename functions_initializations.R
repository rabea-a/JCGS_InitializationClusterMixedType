


### help functions

lambda_kproto <- function(x, lambda = NULL, verbose = FALSE){
  
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  
  # determination of lambda
  if(anynum & anyfact){
    vnum <- mean(sapply(x[,numvars, drop = FALSE], var, na.rm = TRUE))
    vcat <- mean(sapply(x[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
    if (vnum == 0){
      if(verbose) warning("All numerical variables have zero variance.")
      anynum <- FALSE
    } 
    if (vcat == 0){
      if(verbose) warning("All categorical variables have zero variance.")
      anyfact <- FALSE
    } 
    if(anynum & anyfact){
      lambda <- vnum/vcat
      if(verbose) cat("Estimated lambda:", lambda, "\n\n")
    }else{
      lambda <- 1
    }
  }
  
  return(lambda)
}



dists_kproto <- function(x, y = NULL, lambda = NULL, verbose = FALSE){
  
  if(is.null(y)){
    dat_part1 <- x[rep(c(1:nrow(x)), each = nrow(x)),]
    dat_part2 <- x[rep(c(1:nrow(x)), times = nrow(x)),]
  }else{
    
    if(nrow(x) == 1 & nrow(y) == 1){
      return(cbind(x, y, dist = dist_kproto(x,y,lambda = lambda)))
    }else{
      dat_part1 <- x[rep(c(1:nrow(x)), each = nrow(y)),]
      dat_part2 <- y[rep(c(1:nrow(y)), times = nrow(x)),]
    }
    
  }
  
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  
  # determination of lambda
  if(length(lambda) > 1) {if(length(lambda) != sum(c(numvars,catvars))) stop("If lambda is a vector, its length should be the sum of numeric and factor variables in the data frame!")}
  if(is.null(lambda)){
    if(anynum & anyfact){
      vnum <- mean(sapply(x[,numvars, drop = FALSE], var, na.rm = TRUE))
      vcat <- mean(sapply(x[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
      if (vnum == 0){
        if(verbose) warning("All numerical variables have zero variance.")
        anynum <- FALSE
      } 
      if (vcat == 0){
        if(verbose) warning("All categorical variables have zero variance.")
        anyfact <- FALSE
      } 
      if(anynum & anyfact){
        lambda <- vnum/vcat
        if(verbose) cat("Estimated lambda:", lambda, "\n\n")
      }else{
        lambda <- 1
      }
    }
  }
  
  # compute distances 
  nrows <- nrow(x)
  d1 <- (dat_part1[,numvars, drop = FALSE] - dat_part2[,numvars, drop = FALSE])^2
  if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
  if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
  d2 <- sapply(which(catvars), function(j) return(dat_part1[,j] != dat_part2[,j]))
  d2[is.na(d2)] <- FALSE
  if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
  if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
  
  return(cbind(dat_part1, dat_part2, dist = as.vector(d1 + d2)))
}



dist_kproto <- function(x, y, lambda, verbose = FALSE){
  
  x <- rbind(x,y)
  
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  
  # compute distances 
  nrows <- nrow(x)
  d1 <- (x[1, numvars, drop = FALSE] - x[2, numvars, drop = FALSE])^2
  if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
  if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
  d2 <- sapply(which(catvars), function(j) return(x[1,j] != x[2,j]))
  d2[is.na(d2)] <- FALSE
  if(length(lambda) == 1) d2 <- lambda * rowSums(matrix(d2, nrow = 1))
  if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]

  return(as.numeric(d1 + d2))
}





# determine number of randomly choosen sets, 
#   before first one is successful with a desirable set:
nstart.m <- function(p_start, N, k_true){
  # N = number of objects
  # k_true = number of true clusters/cluster to be determined
  # p_start = probability, that minimum one desirable set is randomly choosen
  
  prob <- 1
  for(i in 0:(k_true-1)){
    prob <- prob * (N - i * N/k_true)/(N - i)
  }
  
  # AnmR: pgeom describes x number of failures and x+1 number of failures until 1 success..
  #       => therefore: qgeom + 1
  return(qgeom(p = p_start, prob = prob)+1)
  
}





# based on Kaufman, Rousseeuw (1990)
inits_sel.dist <- function(dat, k){
  
  # determine the first prototype, which is the most centrally located object
  # therefor determine the sum of distances to all other objects 
  # select the object with the smallest sum
  lambda <- lambda_kproto(x = dat)
  all_dists <- dists_kproto(x = dat, lambda = lambda)
  dists_matrix <- matrix(all_dists[,"dist"], ncol = nrow(dat), nrow = nrow(dat)) - diag(NA, nrow = nrow(dat))
  protos_initial <- dat[which.min(colSums(dists_matrix, na.rm = TRUE)),]
  
  while(nrow(protos_initial) < k){
    # determine distances to nearest prototype
    dists_protos <- dists_kproto(x = dat, y = protos_initial, lambda = lambda)
    dists_protos <- apply(X = matrix(dists_protos[,"dist"], ncol = nrow(dat), nrow = nrow(protos_initial)),
                          MARGIN = 2, FUN = min)
    
    # sum Cij in notation of Kaufman Rousseeuw
    sum_changes <- apply(X = cbind(dists_matrix, dists_protos), MARGIN = 1, 
                         FUN = function(x) sum((x[nrow(dat)+1] - x[-(nrow(dat)+1)])[x[-(nrow(dat)+1)] < x[nrow(dat)+1]], na.rm = TRUE))
    # determine the next prototype
    protos_initial <- rbind(protos_initial, dat[which.max(sum_changes),])
  }
  
  return(protos_initial)
  
}






# based on Katsavounidis et al. (1994) 
inits_sel.cen <- function(dat, k){
  
  # determine object which has the highest 
  # sum of distances to all other objects as first prototypes
  lambda <- lambda_kproto(x = dat)
  all_dists <- dists_kproto(x = dat, lambda = lambda)
  dists_matrix <- matrix(all_dists[,"dist"], ncol = nrow(dat), nrow = nrow(dat)) - diag(NA, nrow = nrow(dat))
  protos_initial <- dat[which.max(colSums(dists_matrix, na.rm = TRUE)),]
  
  # determine minimal distance to initial prototypes
  while(nrow(protos_initial) < k){
    Dmin <- numeric()
    for(j in 1:nrow(dat)){
      Dmin[j] <- min(dists_kproto(x = dat[j,, drop = FALSE], y = protos_initial, lambda = lambda)$dist)
    }
    
    # choose next initial prototype by highest distance to already choosen prototypes
    protos_initial <- rbind(protos_initial, dat[which.max(Dmin),])
    
  }
  
  return(protos_initial)
  
}





# based on Arthur, Vassilvitski (2007)
inits_prob.cen <- function(dat, k){
  
  # determine object to save initial prototypes
  protos_initial <- dat[sample(1:nrow(dat),1),]
  
  # determine lambda for distance calculation
  lambda <- lambda_kproto(x = dat)
  
  # determine minimal distance to initial prototypes
  while(nrow(protos_initial) < k){
    Dmin <- numeric()
    for(j in 1:nrow(dat)){
      Dmin[j] <- min(dists_kproto(x = dat[j,, drop = FALSE], y = protos_initial, lambda = lambda)$dist)
    }
    
    # choose next initial prototype by probability (Dmin^2/sum(Dmin^2))
    protos_initial <- rbind(protos_initial, 
                            dat[sample(x = 1:nrow(dat), size = 1, prob = Dmin^2/sum(Dmin^2)),])
    
  }
  
  return(protos_initial)
  
}





# based Xu et al. (2009)
inits_rnn <- function(dat, k){
  
  # determine distances between all pairs of objects:
  lambda <- lambda_kproto(x = dat)
  all_dists <- dists_kproto(x = dat, lambda = lambda)
  dists_matrix <- matrix(all_dists[,"dist"], ncol = nrow(dat), nrow = nrow(dat)) - diag(NA, nrow = nrow(dat))
  
  # selected prototype candidates
  protos_cs <- dat
  
  while(nrow(protos_cs) > (5*k)){
    
    # determine object to save initial prototypes
    protos_initial <- dat[0,]
    
    # determine RNN
    RNN <- vector("list", nrow(protos_cs))
    names(RNN) <- row.names(protos_cs)
    for(j in row.names(protos_cs)){
      # which are the objects of which the object is the nearest neighbor
      RNN[[j]] <- as.numeric(row.names(protos_cs))[which(apply(X = dists_matrix[as.numeric(row.names(protos_cs)), as.numeric(row.names(protos_cs))], 
                                                               MARGIN = 1, FUN = function(x) x[row.names(protos_cs) %in% j] == min(x, na.rm = TRUE)))]
    }
    
    RNN_del <- RNN
    RNN_del <- RNN_del[(names(RNN_del) %in% names(RNN_del)[lapply(RNN, length) > 0])]  
    repeat{
      obj_index <- names(RNN_del)[sort.int(unlist(lapply(RNN_del, length)), decreasing = TRUE, index.return = TRUE)$ix[1]]
      
      protos_initial <- rbind(protos_initial, protos_cs[obj_index,])
      protos_cs <- protos_cs[-which(row.names(protos_cs) %in% c(obj_index, RNN_del[[obj_index]])),]
      
      RNN_del <- RNN_del[!(names(RNN_del) %in% c(obj_index, RNN_del[[obj_index]]))]
      
      if(length(RNN_del) == 0){break}
    }
    
    # next iteration with only the protos candidates
    protos_cs <- protos_initial
  }
  
  if(nrow(protos_cs) == k){
    protos_selected <- protos_cs
  }else{
    # determination of initial prototypes by application of RFN
    protos_selected <- dat[0,]
    for(j in 1:k){
      
      # determine RNN
      RNN <- vector("list", nrow(protos_initial))
      names(RNN) <- row.names(protos_initial)
      for(l in row.names(protos_initial)){
        # which are the objects of which the object is the nearest neighbor
        RNN[[l]] <- as.numeric(row.names(protos_initial))[which(apply(X = dists_matrix[as.numeric(row.names(protos_initial)), as.numeric(row.names(protos_initial))],
                                                                      MARGIN = 1, FUN = function(x) x[row.names(protos_initial) %in% l] == min(x, na.rm = TRUE)))]
      }
      
      # determine RFN
      RFN <- vector("list", nrow(protos_initial))
      names(RFN) <- row.names(protos_initial)
      for(h in row.names(protos_initial)){
        # which are the objects of which the object is the nearest neighbor
        RFN[[h]] <- as.numeric(row.names(protos_initial))[which(apply(X = dists_matrix[as.numeric(row.names(protos_initial)), as.numeric(row.names(protos_initial))], 
                                                                      MARGIN = 1, FUN = function(x) x[row.names(protos_initial) %in% h] == max(x, na.rm = TRUE)))]
      }
      
      # choose object with highest RNF
      obj_index <- names(RFN)[sort.int(unlist(lapply(RFN, length)), decreasing = TRUE, index.return = TRUE)$ix[1]]
      protos_selected <- rbind(protos_selected, protos_initial[obj_index,])
      
      # delete the choosen object and the respective RNNs
      protos_initial <- protos_initial[-which(row.names(protos_initial) %in% c(obj_index, RNN[[obj_index]])),]
      
      # break, if number of possible objects is equal to k
      if((nrow(protos_initial) + nrow(protos_selected)) == k){
        protos_selected <- rbind(protos_selected, protos_initial)
        break
      }
      
      # if not enough protos_intial are avaiable to choose
      if(nrow(protos_initial) == 0){
        
        if((k - nrow(protos_selected)) != 0){
          k_needed <- k - nrow(protos_selected)
          protos_selected <- rbind(protos_selected, protos_cs[row.names(protos_cs)[-which(row.names(protos_cs) %in% row.names(protos_selected))][1:k_needed],])
        }
        
        break
      }
    }
  }
  
  return(protos_selected)
  
}





# based on Ji et al. (2015)
inits_nbh.cen <- function(dat, k){
  
  # determine distances between all pairs of objects:
  lambda <- lambda_kproto(x = dat)
  all_dists <- dists_kproto(x = dat, lambda = lambda)
  
  # determine sigma (as the average distance between all pairs of objects)
  sigma <- (sum(all_dists$dist)/(nrow(all_dists)-nrow(dat)))
  
  # determine NborS, which is the number of neighbors in the neighborhood
  NborS <- numeric()
  for(i in 1:nrow(dat)){
    # noteR: -1 because the distance between X_i and X_i is included and always < nbh_radius...
    NborS[i] <- sum(all_dists[which(rowSums(mapply("==", dat[i,], all_dists[, 1:ncol(dat)])) == ncol(dat)), ]$dist <= sigma) - 1
  }
  
  # save centrality (Cen) for every object
  Cen <- NborS/max(NborS)
  
  # determine object to save initial prototypes
  protos_initial <- dat[sort.int(Cen, index.return = TRUE, decreasing = TRUE)$ix[1],]
  
  while(nrow(protos_initial) < k){
    
    # determine the highest probability of being the next initial prototype
    # (which is minimum dist to the nearest initial prototype multiplied by object's centrality)
    Pro <- numeric()
    for(j in 1:nrow(dat)){
      Pro[j] <- min(dists_kproto(x = dat[j,, drop = FALSE], y = protos_initial, lambda = lambda)$dist * Cen[j])
    }
    
    # determine the highest probability of being the next initial prototype
    protos_initial <- rbind(protos_initial, dat[sort.int(Pro, index.return = TRUE, decreasing = TRUE)$ix[1],])
    
  }
  
  return(protos_initial[1:k,])
  
}





# based on Guo et al. (2018)
intis_nbh.dens <- function(dat, k){
  
  # determine object to save initial prototypes
  protos_initial <- dat[0,]
  
  # determine distances between all pairs of objects:
  lambda <- lambda_kproto(x = dat)
  all_dists <- dists_kproto(x = dat, lambda = lambda)
  
  # determine \bar{d}/2 (where \bar{d} is the average distance between all pairs of objects)
  nbh_radius <- (sum(all_dists$dist)/(nrow(all_dists)-nrow(dat)))/2
  
  # determine N_e, which is the sum of h() for every object 
  # Ne_nbh says the number of objects with distance smaller than nbh_radius for every object
  Ne_nbh <- numeric()
  for(i in 1:nrow(dat)){
    # noteR: -1 because the distance between X_i and X_i is included and always < nbh_radius...
    Ne_nbh[i] <- sum(all_dists[which(rowSums(mapply("==", dat[i,], all_dists[, 1:ncol(dat)])) == ncol(dat)), ]$dist <= nbh_radius) - 1
  }
  
  # determine first initial prototype
  arrange_index <- order(Ne_nbh, decreasing = TRUE)
  protos_initial[1,] <- dat[arrange_index[1],]
  
  # determination of initial prototypes
  # lowering L until we have enough initial objectss found
  L <- nbh_radius*4*2
  while(nrow(protos_initial) < k){
    # determine the other initial prototypes by comparison of objects to yet known prototypes

    L <- L/2
    
    for(i in arrange_index[-1]){
      dists <- dists_kproto(x = dat[i,, drop = FALSE], y = protos_initial, lambda = lambda)
      # noteR: paper: !"...determine whether any M_i \in M has dist(Y_i,M_i) > L holds.
      #        If satisfied, then Y_i is put into the set M as next prototype..
      #        NOT MEANINGFUL, THE DISTANCE TO ALL PROTOS MUST BE > L...
      if(all(dists$dist > L)){
        protos_initial <- rbind(protos_initial, dat[i,])
      }
    }
  }
  
  return(protos_initial[1:k,])
  
}













# # Reddy, Jana (kMeans): voronoi circles
# inits_vcircles_reddy <- function(dat, k){
#   
# }


### Mail on 20221107: explanation of determination of L (formular 17)
# ### Jia and Song: 
# intis_nbh_jiasong <- function(dat, k){
#   
#   # determine object to save initial prototypes
#   protos_initial <- dat[0,]
#   
#   # determine distances between all pairs of objects:
#   lambda <- lambda_kproto(x = dat)
#   all_dists <- dists_kproto(x = dat, lambda = lambda)
#   
#   # # average distance between all pairs of objects
#   bar_d <- sum(all_dists$dist)/(nrow(all_dists)-nrow(dat))
#   
#   # determine local neighborhood density
#   # rho_i for every object, which is the sum of dists to other objects which are <= d_c
#   # "d_c is a critical value that limits the search scope"
#   d_c <- bar_d #AnmR.: willkuerlich! was ist sinnvoll???
#   rho <- numeric()
#   for(i in 1:nrow(dat)){
#     # noteR: -1 because the distance between X_i and X_i is included and always < d_c...
#     rho[i] <- sum(all_dists[which(rowSums(mapply("==", dat[i,], all_dists[, 1:ncol(dat)])) == ncol(dat)), ]$dist <= d_c) - 1
#   }
#   dat_plus <- cbind(dat, rho)
#   dat_plus <- dat_plus[sort.int(rho, index.return = TRUE, decreasing = TRUE)$ix,]
#   
#   
#   protos_initial <- dat_plus[1,]
#   for(j in 2:nrow(dat)){
#     
#     # determine distance threshold
#     for(q in 1:nrow(protos_initial)){
#       L <- ifelse(###noteR: which condition? iterates x_j and x_i is fixed??,
#                   #min dist to any other object of dat
#                   min(dists_kproto(x = dat_plus[j,1:ncol(dat)], y = dat_plus[-j,1:ncol(dat)], lambda = lambda)[,"dist"]),
#                   #max dist to any other object of dat
#                   max(dists_kproto(x = dat_plus[j,1:ncol(dat)], y = dat_plus[-j,1:ncol(dat)], lambda = lambda)[,"dist"])
#                   )
#       
#       if(all(dists_kproto(x = dat_plus[j,1:ncol(dat)], y = protos_initial[,1:ncol(dat)], lambda = lambda)[,"dist"] < L)){
#         protos_initial <- rbind(protos_initial,dat_plus[j,])
#       }
#       
#     }
#     
#   }
#   
#   
# }





