#### Simulation Functions ####


# Helpful functions -----------------------------------------------------------------
`%nin%` <- Negate(`%in%`)

pkgs <- c("tidyverse", "magrittr", "igraph", "matrixcalc",
          "MASS", "diffusr", "Matrix", "cluster", 'mclust',
          'Spectrum', 'NbClust')

suppressMessages(lapply(pkgs, library, character.only = T))

# Kernels --------------------------------------------------------

Gaussian_kernel <- function(Z, rho){
  exp(-(1/rho)*as.matrix(dist(Z, method = "euclidean", upper = T)^2))
}

Zhang_kernel <- function(Z, p){
  d <- as.matrix( dist(Z) )
  
  si <- apply(d, 2, function(s) sort(s)[p+1] )
  db <- si %o% si
  dt <- -d^2/db
  
  exp(dt)
}

knnGrf <- function(A, k, mutual=FALSE){
  knn <- function(aa, k){
    aa[order(aa)[(k+2):length(aa)]] <- 0
    aa[aa > 0] <- 1
    aa
  }
  
  akn <- apply(A,2,knn,k=k)
  
  if(mutual){
    for(i in 1:nrow(akn)){
      for(j in 1:ncol(akn)){
        akn[i,j] <- ifelse(akn[i,j]==akn[j,i], 
                           akn[i,j],
                           0)
      }
    }
  } else{
    for(i in 1:nrow(akn)){
      for(j in 1:ncol(akn)){
        akn[i,j] <- ifelse(akn[i,j]!=akn[j,i], 
                           max(akn[i,j],akn[j,i]),
                           akn[i,j])
      }
    }
  }
  akn
}

## Laplacian ## USED FOR HEATMAPS
gaussLaplac <- function(dat, rho=NULL, kernel='Gaussian',
                        lap.type = 'sym',
                        grf.type = 'full', k=5, p=5,
                        epsilon=0, mutual=FALSE, 
                        binary.grf=FALSE,
                        plots=TRUE, verbose=TRUE){
  
  if(is.null(rho)) rho <- median(dist(dat))
  if(kernel == 'Gaussian') A <- Gaussian_kernel(dat, rho = rho)
  
  if(kernel == 'Zhang') A <- Zhang_kernel(dat, p)
  if(kernel == 'Spectrum'){
    dtd <- as.data.frame( t(dat) )
    A <- Spectrum::CNN_kernel(dtd)
  }
  
  if(kernel == 'Linear') A <- dat%*%t(dat)
  if(kernel == 'Cor') A <- cor(t(dat))
  
  diag(A) <- 0 ## Adjacency matrices need diag of 0
  
  if(verbose){
    message(paste('Distances calculated with', kernel, 'kernel.'))
    message("Summary of parwise distances:")
    print(summary(A[lower.tri(A)]))
  }
  
  if(grf.type == 'e-graph'){
    A[A < epsilon] <- 0
  }
  
  if(grf.type == 'knn'){
    A <- knnGrf(A, k, mutual=mutual)
  }
  
  if(binary.grf) A[A >0] <- 1
  
  deg <- rowSums(A)
  
  if(lap.type == 'sym'){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = I - D^(-1/2) A D^(-1/2)
    L <- diag(nrow = nrow(A)) - diag(ds) %*% A %*% diag(ds)
  }
  
  ## 'Shifted' Laplacian
  if(lap.type == "shift"){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = I + D^(-1/2) A D^(-1/2)
    L <- diag(nrow = nrow(A)) + diag(ds) %*% A %*% diag(ds)
  }
  
  ## Laplacian used by Ng (2002) 'On Spectral Clustering'
  if(lap.type == "Ng"){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = D^(-1/2) A D^(-1/2)
    L <- diag(ds) %*% A %*% diag(ds)
  }
  
  ## Random walk Laplacian (best one)
  if(lap.type == 'rw'){
    ds <- ifelse(deg>0, 1/deg, 0)
    # L = I - D^-1 A
    L <- diag(nrow = nrow(A)) - diag(ds) %*% A
  } 
  
  if(plots){
    gg <- igraph::graph_from_adjacency_matrix(A,
                                              mode='undirected',
                                              weighted = TRUE)
    # plot(gg, layout=igraph::layout_with_kk,
    #      vertex.color='darkorchid1',
    #      vertex.label=NA)
    
    stats::heatmap(A, symm=TRUE, 
                   main = "Heatmap of Adjacency Matrix")
    
    plot(eigen(L)$values, col='dodgerblue3', pch=20, type="b",
         ylab="Eigen values", main = "Eigen values of the Laplacian")
  }
  
  L
}

# Matrix Functions --------------------------------------------------------

## Trace of a matrix
tr <- function(x){
  if(!is.matrix(x)) stop("x must be a matrix")
  sum(diag(x))
}

## Check if matrix is orthogonal
orthoCheck <- function(x){
  xtx <- crossprod(x)
  I <- diag(nrow = nrow(xtx))
  
  # Same tolerance as 'isSymmetric'
  sum((I-xtx)^2) < 100*.Machine$double.eps
}

## Row normalize a matrix (rows have length 1)
rowNrm <- function(x){
  nrm <- apply(x, 1, function(xi) sqrt(sum(xi^2)))
  sweep(x, 1, nrm, "/")
}

# Flag Mean ---------------------------------------------------------------

## Calculate flag mean subspace
## Marrinan, et al. "Finding the Subspace Mean or Median to Fit Your Need"
flagMean <- function(Ul, plots=TRUE){
  
  ## Check orthogonal
  # if(!all(sapply(Ul, orthoCheck))) stop("All matrices must have orthogonal columns")
  
  r <- max(vapply(Ul, ncol, 0))
  
  ## Concatenating matrices
  X <- do.call(cbind, Ul)
  svdX <- svd(X)
  
  if(plots) plot(svdX$d, col='dodgerblue3', pch=20, type="b",
                 ylab="Singular Values", main="Singular values of the flag mean")
  
  svdX
}

flagMeanClust <- function(Llist, kList, ngrp, lap.type='rw'){
  subSpaceL <- mapply(subSpace, Llist, kList,
                      lap.type=lap.type, SIMPLIFY=FALSE)
  fm <- flagMean(subSpaceL, plots = FALSE)$u
  
  km <- kmeans(fm[,1:ngrp], centers=ngrp, nstart=10)
  list(eLR=fm[,1:ngrp], km=km)
}

subSpace <- function(L,k, lap.type='rw'){
  eig <- eigen(L)
  
  vec <- eig$vector 
  if(lap.type %in% c('rw','sym')) ind <- (ncol(vec)-k+1):ncol(vec)
  else ind <- 1:k
  vec[,ind,drop=F]
}

## Flag mean ignoring constant vector. Changes the 'ind'
flagMeanNoConstant <- function(Llist, kList, ngrp){
  subSpace <- function(L,k){
    eig <- eigen(L)
    
    vec <- eig$vector 
    ind <- (ncol(vec)-k+1):(ncol(vec)-1)
    vec[,ind, drop=FALSE]
  }
  
  subSpaceL <- mapply(subSpace, Llist, kList)
  fm <- flagMean(subSpaceL, plots = FALSE)$u
  
  km <- kmeans(fm[,1:(ngrp-1)], centers=ngrp, nstart=10)
  list(eLR=fm[,1:(ngrp-1)], km=km)
}

# L_approx ----------------------------------------------------------------

## Function to extract k primary eigen vectors and values
## and reconstruct Laplacians from that
LapSig <- function(Ls, ks, lap.type='rw', normalize=FALSE){
  
  eig <- eigen(Ls)
  
  if(lap.type %in% c('rw', 'sym')){
    ind <- (ncol(Ls)-ks+1):ncol(Ls)
  } else ind <- seq(1,ks)
  
  vecs <- eig$vectors[,ind]
  vals <- eig$values[ind]
  
  if(normalize) vecs <- rowNrm(vecs)
  
  vecs %*% diag(vals) %*% t(vecs)
}

## Better way to do that but keeping the above one for old code
Lapprox <- function(L, k, lap.type='rw'){
  eig <- eigen(L)
  
  vec <- eig$vector
  if(lap.type %in% c('rw','sym')) ind <- (ncol(vec)-k+1):ncol(vec)
  else ind <- 1:k
  
  ev <- vec[,ind,drop=F]; val <- eig$values[ind]
  
  sweep(ev, 2, val, "*")%*%t(ev)
}

LapproxClust <- function(Llist, Klist, ngrp, lap.type='rw'){
  # ngrp <- max(Klist)
  approxList <- mapply(Lapprox, Llist, Klist, 
                       lap.type=lap.type, SIMPLIFY=FALSE)
  approxList <- lapply(approxList, function(aL) aL/length(approxList))
  
  LR <- Reduce("+", approxList)
  eLR <- eigen(LR)$vectors
  km <- kmeans(eLR[,1:ngrp], centers=ngrp, nstart=10)
  
  list(eLR=eLR[,1:ngrp], km=km)
}


# ClusterEval -------------------------------------------------------------

silhouette_score <- function(vecs, k, nstart=10){
  km <- kmeans(vecs, centers = k, nstart=nstart)
  
  ss <- silhouette(km$cluster, dist(vecs))
  mean(ss[, 3])
}

## F-score / F-measure / F1
f_measure <- function(ll, grps, k, type = 'simple'){
  
  if(type == 'simple') ev <- eigen(ll)$vectors[,1:k]
  else if(type == 'flag') ev <- ll$u[,1:k]
  else stop("'type' should be 'simple' or 'flag' ")
  
  km <- kmeans(ev, centers = k, nstart = k)
  
  cb <- cbind(km$cluster, grps)
  tt <- table(cb[,1], cb[,2])
  
  ## Randomly distributing wrong groupings
  ## to other incorrect places for square matrix
  while(nrow(tt) != ncol(tt)){
    tm <- which.min(colSums(tt))
    miss <- colSums(tt)[tm]
    
    tt <- tt[, -tm]
    zs <- which(tt==0)
    
    for(mm in 1:miss){ 
      ss <- sample(zs, 1)
      tt[ss] <- tt[ss] + 1
    }
  }
  
  ## Reordering so highest agreement is on diagonal
  ord <- order(apply(tt, 2, which.max))
  tt <- tt[, ord]
  colnames(tt) <- 1:ncol(tt)
  rownames(tt) <- 1:nrow(tt)
  
  ## In case of imperfect groupings
  # F1 score / F-measure
  f1 <- confusionMatrix(tt, mode = 'prec_recall')$byClass[,"F1"]
  f1[is.nan(f1)] <- 0
  mean(f1)
}

fMeas <- function(ll, grps, k, type='simple'){
  tt <- table(ll$km$cluster, grps)
  
  ## Making confusion matrix square
  while(nrow(tt) != ncol(tt)){
    if(nrow(tt) > ncol(tt)){
      tt <- cbind(tt, rep(0, nrow(tt)))
    }
    if(ncol(tt) > nrow(tt)){
      tt <- rbind(tt, rep(0, ncol(tt)))
    }
  }
  
  ## Reordering so highest agreement is on diagonal
  ord <- order(apply(tt, 2, which.max))
  tt <- tt[, ord]
  colnames(tt) <- 1:ncol(tt)
  rownames(tt) <- 1:nrow(tt)
  ## In case of imperfect groupings
  # F1 score / F-measure
  f1 <- confusionMatrix(tt, mode = 'prec_recall')$byClass
  if(is.matrix(f1)) f1 <- f1[,"F1"]
  else f1 <- f1["F1"]
  
  f1[is.nan(f1)] <- 0
  f1[is.na(f1)] <- 0
  mean(f1)
}

# Combo Flag Lapprox ------------------------------------------------------

comboU <- function(L,k){
  vec <- eigen(L)$vector 
  ind <- (ncol(vec)-k+1):ncol(vec)
  U <- vec[,ind]
  
  tcrossprod(U)
}

flagLapproxCombo <- function(Llist, Klist, ngrp){
  # ngrp <- max(Klist)
  approxList <- mapply(Lapprox, Llist, Klist, SIMPLIFY=FALSE)
  approxList <- lapply(approxList, function(aL) aL/length(approxList))
  
  UL <- mapply(comboU, Llist, Klist, SIMPLIFY=FALSE)
  
  LR <- Reduce("+", approxList)
  Uks <- Reduce("+", UL)
  
  eig <- eigen(LR+Uks)$vectors[,1:ngrp]
  km <- kmeans(eig, centers=ngrp, nstart=10)
  
  list(eLR=eig, km=km)
}

# SPECTRUM functions ------------------------------------------------------

specSil <- function(sp){
  vecs <- sp$eigensystem$vectors[,1:sp$K]
  ss <- silhouette(sp$assignments, dist(vecs))
  mean(ss[, 3])
}

specARI <- function(sp, grps){
  adjustedRandIndex(sp$assignments, grps)
}

specFMeas <- function(sp, grps, k){
  tt <- table(sp$assignments, grps)
  
  ## Making confusion matrix square
  while(nrow(tt) != ncol(tt)){
    if(nrow(tt) > ncol(tt)){
      tt <- cbind(tt, rep(0, nrow(tt)))
    }
    if(ncol(tt) > nrow(tt)){
      tt <- rbind(tt, rep(0, ncol(tt)))
    }
  }
  
  ## Reordering so highest agreement is on diagonal
  ord <- order(apply(tt, 2, which.max))
  tt <- tt[, ord]
  colnames(tt) <- 1:ncol(tt)
  rownames(tt) <- 1:nrow(tt)
  ## In case of imperfect groupings
  # F1 score / F-measure
  f1 <- confusionMatrix(tt, mode = 'prec_recall')$byClass
  if(is.matrix(f1)) f1 <- f1[,"F1"]
  else f1 <- f1["F1"]
  
  f1[is.nan(f1)] <- 0
  f1[is.na(f1)] <- 0
  mean(f1)
}

# NEMO Functions ----------------------------------------------------------

nemoClust <- function(dd, nGrp){
  clust <- nemo.clustering(dd, num.clusters=nGrp)
  aff <- nemo.affinity.graph(dd)
  sim <- eigen(aff)$vectors[, 1:nGrp]
  
  list(clust=clust, sim=sim)
}

nemoSil <- function(ne){
  ss <- silhouette(ne$clust, dist(ne$sim))
  mean(ss[, 3])
}

nemoFMeas <- function(ne, grps){
  tt <- table(ne$clust, grps)
  
  ## Making confusion matrix square
  while(nrow(tt) != ncol(tt)){
    if(nrow(tt) > ncol(tt)){
      tt <- cbind(tt, rep(0, nrow(tt)))
    }
    if(ncol(tt) > nrow(tt)){
      tt <- rbind(tt, rep(0, ncol(tt)))
    }
  }
  
  ## Reordering so highest agreement is on diagonal
  ord <- order(apply(tt, 2, which.max))
  tt <- tt[, ord]
  colnames(tt) <- 1:ncol(tt)
  rownames(tt) <- 1:nrow(tt)
  ## In case of imperfect groupings
  # F1 score / F-measure
  f1 <- confusionMatrix(tt, mode = 'prec_recall')$byClass
  if(is.matrix(f1)) f1 <- f1[,"F1"]
  else f1 <- f1["F1"]
  
  f1[is.nan(f1)] <- 0
  f1[is.na(f1)] <- 0
  mean(f1)
}

nemoARI <- function(ne, grps){
  adjustedRandIndex(ne$clust, grps)
}

ARI <- function(ll, grps, k, type = 'simple'){
  
  if(type == 'simple') ev <- eigen(ll)$vectors[,1:k]
  else if(type == 'flag') ev <- ll$u[,1:k]
  else stop("'type' should be 'simple' or 'flag' ")
  
  km <- kmeans(ev, centers = k, nstart = k)
  adjustedRandIndex(km$cluster, grps)
}

ARI.n <- function(ll, grps){
  adjustedRandIndex(ll$km$cluster, grps)
}

# SNF Functions -----------------------------------------------------------

snfClust <- function(dd, nGrp){
  dL <- lapply(dd, function(d){
    dis <- SNFtool::dist2(as.matrix(d),
                          as.matrix(d))^(1/2)
    SNFtool::affinityMatrix(dis)
  })
  
  W <- SNFtool::SNF(dL)
  sim <- eigen(W)$vectors[,1:nGrp]
  clust <- SNFtool::spectralClustering(W, K=nGrp)
  
  list(clust=clust, sim=sim)
}

# Data Generation Functions -----------------------------------------------


## Generate data in a circle
circle.data <- function(n, r){
  aa <- runif(n=n, max = 2*pi)
  ss <- diag(c(0.05,0.05))
  ee <- mvrnorm(n=n, mu=c(0,0), Sigma=ss)
  
  cbind(x = r*cos(aa), y = r*sin(aa)) + ee
}

clustStruct <- function(n, p, k, noiseDat=NULL, randNoise=2){
  if(any(n%%k!=0)) stop("n must be divisible by k.")
  
  mapply(function(kk, pp, rn){
    
    if(kk == 1){
      means <- 0
    } else means <- c(0, 2^( 1:(kk-1) ) )
    
    datL <- lapply(means, function(mm){
      mvrnorm(n/kk, rep(mm,pp), diag(pp))
    })
    
    dat <- do.call(rbind, datL)
    
    if(!is.null(noiseDat)){
      if(is.character(noiseDat)){
        S <- rn*diag(pp)
        noiseDat <- mvrnorm(n=n, mu=rep(0,pp), Sigma=S)
      } 
      dat <- dat+noiseDat
    } 
    
    dat
  }, kk=k, pp=p, rn=randNoise, SIMPLIFY=FALSE)
  
}

getGroups <- function(k){
  ttt <- numeric(0)
  for(i in 1:length(k)){
    tk <- n/k[i]
    ttt <- cbind(ttt, rep(1:k[i], each=tk))
  }
  as.data.frame(ttt)
}

trueGroups <- function(k){
  grps <- getGroups(k)
  uGrp <- unique(grps)
  
  grpN <- cbind(uGrp, Grps=1:nrow(uGrp))
  merge(grps, grpN)
}


# Estimating K ------------------------------------------------------------


## Looking accross eigen vectors (changing space)
slopeStat <- function(L, kMax=20, pp=1, nstart=10){
  vecs <- eigen(L)$vectors
  
  Ks <- vapply(2:kMax, function(K){
    ind <- (ncol(vecs)-K+1):ncol(vecs)
    silhouette_score(vecs[,ind], k=K, nstart=nstart)
  }, numeric(1))
  
  -diff(Ks)*Ks[-length(Ks)]^pp
}

slopeStatPlot <- function(Llist, kMax=20, pp=1, nstart=10){
  
  slpSts <- lapply(Llist, slopeStat, 
                   kMax=kMax, pp=pp, nstart=nstart)
  
  Ll <- length(Llist)
  vv <- paste("View", 1:Ll)
  sm <- lapply(slpSts, function(x) ifelse(x==max(x), T, F))
  
  data <- data.frame(y = unlist(slpSts),
                     x = rep(2:(kMax-1), Ll),
                     v = rep(vv, each=kMax-2),
                     sm = unlist(sm))
  
  ggplot(data, aes(x=x, y=y)) +
    geom_segment( aes(x=x, xend=x, y=0, yend=y), color="grey") +
    geom_point(aes(color=sm), size=1, show.legend = FALSE) +
    facet_wrap(~v, nrow = Ll, scales = 'free') +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) + 
    scale_colour_manual(values = c('dodgerblue3', 'firebrick2')) +
    xlab("Number of clusters") +
    ylab("Slope Statistic") +
    geom_text(data = subset(data, sm), 
              aes(x=x+0.5,y=y,label=x))
  
}

## Clustering space stays the same
slopeStatistic <- function(x, kMax=20, pp=1, nstart=10){
  Ks <- vapply(2:kMax, function(K){
    silhouette_score(x, k=K, nstart=nstart)
  }, numeric(1))
  
  -diff(Ks)*Ks[-length(Ks)]^pp
}

slopeStatisticPlot <- function(X, kMax=20, pp=1, nstart=10){
  
  slpSts <- slopeStatistic(X, kMax=kMax, pp=pp, nstart=nstart)
  sm <- ifelse(slpSts==max(slpSts), T, F)
  
  data <- data.frame(y = slpSts,
                     x = 2:(kMax-1),
                     sm = sm)
  
  ggplot(data, aes(x=x, y=y)) +
    geom_segment( aes(x=x, xend=x, y=0, yend=y), color="grey") +
    geom_point(aes(color=sm), size=1, show.legend = FALSE) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) + 
    scale_colour_manual(values = c('dodgerblue3', 'firebrick2')) +
    xlab("Number of clusters") +
    ylab("Slope Statistic") +
    geom_text(data = subset(data, sm), 
              aes(x=x+0.5,y=y,label=x))
  
}


eigenGap <- function(L, lap.type='shift'){
  ev <- eigen(L)$values
  ## Suggested k
  
  if(lap.type %in% c('sym', 'rw')){
    eg <- length(ev) - which.min(diff(ev))
  } else eg <- which.min(diff(ev[-1]))+1
  
  eg
}

svdGap <- function(X){
  X <- do.call(cbind, X)
  sv <- svd(X)$d
  which.min(diff(sv[-1]))+1
}

modalityK <- function(dat, kernel='Gaussian', lap.type='sym', kMax=20){
  
  rho <- median(dist(dat))
  if(kernel == 'Gaussian') A <- Gaussian_kernel(dat, rho = rho)
  
  if(kernel == 'Spectrum'){
    dtd <- as.data.frame( t(dat) )
    A <- Spectrum::CNN_kernel(dtd)
  }
  
  diag(A) <- 0 ## Adjacency matrices need diag of 0
  deg <- rowSums(A)
  gg <- graph_from_adjacency_matrix(A, mode='undirected', weighted=TRUE)
  
  if(lap.type == 'sym'){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = I - D^(-1/2) A D^(-1/2)
    L <- diag(nrow = nrow(A)) - diag(ds) %*% A %*% diag(ds)
  }
  
  ## Random walk Laplacian (best one)
  if(lap.type == 'rw'){
    ds <- ifelse(deg>0, 1/deg, 0)
    # L = I - D^-1 A
    L <- diag(nrow = nrow(A)) - diag(ds) %*% A
  } 
  
  vecs <- eigen(L)$vectors
  
  mod <- vapply(2:kMax, function(K){
    ind <- (ncol(vecs)-K+1):ncol(vecs)
    kmem <- kmeans(vecs[,ind], K, nstart=10)$cluster
    igraph::modularity(gg, kmem)
  }, numeric(1))
  
  which.max(mod)+1
}

SpectrumK <- function(dd, m=1){
  dt <- as.data.frame(t(dd))
  Spectrum(dt, method=m, silent=T, showres=F)$K
}
