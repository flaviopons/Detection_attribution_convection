## Calcul de dimensions locales (d) et persistence (theta) sur un ensemble de points X.0, a
## partir d'une trajectoire de reference X.ref.
## Utilisation de l'estimateur de Sueveges ou de Ferro pour la persistence theta
## Inspiration du code matlab de Davide Faranda
## C'est a priori parallelisable
## Pascal Yiou (LSCE), Nov. 2020

Ferro <- function(X, rho){
  u <- quantile(X, probs = rho, na.rm=TRUE)
  Li <- which(X > u)
  Ti <- diff(Li)
  Nc <- length(Ti)
  theta <- 2 * (sum(Ti-1))^2 / (Nc * sum((Ti-1)*(Ti-2)))
  return(theta)
}

likelihood <- function(N, Nc, theta, sqSi){
  if(N == Nc){
    ll <- 2*Nc*log(theta) - theta*sqSi
  } else {
    ll <- (N-Nc)*log(1-theta) + 2*Nc*log(theta) - theta*sqSi
  }
  return(ll)
}

Sueveges <- function(X, rho){
  # browser()
  q <- 1 - rho
  u <- quantile(X, probs = rho, na.rm=TRUE)
  Li <- which(X > u)
  Ti <- diff(Li)
  Si = Ti-1
  sqSi = sum(q*Si)
  Nc = length(which(Si > 0))
  N = length(Ti)
  delta <- (sqSi + N + Nc)^2 - 8*Nc*sqSi
  theta <- (sqSi + N + Nc - sqrt(delta))/(2*sqSi)
  return(theta)
}

Caby <- function(X, rho, m){
  # browser()
  # The threshold is the p-quantile of the distribution
  u <- quantile(X, rho, na.rm=TRUE)
  # Number of exceedences of Y over the threshold u
  iexcess <- (X > u)
  N <- length(which(iexcess))
  p <- numeric(m)
  T <- length(X)
  # Computation of p_k
  for (k in seq.int(m)){
    p[k] <- sum(
      sapply(
        which(iexcess[1 : (T-k)]),
        function(i){
          pk <- iexcess[i+k]
          if(k > 1) pk <- pk & all(!iexcess[(i+1) : (i+k-1)])
          return(pk)
        }
      ),
      na.rm = TRUE
    )
    p[k] <- p[k]/N;
  }
  # Computation of theta
  # print(p)
  theta <- 1 - sum(p)
  return(theta)
}

Caby2 <- function(X, rho, m){
  # The threshold is the p-quantile of the distribution
  u <- quantile(X, rho, na.rm=TRUE)
  # Number of exceedences of Y over the threshold u
  iexcess <- (X > u)
  N <- length(which(iexcess))
  p <- numeric(m)
  T <- length(X)
  # Computation of p_0
  for (i in 1:(length(X)-1)){
    if(X[i]>u & X[i+1]>u){
      p[1] <- p[1]+1;
    }
  }
  p[1] <- p[1]/N;
  # Computation of p_k
  if(m > 1){
    for (k in 2:m){
      for(i in 1:(length(X)-k)){
        if(X[i] > u){ 
          if (max(X[(i+1):(i+k-1)])<u & X[i+k]>u){
            p[k] <- p[k]+1
          }
        }
      }
      p[k] <- p[k]/N;
    }
  }
  # print(p)
  # Computation of theta
  theta <- 1 - sum(p)
  return(theta)
}

localdim <- function(X, rho){
  # X <- X[!is.infinite(X)]
  # u <- quantile(X, probs = rho, na.rm=TRUE)
  # mean_excess <- mean(X[X > u], na.rm=TRUE) - u
  # dim <- 1 / mean_excess
  # return(dim)
 X <- X[!is.infinite(X)]
 u <- quantile(X, probs = rho, na.rm=TRUE)
 Y <- as.matrix(X[X > u]-u)
 mean_excess <- mean(Y, na.rm=TRUE)
 var_excess <- var(Y, na.rm = TRUE)
 dim <- 1/(0.5*mean_excess*(1+mean_excess^2/var_excess))
 return(dim)
}

dtheta_onepoint <- function(dX, rho = 0.98, method = "Ferro", method.args = list()){
    # dim(X0) <- c(1, length(X0))
    # if(is.vector(Xref)){
    #    dim(Xref) <- c(length(Xref), 1)
    # }
    # dX <- -log(pracma::pdist2(Xref, X0))
    # dX <- -log(colMeans((t(Xref) - X0)^2))
    # dX[is.infinite(dX)] <- NA
    dim <- localdim(dX, rho)
    theta <- switch(method,
        Ferro = Ferro(dX, rho),
        Sueveges = Sueveges(dX, rho),
        Caby = do.call(Caby, c(X = list(dX), rho = list(rho), method.args)),
        Caby2 = do.call(Caby2, c(X = list(dX), rho = list(rho), method.args))
    )
    dtheta <- c(dim, theta)
    invisible(dtheta)
}


dtheta_allpoints <- function(X0, Xref = X0, rho = 0.98, method = "Ferro", method.args = list(), ...){
    q <- 1 - rho
    if(is.vector(Xref) | is.vector(X0)){
        dim(X0) <- c(length(X0), 1)
        dim(Xref) <- c(length(Xref), 1)
    }
    # wordspace.openmp(threads = 1) # number of threads for parallel computing
    distmat <- wordspace::dist.matrix(X0, Xref, method = "euclidean")
    distmat <- -log(distmat)
    dtheta <- matrix(NA, nrow = 2, ncol = nrow(X0))
    for (irow in seq.int(nrow(X0))) {
      row <- distmat[irow, -irow]
      dtheta[, irow] <- dtheta_onepoint(row, rho, method, method.args)
    }
    invisible(dtheta)
}
