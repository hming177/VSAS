library(Runuran)
lpre.vsas <- function(x, y, size=NULL, intercept=TRUE, normalize=TRUE, tau=0.5){
  x <- as.matrix(x)
  y <- as.matrix(y)
  np <- dim(x)
  n <- np[1]
  p <- np[2]
  
  if(intercept){
    meanx <- colMeans(x)
    mprody <- prod(y^(1/n))
  }else{
    meanx <- rep(0, p)
    mprody <- 1
  }
  x <- scale(x, meanx, FALSE)
  y <- y/mprody
  if(normalize){normx <- sqrt(colSums(x^2)) }else{ normx <- rep(1, p)}
  x <- scale(x, FALSE, normx)
  x2 <- x^2
  tx <- t(x)
  tx2 <- t(x2)
  yinv <- 1/y
  
  if(is.null(size)) {size <- seq(2, ceiling(n/log(n)), 1)}
  
  
  bk_fun <- function(x, y){
    x <- as.matrix(x)
    y <- as.vector(y)
    tx <- t(x)
    yinv <- 1/y
    par <- coef(lm(log(y)~x+0))
    b.est <- par+100
    for(k in 1:1000){
      f <- tx %*% (yinv * exp(x %*% par) - y*exp(-x %*% par)) #一阶导???
      tt <- yinv*exp(x %*% par)+y*exp(-x %*% par) 
      f_diff <- matrix(0, ncol(x), ncol(x))
      for (i in 1:n) {
        f_diff <- f_diff + as.matrix(x[i, ]) %*% tx[, i] * tt[i] 
      }
      b.est <- par-solve(f_diff)%*%f
      if(norm(b.est - par,"2")<=1e-4) break
      par <- b.est
    }
    return(b.est)
  }
  
  vsas <- function(ms){
    bk <- rep(0, p)
    gk <- tx2 %*% (yinv + y)
    dk <- -tx %*% (yinv - y)/gk
    Hk <- sqrt(gk)*abs(tau*bk+(1-tau)*dk)
    Ak <- which(Hk >= sort(Hk, decreasing =TRUE)[ms])
    maxiter <- min(2*ms, 20)
    for (k in 1:maxiter) {
      Ik <- setdiff(1:p, Ak)
      bk[Ik] <- 0
      bk[Ak] <- bk_fun(x[, Ak], y)
      dk[Ak] <- 0
      pe <- as.matrix(x[, Ak]) %*% bk[Ak]
      gk <- tx2 %*% (yinv * exp(pe) + y*exp(-pe))
      dk[Ik] <- -tx[Ik, ] %*% (yinv * exp(pe) - y*exp(-pe))/gk[Ik]
      Hk <- sqrt(gk)*abs(tau*bk+(1-tau)*dk)
      Anew <- which(Hk >= sort(Hk, decreasing =TRUE)[ms])
      if(setequal(Ak, Anew)){
        return(c(k, bk))
      }
      Ak <- Anew
    }
    Ak
    return(c(k, bk))
  }
  
  para_beta <- sapply(size, function(ms){vsas(ms)})
  iter <- para_beta[1, ]
  b <- scale(t(para_beta[-1, ]), FALSE, normx)
  beta <- rbind(log(mprody) - meanx %*% t(b), t(b))
  
  return(list(beta=beta, size=size, iter=iter))
}


DGPfun2 <- function(n=200, q=12, p=500, rho=0.5, error_type=1, R=3){
  ntrain <- n-100
  minbeta <- sqrt(2*log(p)/ntrain)
  beta <- rep(0,p)
  poi <- sort(sample(1:p, q))
  beta[poi] <- runif(q, minbeta, R*minbeta)
  X <- xx <- matrix(rnorm(n*p), n, p)
  X[, 1] <- xx[, 1]
  for (j in 2:(p-1)) {
    X[, j] <- xx[, j] + rho*(xx[, j-1]+xx[, j+1])
  }
  
  lpdf1 = function(x){ -x-1/x-log(x)+2}              
  gen1 = ars.new(logpdf=lpdf1, lb=0.0000001, ub=Inf)
  f = function(x){ (1/2)*x^2-(1/2)*0.25-log(x)+log(0.5) }
  uni = uniroot(f, c(0.51,2))$root
  epsilon = switch(error_type, ur(gen1, n), rlnorm(n), exp(runif(n, -2, 2)), runif(n, 0.5, uni))
  
  Y = exp(X%*%beta)*epsilon
  Xtrain <- X[1:ntrain, ]
  Xtest <-  X[(ntrain+1):n, ]
  ytrain <- Y[1:ntrain]
  ytest <- Y[(ntrain+1):n]
  return(list(xtrain=Xtrain, ytrain=ytrain, xtest=Xtest, ytest=ytest, beta=beta))
}


dat <- DGPfun2(n=300, q=10, p=5000, rho=0.2, error_type=3, R=3)
beta <- c(0, dat$beta)
xtrain <- dat$xtrain
ytrain <- dat$ytrain

beta[which(beta!=0)]
lpre.l0 <- lpre.vsas(xtrain, ytrain, size=10, intercept=FALSE)
lpre.l0$iter
lpre.l0$beta[which(beta!=0)]

norm(lpre.l0$beta-beta, "2")

