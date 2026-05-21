library(extraDistr)

gbparma.fit <- function(y, ar, ma, mu.link, prec.link, names_phi, names_theta, 
                        names_beta, h1, X, X_hat, resid, diag){
  
  linkstr1 <- mu.link
  linkobj1 <- make.link(linkstr1)
  linkfun1 <- linkobj1$linkfun
  linkinv1 <- linkobj1$linkinv
  mu.eta1 <-  linkobj1$mu.eta

  linkstr2 <- prec.link
  linkobj2 <- make.link(linkstr2)
  linkfun2 <- linkobj2$linkfun
  linkinv2 <- linkobj2$linkinv
  mu.eta2 <-  linkobj2$mu.eta
  
  ynew <- linkfun1(y)
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,2,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  y_prev <- c(rep(NA,(n+h1)))
  prec_prev <- c(rep(NA,(n+h1)))
  
  # initial values
  if(any(is.na(ar)==F))
  {
    P <- matrix(rep(NA,(n-(m+1))*p1),ncol=p1)
    
    for(i in 1:(n-(m+1)))
    {
      P[i,] <- ynew[i+(m+1)-ar]
    }
    
    Z <- cbind(rep(1,(n-(m+1))),P)
  }else{
    Z <- as.matrix(rep(1,(n-(m+1))))
  }
  
  if(any(is.na(X)==F)){
    X_hat <- as.matrix(X_hat)
    X <- as.matrix(X)
    x <- cbind(as.matrix(Z),X[(m+2):n,])
    Y <- y[(m+2):n]
    Ynew <- linkfun1(Y)
    Ystar <- log(Y/(1+Y))
    ajuste <- lm.fit(x, Ynew)
    mqo <- c(ajuste$coef)
  }
  
  # BPAR --------------------------------------------------------------------
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==T)){
    
    loglik <- function(z){
      alpha1 <- z[1]
      phi <- z[2:(p1+1)]
      alpha2 <- z[p1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[p1+3]
      delta <- min(delta, -0.001)
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + (phi %*% ynew[i-ar]) 
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- linkinv1(eta1[seq1])
      prec <- linkinv2(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      a.mais.b <- a + b
      
      ll <- (a-1)*log(y1) - (a.mais.b)*log1p(y1) - lgamma(a) - lgamma(b) + lgamma(a.mais.b)
      
      return(sum(ll))
    }
    
    escore <- function(z){
      alpha1 <- z[1]
      phi <- z[2:(p1+1)]
      alpha2 <- z[p1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[p1+3]
      delta <- min(delta, -0.001)
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + (phi %*% ynew[i-ar]) 
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- linkinv1(eta1[seq1])
      prec <- linkinv2(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      n.menos.m <- n - (m+1)
      
      rP <- P
      vI <- rep(1, n.menos.m)
      c1 <- digamma(a+b)
      ystar <- log(y1/(1+y1))
      mustar <- digamma(a) - c1
      ydag <- mu*log(y1) - (1+mu)*log(1+y1)
      Delta <- c1 - digamma(b)
      mudag <- mu*(mustar - Delta/mu)
      f <- mu.eta1(eta1[seq1]) * ((prec+1) * (ystar-mustar))
      g <- mu.eta2(eta2[seq1]) * (ydag - mudag)
      
      Ualpha1 <- crossprod(vI, f)
      Uphi <<- crossprod(rP, f)
      Ualpha2 <- crossprod(vI, g)
      Udelta <<-  t(y[(m+1):(n-1)]/y[(m):(n-2)]) %*% g
      
      rval <- cbind(Ualpha1, t(Uphi), Ualpha2, Udelta)
      return(rval)
      
    }
    
    alpha1.ini <- log(mean(y))
    phi.ini <- rep(0,p1)
    alpha2.ini <- 3
    delta.ini <- -1
    
    reg <- c(alpha1.ini, phi.ini, alpha2.ini, delta.ini)
    
    names_par <- c("alpha1",names_phi,"alpha2","delta")
    
    opt <- optim(par = reg, fn = loglik, gr = escore,
                 method = "L-BFGS-B",
                 lower = c(-10,rep(-0.999,p1),0.01,-Inf),
                 upper = c(10,rep(0.999,p1),15,-0.01),
                 control = list(fnscale = -1, maxit = 1000))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")

      opt <- optim(par = reg, fn = loglik, #gr = escore,
                   method = "L-BFGS-B",
                   lower = c(-10,rep(-0.999,p1),0.01,-Inf),
                   upper = c(10,rep(0.999,p1),15,-0.01),
                   control = list(fnscale = -1, maxit = 1000))

      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+3)]
    names(coef) <- names_par
    z$coeff <- coef
    
    alpha1 <- coef[1]
    phi <- coef[2:(p1+1)]
    alpha2 <- coef[p1+2]
    delta <- coef[p1+3]
    
    z$alpha1 <- alpha1
    z$phi <- phi
    z$alpha2 <- alpha2
    z$delta <- delta
    
    errorhat <- rep(0,n)
    eta1hat <- eta2hat <- rep(NA,n)
    
    seq1 <- (m+2):n
    n.menos.m <- n - (m+1)
    
    for(i in seq1)
    {
      eta1hat[i] <- alpha1 + (phi%*%ynew[i-ar]) 
      eta2hat[i] <- alpha2 + delta*(y[i-1]/y[i-2])
    }
    
    muhat <- linkinv1(eta1hat[seq1])
    prechat <- linkinv2(eta2hat[seq1])
    y1 <- y[seq1]
    
    z$fitted <- ts(c(rep(NA,(m+1)),muhat),start=start(y),frequency=frequency(y))
    z$fitted_prec <- ts(c(rep(NA,(m+1)),prechat),start=start(y),frequency=frequency(y))
    z$eta1hat <- eta1hat
    z$errorhat <- errorhat
    z$eta2hat <- eta2hat
    z$mustarhat <- digamma(muhat*(1+prechat)) - digamma(muhat*(1+prechat)+prechat+2)
    
    rP <- P
    vI <- matrix(rep(1,n.menos.m),ncol=1)
    a <- muhat * (1 + prechat)
    b <- 2 + prechat
    c1 <- trigamma(a+b)
    mT <- diag(mu.eta1(eta1hat[seq1]))
    pT <- diag(mu.eta2(eta2hat[seq1]))
    D <- diag(as.vector((1+prechat)^2 * (trigamma(a) - c1))) %*% mT^2
    L <- pT %*% diag(as.vector(-(1+prechat) * (c1 + muhat*(c1-trigamma(a))))) %*% mT
    H <- diag(as.vector(muhat^2 * trigamma(a) - (1+muhat)^2 * c1 + trigamma(b))) %*% pT^2
    y.y2 <- (y[(m+1):(n-1)]/y[(m):(n-2)])
    
    Ka1a1 <- t(vI) %*% D %*% vI
    Ka1p <- t(vI) %*% D %*% rP
    Kpa1 <- t(Ka1p)
    Ka1a2 <- t(vI) %*% L %*% vI
    Ka2a1 <- t(Ka1a2)
    Ka1d <- t(vI) %*% L %*% y.y2
    Kda1 <- t(Ka1d)
    
    Kpp <- t(rP) %*% D %*% rP
    Kpa2 <- t(rP) %*% L %*% vI
    Ka2p <- t(Kpa2)
    Kpd <- t(rP) %*% L %*% y.y2
    Kdp <- t(Kpd)
    
    Ka2a2 <- t(vI) %*% H %*% vI
    Ka2d <- t(vI) %*% H %*% y.y2
    Kda2 <- t(Ka2d)
    
    Kdd <- t(y.y2) %*% H %*% y.y2
    
    K <- rbind(
      cbind(Ka1a1,Ka1p,Ka1a2,Ka1d),
      cbind(Kpa1 ,Kpp ,Kpa2 ,Kpd ),
      cbind(Ka2a1,Ka2p,Ka2a2,Ka2d),
      cbind(Kda1 ,Kdp ,Kda2 ,Kdd )
    )
    
    z$K <- K
    
    # #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    prec_prev[1:n] <- z$fitted_prec
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha1 + (phi%*%ynew_prev[n+i-ar])
      y_prev[n+i] <- linkinv1(ynew_prev[n+i])
      prec_prev[n+i] <- linkinv2(alpha2 + delta%*%(y_prev[n+i-1]/y_prev[n+i-2]))
    }
    
    z$serie <- y
    z$forecast <- y_prev[(n+1):(n+h1)]
    
  }
  
  # BPMA --------------------------------------------------------------------
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==T)){
    
    loglik <- function(z){
      alpha1 <- z[1]
      theta <- z[2:(q1+1)]
      alpha2 <- z[q1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[q1+3]
      delta <- min(delta, -0.001)
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + (theta %*% error[i-ma])
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- linkinv1(eta1[seq1])
      prec <- linkinv2(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      a.mais.b <- a + b
      
      ll <- (a-1)*log(y1) - (a.mais.b)*log1p(y1) - lgamma(a) - lgamma(b) + lgamma(a.mais.b)
      
      return(sum(ll))
    }
    
    escore <- function(z){
      alpha1 <- z[1]
      theta <- z[2:(q1+1)]
      alpha2 <- z[q1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[q1+3]
      delta <- min(delta, -0.001)
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + (theta %*% error[i-ma])
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- linkinv1(eta1[seq1])
      prec <- linkinv2(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      n.menos.m <- n - (m+1)
      
      R <- matrix(rep(NA, n.menos.m*q1), ncol=q1)
      for(i in 1:n.menos.m)
      {
        R[i,] <- error[i+(m+1)-ma]
      }
      
      ### FB recorrences
      deta.dalpha <- rep(0, n)
      deta.dtheta <- matrix(0, ncol=q1, nrow=n)
      
      for(i in seq1)
      {
        i.menos.ma <- i-ma
        deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i.menos.ma]
        deta.dtheta[i,] <- R[(i-(m+1)),] - theta %*% deta.dtheta[i.menos.ma,]
      }
      
      v <- deta.dalpha[seq1]
      rR <- deta.dtheta[seq1,]
      vI <- rep(1, n.menos.m)
      c1 <- digamma(a+b)
      ystar <- log(y1/(1+y1))
      mustar <- digamma(a) - c1
      ydag <- mu*log(y1) - (1+mu)*log(1+y1)
      Delta <- c1 - digamma(b)
      mudag <- mu*(mustar - Delta/mu)
      f <- mu.eta1(eta1[seq1]) * ((prec+1) * (ystar-mustar))
      g <- mu.eta2(eta2[seq1]) * (ydag - mudag)
      
      Ualpha1 <- crossprod(v, f)
      Utheta <<- crossprod(rR, f)
      Ualpha2 <- crossprod(vI, g)
      Udelta <<-  t(y[(m+1):(n-1)]/y[(m):(n-2)]) %*% g
      
      rval <- cbind(Ualpha1, t(Utheta), Ualpha2, Udelta)
      return(rval)
    }
    
    alpha1.ini <- log(mean(y))
    theta.ini <- rep(0,q1)
    alpha2.ini <- 3
    delta.ini <- -1
    
    reg <- c(alpha1.ini, theta.ini, alpha2.ini, delta.ini)
    
    names_par <- c("alpha1",names_theta,"alpha2","delta")
    
    opt <- optim(par = reg, fn = loglik, gr = escore,
                 method = "L-BFGS-B",
                 lower = c(-10,rep(-0.999,q1),0.01,-Inf),
                 upper = c(10,rep(0.999,q1),15,-0.01),
                 control = list(fnscale = -1, maxit = 1000))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      
      opt <- optim(par = reg, fn = loglik, #gr = escore,
                   method = "L-BFGS-B",
                   lower = c(-10,rep(-0.999,q1),0.01,-Inf),
                   upper = c(10,rep(0.999,q1),15,-0.01),
                   control = list(fnscale = -1, maxit = 1000))
      
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+3)]
    names(coef) <- names_par
    z$coeff <- coef
    
    alpha1 <- coef[1]
    theta <- coef[2:(q1+1)]
    alpha2 <- coef[q1+2]
    delta <- coef[q1+3]
    
    z$alpha1 <- alpha1
    z$theta <- theta
    z$alpha2 <- alpha2
    z$delta <- delta
    
    errorhat <- rep(0,n)
    eta1hat <- eta2hat <- rep(NA,n)
    
    seq1 <- (m+2):n
    n.menos.m <- n - (m+1)
    
    for(i in seq1)
    {
      eta1hat[i] <- alpha1 + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i] - eta1hat[i]
      eta2hat[i] <- alpha2 + delta*(y[i-1]/y[i-2])
    }
    
    muhat <- linkinv1(eta1hat[seq1])
    prechat <- linkinv2(eta2hat[seq1])
    y1 <- y[seq1]
    
    z$fitted <- ts(c(rep(NA,(m+1)),muhat),start=start(y),frequency=frequency(y))
    z$fitted_prec <- ts(c(rep(NA,(m+1)),prechat),start=start(y),frequency=frequency(y))
    z$eta1hat <- eta1hat
    z$errorhat <- errorhat
    z$eta2hat <- eta2hat
    z$mustarhat <- digamma(muhat*(1+prechat)) - digamma(muhat*(1+prechat)+prechat+2)
    
    R <- matrix(rep(NA, n.menos.m*q1), ncol=q1)
    for(i in 1:n.menos.m)
    {
      R[i,] <- errorhat[i+(m+1)-ma]
    }
    
    # ##FB  recorrences
    deta.dalpha <- rep(0,n)
    deta.dtheta <- matrix(0, ncol=q1,nrow=n)
    
    for(i in seq1)
    {
      i.menos.ma <- i-ma
      deta.dalpha[i] <- 1 - theta%*%deta.dalpha[i.menos.ma]
      deta.dtheta[i,] <- R[(i-(m+1)),] - theta%*%deta.dtheta[i.menos.ma,]
    }
    
    v <- matrix(as.vector(deta.dalpha[seq1]),ncol=1)
    rR <- deta.dtheta[seq1,]
    vI <- matrix(rep(1,n.menos.m),ncol=1)
    a <- muhat * (1 + prechat)
    b <- 2 + prechat
    c1 <- trigamma(a+b)
    mT <- diag(mu.eta1(eta1hat[seq1]))
    pT <- diag(mu.eta2(eta2hat[seq1]))
    D <- diag(as.vector((1+prechat)^2 * (trigamma(a) - c1))) %*% mT^2
    L <- pT %*% diag(as.vector(-(1+prechat) * (c1 + muhat*(c1-trigamma(a))))) %*% mT
    H <- diag(as.vector(muhat^2 * trigamma(a) - (1+muhat)^2 * c1 + trigamma(b))) %*% pT^2
    y.y2 <- (y[(m+1):(n-1)]/y[(m):(n-2)])
    
    Ka1a1 <- t(v) %*% D %*% v
    Ka1t <- t(v) %*% D %*% rR
    Kta1 <- t(Ka1t)
    Ka1a2 <- t(v) %*% L %*% vI
    Ka2a1 <- t(Ka1a2)
    Ka1d <- t(v) %*% L %*% y.y2
    Kda1 <- t(Ka1d)
    
    Ktt <- t(rR) %*% D %*% rR
    Kta2 <- t(rR) %*% L %*% vI
    Ka2t <- t(Kta2)
    Ktd <- t(rR) %*% L %*% y.y2
    Kdt <- t(Ktd)
    
    Ka2a2 <- t(vI) %*% H %*% vI
    Ka2d <- t(vI) %*% H %*% y.y2
    Kda2 <- t(Ka2d)
    
    Kdd <- t(y.y2) %*% H %*% y.y2
    
    K <- rbind(
      cbind(Ka1a1,Ka1t,Ka1a2,Ka1d),
      cbind(Kta1 ,Ktt ,Kta2 ,Ktd ),
      cbind(Ka2a1,Ka2t,Ka2a2,Ka2d),
      cbind(Kda1 ,Kdt ,Kda2 ,Kdd )
    )
    
    z$K <- K
    
    # #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    prec_prev[1:n] <- z$fitted_prec
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha1 + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv1(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
      
      prec_prev[n+i] <- linkinv2(alpha2 + delta%*%(y_prev[n+i-1]/y_prev[n+i-2]))
    }
    
    z$serie <- y
    # z$bparma <- names_par
    z$forecast <- y_prev[(n+1):(n+h1)]
    
  }
  
  # BPARMA ------------------------------------------------------------------
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T)){
    
    loglik <- function(z){
      p1.mais.q1 <- p1+q1
      alpha1 <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1.mais.q1+1)]
      alpha2 <- z[p1.mais.q1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[p1.mais.q1+3]
      delta <- min(delta, -0.001)
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + (phi %*% ynew[i-ar]) + (theta %*% error[i-ma])
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- linkinv1(eta1[seq1])
      prec <- linkinv2(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      a.mais.b <- a + b
      
      ll <- (a-1)*log(y1) - (a.mais.b)*log1p(y1) - lgamma(a) - lgamma(b) + lgamma(a.mais.b)
      
      return(sum(ll))
    }
    
    escore <- function(z){
      p1.mais.q1 <- p1+q1
      alpha1 <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1.mais.q1+1)]
      alpha2 <- z[p1.mais.q1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[p1.mais.q1+3]
      delta <- min(delta, -0.001)
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + (phi %*% ynew[i-ar]) + (theta %*% error[i-ma])
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- linkinv1(eta1[seq1])
      prec <- linkinv2(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      n.menos.m <- n - (m+1)
      
      R <- matrix(rep(NA, n.menos.m*q1), ncol=q1)
      for(i in 1:n.menos.m)
      {
        R[i,] <- error[i+(m+1)-ma]
      }
      
      ### FB recorrences
      deta.dalpha <- rep(0, n)
      deta.dphi <- matrix(0, ncol=p1, nrow=n)
      deta.dtheta <- matrix(0, ncol=q1, nrow=n)
      
      for(i in seq1)
      {
        i.menos.ma <- i-ma
        deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i.menos.ma]
        deta.dphi[i,] <- P[(i-(m+1)),] - theta %*% deta.dphi[i.menos.ma,]
        deta.dtheta[i,] <- R[(i-(m+1)),] - theta %*% deta.dtheta[i.menos.ma,]
      }
      
      v <- deta.dalpha[seq1]
      rP <- deta.dphi[seq1,]
      rR <- deta.dtheta[seq1,]
      vI <- rep(1, n.menos.m)
      c1 <- digamma(a+b)
      ystar <- log(y1/(1+y1))
      mustar <- digamma(a) - c1
      ydag <- mu*log(y1) - (1+mu)*log(1+y1)
      Delta <- c1 - digamma(b)
      mudag <- mu*(mustar - Delta/mu)
      f <- mu.eta1(eta1[seq1]) * ((prec+1) * (ystar-mustar))
      g <- mu.eta2(eta2[seq1]) * (ydag - mudag)
      
      Ualpha1 <- crossprod(v, f)
      Uphi <<- crossprod(rP, f)
      Utheta <<- crossprod(rR, f)
      Ualpha2 <- crossprod(vI, g)
      Udelta <<-  t(y[(m+1):(n-1)]/y[(m):(n-2)]) %*% g
      
      rval <- cbind(Ualpha1, t(Uphi), t(Utheta), Ualpha2, Udelta)
      return(rval)
    }
    
    alpha1.ini <- log(mean(y))
    phi.ini <- rep(0,p1)
    theta.ini <- rep(0,q1)
    alpha2.ini <- 3
    delta.ini <- -1
    
    reg <- c(alpha1.ini, phi.ini, theta.ini, alpha2.ini, delta.ini)
    
    names_par <- c("alpha1",names_phi,names_theta,"alpha2","delta")
    
    opt <- optim(par = reg, fn = loglik, gr = escore,
                 method = "L-BFGS-B",
                 lower = c(-10,rep(-0.999,p1),rep(-0.999,q1),0.01,-Inf),
                 upper = c(10,rep(0.999,p1),rep(0.999,q1),15,-0.01),
                 control = list(fnscale = -1, maxit = 1000))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      
      opt <- optim(par = reg, fn = loglik, #gr = escore,
                   method = "L-BFGS-B",
                   lower = c(-10,rep(-0.999,p1),rep(-0.999,q1),0.01,-Inf),
                   upper = c(10,rep(0.999,p1),rep(0.999,q1),15,-0.01),
                   control = list(fnscale = -1, maxit = 1000))
      
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+3)]
    names(coef) <- names_par
    z$coeff <- coef
    
    alpha1 <- coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    alpha2 <- coef[p1+q1+2]
    delta <- coef[p1+q1+3]
    
    z$alpha1 <- alpha1
    z$phi <- phi
    z$theta <- theta
    z$alpha2 <- alpha2
    z$delta <- delta
    
    errorhat <- rep(0,n)
    eta1hat <- eta2hat <- rep(NA,n)
    
    seq1 <- (m+2):n
    n.menos.m <- n - (m+1)
    
    for(i in seq1)
    {
      eta1hat[i] <- alpha1 + (phi%*%ynew[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i] - eta1hat[i]
      eta2hat[i] <- alpha2 + delta*(y[i-1]/y[i-2])
    }
    
    muhat <- linkinv1(eta1hat[seq1])
    prechat <- linkinv2(eta2hat[seq1])
    y1 <- y[seq1]
    
    z$fitted <- ts(c(rep(NA,(m+1)),muhat),start=start(y),frequency=frequency(y))
    z$fitted_prec <- ts(c(rep(NA,(m+1)),prechat),start=start(y),frequency=frequency(y))
    z$eta1hat <- eta1hat
    z$errorhat <- errorhat
    z$eta2hat <- eta2hat
    z$mustarhat <- digamma(muhat*(1+prechat)) - digamma(muhat*(1+prechat)+prechat+2)
    
    R <- matrix(rep(NA, n.menos.m*q1), ncol=q1)
    for(i in 1:n.menos.m)
    {
      R[i,] <- errorhat[i+(m+1)-ma]
    }
    
    # ##FB  recorrences
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta <- matrix(0, ncol=q1,nrow=n)
    
    for(i in seq1)
    {
      i.menos.ma <- i-ma
      deta.dalpha[i] <- 1 - theta%*%deta.dalpha[i.menos.ma]
      deta.dphi[i,] <- P[(i-(m+1)),] - theta%*%deta.dphi[i.menos.ma,]
      deta.dtheta[i,] <- R[(i-(m+1)),] - theta%*%deta.dtheta[i.menos.ma,]
    }
    
    v <- matrix(as.vector(deta.dalpha[seq1]),ncol=1)
    rP <- deta.dphi[seq1,]
    rR <- deta.dtheta[seq1,]
    vI <- matrix(rep(1,n.menos.m),ncol=1)
    a <- muhat * (1 + prechat)
    b <- 2 + prechat
    c1 <- trigamma(a+b)
    mT <- diag(mu.eta1(eta1hat[seq1]))
    pT <- diag(mu.eta2(eta2hat[seq1]))
    D <- diag(as.vector((1+prechat)^2 * (trigamma(a) - c1))) %*% mT^2
    L <- pT %*% diag(as.vector(-(1+prechat) * (c1 + muhat*(c1-trigamma(a))))) %*% mT
    H <- diag(as.vector(muhat^2 * trigamma(a) - (1+muhat)^2 * c1 + trigamma(b))) %*% pT^2
    y.y2 <- (y[(m+1):(n-1)]/y[(m):(n-2)])
    
    Ka1a1 <- t(v) %*% D %*% v
    Ka1p <- t(v) %*% D %*% rP
    Kpa1 <- t(Ka1p)
    Ka1t <- t(v) %*% D %*% rR
    Kta1 <- t(Ka1t)
    Ka1a2 <- t(v) %*% L %*% vI
    Ka2a1 <- t(Ka1a2)
    Ka1d <- t(v) %*% L %*% y.y2
    Kda1 <- t(Ka1d)
    
    Kpp <- t(rP) %*% D %*% rP
    Kpt <- t(rP) %*% D %*% rR
    Ktp <- t(Kpt)
    Kpa2 <- t(rP) %*% L %*% vI
    Ka2p <- t(Kpa2)
    Kpd <- t(rP) %*% L %*% y.y2
    Kdp <- t(Kpd)
    
    Ktt <- t(rR) %*% D %*% rR
    Kta2 <- t(rR) %*% L %*% vI
    Ka2t <- t(Kta2)
    Ktd <- t(rR) %*% L %*% y.y2
    Kdt <- t(Ktd)
    
    Ka2a2 <- t(vI) %*% H %*% vI
    Ka2d <- t(vI) %*% H %*% y.y2
    Kda2 <- t(Ka2d)
    
    Kdd <- t(y.y2) %*% H %*% y.y2
    
    K <- rbind(
      cbind(Ka1a1,Ka1p,Ka1t,Ka1a2,Ka1d),
      cbind(Kpa1 ,Kpp ,Kpt ,Kpa2 ,Kpd ),
      cbind(Kta1 ,Ktp ,Ktt ,Kta2 ,Ktd ),
      cbind(Ka2a1,Ka2p,Ka2t,Ka2a2,Ka2d),
      cbind(Kda1 ,Kdp ,Kdt ,Kda2 ,Kdd )
    )
    
    z$K <- K
    
    # #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    prec_prev[1:n] <- z$fitted_prec
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha1 + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv1(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
      
      prec_prev[n+i] <- linkinv2(alpha2 + delta%*%(y_prev[n+i-1]/y_prev[n+i-2]))
    }
    
    z$serie <- y
    # z$bparma <- names_par
    z$forecast <- y_prev[(n+1):(n+h1)]
    
  }
  
  # BPARX -------------------------------------------------------------------
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==F)){
    
    loglik <- function(z){
      alpha1 <- z[1]
      phi <- z[2:(p1+1)]
      alpha2 <- z[(p1+2)]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[(p1+3)]
      delta <- min(delta, -0.001)
      beta <- z[(p1+4):length(z)]
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + X[i,]%*%as.matrix(beta) + 
          (phi %*% (ynew[i-ar]-X[i-ar,]%*%as.matrix(beta)))
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- exp(eta1[seq1])
      prec <- exp(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      a.mais.b <- a + b
      
      ll <- (a-1)*log(y1) - (a.mais.b)*log1p(y1) - lgamma(a) - lgamma(b) + lgamma(a.mais.b)
      
      return(sum(ll))
    }
    
    escore <- function(z){
      alpha1 <- z[1]
      phi <- z[2:(p1+1)]
      alpha2 <- z[(p1+2)]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[(p1+3)]
      delta <- min(delta, -0.001)
      beta <- z[(p1+4):length(z)]
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + X[i,]%*%as.matrix(beta) + 
          (phi %*% (ynew[i-ar]-X[i-ar,]%*%as.matrix(beta)))
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- exp(eta1[seq1])
      prec <- exp(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      n.menos.m <- n - (m+1)
      
      P <- matrix(rep(NA, n.menos.m*p1), ncol=p1)
      for(i in 1:n.menos.m)
      {
        P[i,] <- ynew[i+(m+1)-ar] - X[i+(m+1)-ar,]%*%as.matrix(beta)
      }
      
      k1<- length(beta)
      M <- matrix(rep(NA, n.menos.m*k1),ncol=k1)
      for(i in 1:n.menos.m)
      {
        for(j in 1:k1)
          M[i,j] <- X[i+(m+1),j]-sum(phi*X[i+(m+1)-ar,j])
      }
      
      v <- matrix(rep(1,(n-(m+1))),ncol=1)
      rP <- P
      rM <- M
      vI <- rep(1, n-(m+1))
      c1 <- digamma(a+b)
      ystar <- log(y1/(1+y1))
      mustar <- digamma(a) - c1
      ydag <- mu*log(y1) - (1+mu)*log(1+y1)
      Delta <- c1 - digamma(b)
      mudag <- mu*(mustar - Delta/mu)
      f <- exp(eta1[seq1]) * ((prec+1) * (ystar-mustar))
      g <- exp(eta2[seq1]) * (ydag - mudag)
      
      Ualpha1 <- crossprod(v, f)
      Uphi <<- crossprod(rP, f)
      Ualpha2 <- crossprod(vI, g)
      Udelta <<-  t(y[(m+1):(n-1)]/y[(m):(n-2)]) %*% g
      Ubeta <<- crossprod(rM, f)
      
      rval <- cbind(Ualpha1, t(Uphi), Ualpha2, Udelta, t(Ubeta))
      return(rval)
    }
    
    alpha1.ini <- log(mean(y))
    phi.ini <- rep(0,p1)
    alpha2.ini <- 3
    delta.ini <- -1
    beta.ini <- rep(0,ncol(X))
    
    reg <- c(alpha1.ini, phi.ini, alpha2.ini, delta.ini, beta.ini)
    
    names_par <- c("alpha1",names_phi,"alpha2","delta",names_beta)
    
    opt <- optim(par = reg, fn = loglik, gr = escore,
                 method = "L-BFGS-B",
                 lower = c(-10,rep(-0.999,p1),0.01,-Inf, rep(-10, ncol(X))),
                 upper = c(10,rep(0.999,p1),15,-0.01, rep(10, ncol(X))),
                 control = list(fnscale = -1, maxit = 1000))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")

      opt <- optim(par = reg, fn = loglik, #gr = escore,
                   method = "L-BFGS-B",
                   lower = c(-10,rep(-0.999,p1),0.01,-Inf, rep(-10, ncol(X))),
                   upper = c(10,rep(0.999,p1),15,-0.01, rep(10, ncol(X))),
                   control = list(fnscale = -1, maxit = 1000))

      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+3+ncol(X))]
    names(coef) <- names_par
    z$coeff <- coef
    
    alpha1 <- coef[1]
    phi <- coef[2:(p1+1)]
    alpha2 <- coef[p1+2]
    delta <- coef[p1+3]
    beta <- coef[(p1+4):length(coef)]
    
    z$alpha1 <- alpha1
    z$phi <- phi
    z$alpha2 <- alpha2
    z$delta <- delta
    z$beta <- beta
    
    errorhat <- rep(0,n)
    eta1hat <- eta2hat <- rep(NA,n)
    
    seq1 <- (m+2):n
    n.menos.m <- n - (m+1)
    
    for(i in seq1)
    {
      eta1hat[i] <- alpha1 + X[i,]%*%as.matrix(beta) + 
        (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta)))
      errorhat[i] <- ynew[i] - eta1hat[i]
      eta2hat[i] <- alpha2 + delta*(y[i-1]/y[i-2])
    }
    
    muhat <- exp(eta1hat[seq1])
    prechat <- exp(eta2hat[seq1])
    y1 <- y[seq1]
    
    z$fitted <- ts(c(rep(NA,(m+1)),muhat),start=start(y),frequency=frequency(y))
    z$fitted_prec <- ts(c(rep(NA,(m+1)),prechat),start=start(y),frequency=frequency(y))
    z$eta1hat <- eta1hat
    z$errorhat <- errorhat
    z$eta2hat <- eta2hat
    z$mustarhat <- digamma(muhat*(1+prechat)) - digamma(muhat*(1+prechat)+prechat+2)
    
    P <- matrix(rep(NA, n.menos.m*p1), ncol=p1)
    for(i in 1:n.menos.m)
    {
      P[i,] <- ynew[i+(m+1)-ar] - X[i+(m+1)-ar,]%*%as.matrix(beta)
    }
    
    k1<- length(beta)
    M <- matrix(rep(NA, n.menos.m*k1),ncol=k1)
    for(i in 1:n.menos.m)
    {
      for(j in 1:k1)
        M[i,j] <- X[i+(m+1),j]-sum(phi*X[i+(m+1)-ar,j])
    }
    
    v <- matrix(rep(1,n.menos.m),ncol=1)
    rP <- P
    rM <- M
    vI <- matrix(rep(1,n.menos.m),ncol=1)
    a <- muhat * (1 + prechat)
    b <- 2 + prechat
    c1 <- trigamma(a+b)
    mT <- diag(exp(eta1hat[seq1]))
    pT <- diag(exp(eta2hat[seq1]))
    D <- diag(as.vector((1+prechat)^2 * (trigamma(a) - c1))) %*% mT^2
    L <- pT %*% diag(as.vector(-(1+prechat) * (c1 + muhat*(c1-trigamma(a))))) %*% mT
    H <- diag(as.vector(muhat^2 * trigamma(a) - (1+muhat)^2 * c1 + trigamma(b))) %*% pT^2
    y.y2 <- (y[(m+1):(n-1)]/y[(m):(n-2)])
    
    Ka1a1 <- t(v) %*% D %*% v
    Ka1p <- t(v) %*% D %*% rP
    Kpa1 <- t(Ka1p)
    Ka1a2 <- t(v) %*% L %*% vI
    Ka2a1 <- t(Ka1a2)
    Ka1d <- t(v) %*% L %*% y.y2
    Kda1 <- t(Ka1d)
    
    Kpp <- t(rP) %*% D %*% rP
    Kpa2 <- t(rP) %*% L %*% vI
    Ka2p <- t(Kpa2)
    Kpd <- t(rP) %*% L %*% y.y2
    Kdp <- t(Kpd)
    
    Ka2a2 <- t(vI) %*% H %*% vI
    Ka2d <- t(vI) %*% H %*% y.y2
    Kda2 <- t(Ka2d)
    
    Kdd <- t(y.y2) %*% H %*% y.y2
    
    Ka1b <- t(v) %*% D %*% rM
    Kba1 <- t(Ka1b)
    Kpb <- t(rP) %*% D %*% rM
    Kbp <- t(Kpb)
    Ka2b <- t(vI) %*% L %*% rM
    Kba2 <- t(Ka2b)
    Kdb <- t(y.y2) %*% L %*% rM
    Kbd <- t(Kdb)
    Kbb <- t(rM) %*% D %*% rM
    
    K <- rbind(
      cbind(Ka1a1,Ka1p,Ka1a2,Ka1d,Ka1b),
      cbind(Kpa1 ,Kpp ,Kpa2 ,Kpd ,Kpb ),
      cbind(Ka2a1,Ka2p,Ka2a2,Ka2d,Ka2b),
      cbind(Kda1 ,Kdp ,Kda2 ,Kdd ,Kdb ),
      cbind(Kba1 ,Kbp ,Kba2 ,Kbd ,Kbb )
    )
    
    z$K <- K
    
    # #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    prec_prev[1:n] <- z$fitted_prec
    
    X_prev<- rbind(X,X_hat)
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha1 + X_prev[n+i,]%*%as.matrix(beta) +
        (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta))) 
      y_prev[n+i] <- exp(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
      
      prec_prev[n+i] <- exp(alpha2 + delta%*%(y_prev[n+i-1]/y_prev[n+i-2]))
    }
    
    z$serie <- y
    z$forecast <- y_prev[(n+1):(n+h1)]
    z$forecast_prev <- prec_prev[(n+1):(n+h1)]
    
  }
  
  # BPMAX -------------------------------------------------------------------
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==F)){
    
    loglik <- function(z){
      alpha1 <- z[1]
      theta <- z[(2):(q1+1)]
      alpha2 <- z[q1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[q1+3]
      delta <- min(delta, -0.001)
      beta <- z[(q1+4):length(z)]
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + X[i,]%*%as.matrix(beta) + (theta %*% error[i-ma])
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- exp(eta1[seq1])
      prec <- exp(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      a.mais.b <- a + b
      
      ll <- (a-1)*log(y1) - (a.mais.b)*log1p(y1) - lgamma(a) - lgamma(b) + lgamma(a.mais.b)
      
      return(sum(ll))
    }
    
    escore <- function(z){
      alpha1 <- z[1]
      theta <- z[(2):(q1+1)]
      alpha2 <- z[q1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[q1+3]
      delta <- min(delta, -0.001)
      beta <- z[(q1+4):length(z)]
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + X[i,]%*%as.matrix(beta) + (theta %*% error[i-ma])
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- exp(eta1[seq1])
      prec <- exp(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      n.menos.m <- n - (m+1)
      
      R <- matrix(rep(NA, n.menos.m*q1), ncol=q1)
      for(i in 1:n.menos.m)
      {
        R[i,] <- error[i+(m+1)-ma]
      }
      
      k1<- length(beta)
      M <- matrix(rep(NA, n.menos.m*k1),ncol=k1)
      for(i in 1:n.menos.m)
      {
        for(j in 1:k1)
          M[i,j] <- X[i+(m+1),j]
      }
      
      ### FB recorrences
      deta.dalpha <- rep(0, n)
      deta.dtheta <- matrix(0, ncol=q1, nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      for(i in seq1)
      {
        i.menos.ma <- i-ma
        deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i.menos.ma]
        deta.dtheta[i,] <- R[(i-(m+1)),] - theta %*% deta.dtheta[i.menos.ma,]
        deta.dbeta[i,]<- M[(i-(m+1)),] - theta%*%deta.dbeta[i.menos.ma,]
      }
      
      v <- deta.dalpha[seq1]
      rR <- deta.dtheta[seq1,]
      rM <- deta.dbeta[seq1,]
      vI <- rep(1, n-(m+1))
      c1 <- digamma(a+b)
      ystar <- log(y1/(1+y1))
      mustar <- digamma(a) - c1
      ydag <- mu*log(y1) - (1+mu)*log(1+y1)
      Delta <- c1 - digamma(b)
      mudag <- mu*(mustar - Delta/mu)
      f <- exp(eta1[seq1]) * ((prec+1) * (ystar-mustar))
      g <- exp(eta2[seq1]) * (ydag - mudag)
      
      Ualpha1 <- crossprod(v, f)
      Utheta <<- crossprod(rR, f)
      Ualpha2 <- crossprod(vI, g)
      Udelta <<-  t(y[(m+1):(n-1)]/y[(m):(n-2)]) %*% g
      Ubeta <<- crossprod(rM, f)
      
      rval <- cbind(Ualpha1, t(Utheta), Ualpha2, Udelta, t(Ubeta))
      return(rval)
    }
    
    alpha1.ini <- log(mean(y))
    theta.ini <- rep(0,q1)
    alpha2.ini <- 3
    delta.ini <- -1
    beta.ini <- mqo[(2):length(mqo)]
    
    reg <- c(alpha1.ini, theta.ini, alpha2.ini, delta.ini, beta.ini)
    
    names_par <- c("alpha1",names_theta,"alpha2","delta",names_beta)
    
    opt <- optim(par = reg, fn = loglik, gr = escore,
                 method = "L-BFGS-B", 
                 lower = c(-10,rep(-0.999,q1),0.01,-Inf, rep(-10, ncol(X))),
                 upper = c(10,rep(0.999,q1),15,-0.01, rep(10, ncol(X))),
                 control = list(fnscale = -1, maxit = 1000))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      
      opt <- optim(par = reg, fn = loglik, #gr = escore,
                   method = "L-BFGS-B", 
                   lower = c(-10,rep(-0.999,q1),0.01,-Inf, rep(-10, ncol(X))),
                   upper = c(10,rep(0.999,q1),15,-0.01, rep(10, ncol(X))),
                   control = list(fnscale = -1, maxit = 1000))
      
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+3+ncol(X))]
    names(coef) <- names_par
    z$coeff <- coef
    
    alpha1 <- coef[1]
    theta <- coef[(2):(q1+1)]
    alpha2 <- coef[q1+2]
    delta <- coef[q1+3]
    beta <- coef[(q1+4):length(coef)]
    
    z$alpha1 <- alpha1
    z$theta <- theta
    z$alpha2 <- alpha2
    z$delta <- delta
    z$beta <- beta
    
    errorhat <- rep(0,n)
    eta1hat <- eta2hat <- rep(NA,n)
    
    seq1 <- (m+2):n
    n.menos.m <- n - (m+1)
    
    for(i in seq1)
    {
      eta1hat[i] <- alpha1 + X[i,]%*%as.matrix(beta) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i] - eta1hat[i]
      eta2hat[i] <- alpha2 + delta*(y[i-1]/y[i-2])
    }
    
    muhat <- exp(eta1hat[seq1])
    prechat <- exp(eta2hat[seq1])
    y1 <- y[seq1]
    
    z$fitted <- ts(c(rep(NA,(m+1)),muhat),start=start(y),frequency=frequency(y))
    z$fitted_prec <- ts(c(rep(NA,(m+1)),prechat),start=start(y),frequency=frequency(y))
    z$eta1hat <- eta1hat
    z$errorhat <- errorhat
    z$eta2hat <- eta2hat
    z$mustarhat <- digamma(muhat*(1+prechat)) - digamma(muhat*(1+prechat)+prechat+2)
    
    R <- matrix(rep(NA, n.menos.m*q1), ncol=q1)
    for(i in 1:n.menos.m)
    {
      R[i,] <- errorhat[i+(m+1)-ma]
    }
    
    k1<- length(beta)
    M <- matrix(rep(NA, n.menos.m*k1),ncol=k1)
    for(i in 1:n.menos.m)
    {
      for(j in 1:k1)
        M[i,j] <- X[i+(m+1),j]
    }
    
    # ##FB  recorrences
    deta.dalpha <- rep(0,n)
    deta.dtheta <- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
    
    for(i in seq1)
    {
      i.menos.ma <- i-ma
      deta.dalpha[i] <- 1 - theta%*%deta.dalpha[i.menos.ma]
      deta.dtheta[i,] <- R[(i-(m+1)),] - theta%*%deta.dtheta[i.menos.ma,]
      deta.dbeta[i,]<- M[(i-(m+1)),] - theta%*%deta.dbeta[i.menos.ma,]
    }
    
    v <- matrix(as.vector(deta.dalpha[seq1]),ncol=1)
    rR <- deta.dtheta[seq1,]
    rM <- deta.dbeta[seq1,]
    vI <- matrix(rep(1,n.menos.m),ncol=1)
    a <- muhat * (1 + prechat)
    b <- 2 + prechat
    c1 <- trigamma(a+b)
    mT <- diag(exp(eta1hat[seq1]))
    pT <- diag(exp(eta2hat[seq1]))
    D <- diag(as.vector((1+prechat)^2 * (trigamma(a) - c1))) %*% mT^2
    L <- pT %*% diag(as.vector(-(1+prechat) * (c1 + muhat*(c1-trigamma(a))))) %*% mT
    H <- diag(as.vector(muhat^2 * trigamma(a) - (1+muhat)^2 * c1 + trigamma(b))) %*% pT^2
    y.y2 <- (y[(m+1):(n-1)]/y[(m):(n-2)])
    
    Ka1a1 <- t(v) %*% D %*% v
    Ka1t <- t(v) %*% D %*% rR
    Kta1 <- t(Ka1t)
    Ka1a2 <- t(v) %*% L %*% vI
    Ka2a1 <- t(Ka1a2)
    Ka1d <- t(v) %*% L %*% y.y2
    Kda1 <- t(Ka1d)
    
    Ktt <- t(rR) %*% D %*% rR
    Kta2 <- t(rR) %*% L %*% vI
    Ka2t <- t(Kta2)
    Ktd <- t(rR) %*% L %*% y.y2
    Kdt <- t(Ktd)
    
    Ka2a2 <- t(vI) %*% H %*% vI
    Ka2d <- t(vI) %*% H %*% y.y2
    Kda2 <- t(Ka2d)
    
    Kdd <- t(y.y2) %*% H %*% y.y2
    
    Ka1b <- t(v) %*% D %*% rM
    Kba1 <- t(Ka1b)
    Ktb <- t(rR) %*% D %*% rM
    Kbt <- t(Ktb)
    Ka2b <- t(vI) %*% L %*% rM
    Kba2 <- t(Ka2b)
    Kdb <- t(y.y2) %*% L %*% rM
    Kbd <- t(Kdb)
    Kbb <- t(rM) %*% D %*% rM
    
    K <- rbind(
      cbind(Ka1a1,Ka1t,Ka1a2,Ka1d,Ka1b),
      cbind(Kta1 ,Ktt ,Kta2 ,Ktd ,Ktb ),
      cbind(Ka2a1,Ka2t,Ka2a2,Ka2d,Ka2b),
      cbind(Kda1 ,Kdt ,Kda2 ,Kdd ,Kdb ),
      cbind(Kba1 ,Kbt ,Kba2 ,Kbd ,Kbb )
    )
    
    z$K <- K
    
    # #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    prec_prev[1:n] <- z$fitted_prec
    
    X_prev<- rbind(X,X_hat)
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha1 + X_prev[n+i,]%*%as.matrix(beta) +
        (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- exp(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
      
      prec_prev[n+i] <- exp(alpha2 + delta%*%(y_prev[n+i-1]/y_prev[n+i-2]))
    }
    
    z$serie <- y
    # z$bparma <- names_par
    z$forecast <- y_prev[(n+1):(n+h1)]
    z$forecast_prev <- prec_prev[(n+1):(n+h1)]
    
  }
  
  # BPARMAX -----------------------------------------------------------------
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F)){
    
    loglik <- function(z){
      p1.mais.q1 <- p1+q1
      alpha1 <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1.mais.q1+1)]
      alpha2 <- z[p1.mais.q1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[p1.mais.q1+3]
      delta <- min(delta, -0.001)
      beta <- z[(p1.mais.q1+4):length(z)]
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + X[i,]%*%as.matrix(beta) + 
          (phi %*% (ynew[i-ar]-X[i-ar,]%*%as.matrix(beta))) + (theta %*% error[i-ma])
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- exp(eta1[seq1])
      prec <- exp(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      a.mais.b <- a + b
      
      ll <- (a-1)*log(y1) - (a.mais.b)*log1p(y1) - lgamma(a) - lgamma(b) + lgamma(a.mais.b)
      
      return(sum(ll))
    }
    
    escore <- function(z){
      p1.mais.q1 <- p1+q1
      alpha1 <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1.mais.q1+1)]
      alpha2 <- z[p1.mais.q1+2]
      alpha2 <- max(alpha2, 0.001)
      delta <- z[p1.mais.q1+3]
      delta <- min(delta, -0.001)
      beta <- z[(p1.mais.q1+4):length(z)]
      
      eta1 <- eta2 <- error <- rep(0, n)
      
      seq1 <- (m+2):n
      
      for(i in seq1){
        eta1[i] <- alpha1 + X[i,]%*%as.matrix(beta) + 
          (phi %*% (ynew[i-ar]-X[i-ar,]%*%as.matrix(beta))) + (theta %*% error[i-ma])
        error[i] <- ynew[i] - eta1[i]
        eta2[i] <- alpha2 + delta*(y[i-1]/y[i-2])
      }
      
      mu <- exp(eta1[seq1])
      prec <- exp(eta2[seq1])
      
      mu[mu < 0.01] = 0.01
      prec[prec < 0.01] = 0.01
      
      y1 <- y[seq1]
      
      a <- mu * (1 + prec)
      b <- 2 + prec
      
      n.menos.m <- n - (m+1)
      
      P <- matrix(rep(NA, n.menos.m*p1), ncol=p1)
      for(i in 1:n.menos.m)
      {
        P[i,] <- ynew[i+(m+1)-ar] - X[i+(m+1)-ar,]%*%as.matrix(beta)
      }
      
      R <- matrix(rep(NA, n.menos.m*q1), ncol=q1)
      for(i in 1:n.menos.m)
      {
        R[i,] <- error[i+(m+1)-ma]
      }
      
      k1<- length(beta)
      M <- matrix(rep(NA, n.menos.m*k1),ncol=k1)
      for(i in 1:(n-(m+1)))
      {
        for(j in 1:k1)
          M[i,j] <- X[i+(m+1),j]-sum(phi*X[i+(m+1)-ar,j])
      }
      
      ### FB recorrences
      deta.dalpha <- rep(0, n)
      deta.dphi <- matrix(0, ncol=p1, nrow=n)
      deta.dtheta <- matrix(0, ncol=q1, nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      for(i in seq1)
      {
        i.menos.ma <- i-ma
        deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i.menos.ma]
        deta.dphi[i,] <- P[(i-(m+1)),] - theta %*% deta.dphi[i.menos.ma,]
        deta.dtheta[i,] <- R[(i-(m+1)),] - theta %*% deta.dtheta[i.menos.ma,]
        deta.dbeta[i,]<- M[(i-(m+1)),] - theta%*%deta.dbeta[i.menos.ma,]
      }
      
      v <- deta.dalpha[seq1]
      rP <- deta.dphi[seq1,]
      rR <- deta.dtheta[seq1,]
      rM <- deta.dbeta[seq1,]
      vI <- rep(1, n-(m+1))
      c1 <- digamma(a+b)
      ystar <- log(y1/(1+y1))
      mustar <- digamma(a) - c1
      ydag <- mu*log(y1) - (1+mu)*log(1+y1)
      Delta <- c1 - digamma(b)
      mudag <- mu*(mustar - Delta/mu)
      f <- exp(eta1[seq1]) * ((prec+1) * (ystar-mustar))
      g <- exp(eta2[seq1]) * (ydag - mudag)
      
      Ualpha1 <- crossprod(v, f)
      Uphi <<- crossprod(rP, f)
      Utheta <<- crossprod(rR, f)
      Ualpha2 <- crossprod(vI, g)
      Udelta <<-  t(y[(m+1):(n-1)]/y[(m):(n-2)]) %*% g
      Ubeta <<- crossprod(rM, f)
      
      rval <- cbind(Ualpha1, t(Uphi), t(Utheta), Ualpha2, Udelta, t(Ubeta))
      return(rval)
    }
    
    alpha1.ini <- log(mean(y))
    phi.ini <- rep(0,p1)
    theta.ini <- rep(0,q1)
    alpha2.ini <- 3
    delta.ini <- -1
    beta.ini <- rep(0,ncol(X))
    
    reg <- c(alpha1.ini, phi.ini, theta.ini, alpha2.ini, delta.ini, beta.ini)
    
    names_par <- c("alpha1",names_phi,names_theta,"alpha2","delta",names_beta)
    
    opt <- optim(par = reg, fn = loglik, gr = escore,
                 method = "L-BFGS-B",
                 lower = c(-10,rep(-0.999,p1),rep(-0.999,q1),0.01,-Inf, rep(-10, ncol(X))),
                 upper = c(10,rep(0.999,p1),rep(0.999,q1),15,-0.01, rep(10, ncol(X))),
                 control = list(fnscale = -1, maxit = 1000))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      
      opt <- optim(par = reg, fn = loglik, #gr = escore,
                   method = "L-BFGS-B", 
                   lower = c(-10,rep(-0.999,p1),rep(-0.999,q1),0.01,-Inf, rep(-10, ncol(X))),
                   upper = c(10,rep(0.999,p1),rep(0.999,q1),15,-0.01, rep(10, ncol(X))),
                   control = list(fnscale = -1, maxit = 1000))
      
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+3+ncol(X))]
    names(coef) <- names_par
    z$coeff <- coef
    
    alpha1 <- coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    alpha2 <- coef[p1+q1+2]
    delta <- coef[p1+q1+3]
    beta <- coef[(p1+q1+4):length(coef)]
    
    z$alpha1 <- alpha1
    z$phi <- phi
    z$theta <- theta
    z$alpha2 <- alpha2
    z$delta <- delta
    z$beta <- beta
    
    errorhat <- rep(0,n)
    eta1hat <- eta2hat <- rep(NA,n)
    
    seq1 <- (m+2):n
    n.menos.m <- n - (m+1)
    
    for(i in seq1)
    {
      eta1hat[i] <- alpha1 + X[i,]%*%as.matrix(beta) + 
        (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta))) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i] - eta1hat[i]
      eta2hat[i] <- alpha2 + delta*(y[i-1]/y[i-2])
    }
    
    muhat <- exp(eta1hat[seq1])
    prechat <- exp(eta2hat[seq1])
    y1 <- y[seq1]
    
    z$fitted <- ts(c(rep(NA,(m+1)),muhat),start=start(y),frequency=frequency(y))
    z$fitted_prec <- ts(c(rep(NA,(m+1)),prechat),start=start(y),frequency=frequency(y))
    z$eta1hat <- eta1hat
    z$errorhat <- errorhat
    z$eta2hat <- eta2hat
    z$mustarhat <- digamma(muhat*(1+prechat)) - digamma(muhat*(1+prechat)+prechat+2)
    
    P <- matrix(rep(NA, n.menos.m*p1), ncol=p1)
    for(i in 1:n.menos.m)
    {
      P[i,] <- ynew[i+(m+1)-ar] - X[i+(m+1)-ar,]%*%as.matrix(beta)
    }
    
    R <- matrix(rep(NA, n.menos.m*q1), ncol=q1)
    for(i in 1:n.menos.m)
    {
      R[i,] <- errorhat[i+(m+1)-ma]
    }
    
    k1<- length(beta)
    M <- matrix(rep(NA, n.menos.m*k1),ncol=k1)
    for(i in 1:n.menos.m)
    {
      for(j in 1:k1)
        M[i,j] <- X[i+(m+1),j]-sum(phi*X[i+(m+1)-ar,j])
    }
    
    # ##FB  recorrences
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta <- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
    
    for(i in seq1)
    {
      i.menos.ma <- i-ma
      deta.dalpha[i] <- 1 - theta%*%deta.dalpha[i.menos.ma]
      deta.dphi[i,] <- P[(i-(m+1)),] - theta%*%deta.dphi[i.menos.ma,]
      deta.dtheta[i,] <- R[(i-(m+1)),] - theta%*%deta.dtheta[i.menos.ma,]
      deta.dbeta[i,]<- M[(i-(m+1)),] - theta%*%deta.dbeta[i.menos.ma,]
    }
    
    v <- matrix(as.vector(deta.dalpha[seq1]),ncol=1)
    rP <- deta.dphi[seq1,]
    rR <- deta.dtheta[seq1,]
    rM <- deta.dbeta[seq1,]
    vI <- matrix(rep(1,n.menos.m),ncol=1)
    a <- muhat * (1 + prechat)
    b <- 2 + prechat
    c1 <- trigamma(a+b)
    mT <- diag(exp(eta1hat[seq1]))
    pT <- diag(exp(eta2hat[seq1]))
    D <- diag(as.vector((1+prechat)^2 * (trigamma(a) - c1))) %*% mT^2
    L <- pT %*% diag(as.vector(-(1+prechat) * (c1 + muhat*(c1-trigamma(a))))) %*% mT
    H <- diag(as.vector(muhat^2 * trigamma(a) - (1+muhat)^2 * c1 + trigamma(b))) %*% pT^2
    y.y2 <- (y[(m+1):(n-1)]/y[(m):(n-2)])
    
    Ka1a1 <- t(v) %*% D %*% v
    Ka1p <- t(v) %*% D %*% rP
    Kpa1 <- t(Ka1p)
    Ka1t <- t(v) %*% D %*% rR
    Kta1 <- t(Ka1t)
    Ka1a2 <- t(v) %*% L %*% vI
    Ka2a1 <- t(Ka1a2)
    Ka1d <- t(v) %*% L %*% y.y2
    Kda1 <- t(Ka1d)
    
    Kpp <- t(rP) %*% D %*% rP
    Kpt <- t(rP) %*% D %*% rR
    Ktp <- t(Kpt)
    Kpa2 <- t(rP) %*% L %*% vI
    Ka2p <- t(Kpa2)
    Kpd <- t(rP) %*% L %*% y.y2
    Kdp <- t(Kpd)
    
    Ktt <- t(rR) %*% D %*% rR
    Kta2 <- t(rR) %*% L %*% vI
    Ka2t <- t(Kta2)
    Ktd <- t(rR) %*% L %*% y.y2
    Kdt <- t(Ktd)
    
    Ka2a2 <- t(vI) %*% H %*% vI
    Ka2d <- t(vI) %*% H %*% y.y2
    Kda2 <- t(Ka2d)
    
    Kdd <- t(y.y2) %*% H %*% y.y2
    
    Ka1b <- t(v) %*% D %*% rM
    Kba1 <- t(Ka1b)
    Kpb <- t(rP) %*% D %*% rM
    Kbp <- t(Kpb)
    Ktb <- t(rR) %*% D %*% rM
    Kbt <- t(Ktb)
    Ka2b <- t(vI) %*% L %*% rM
    Kba2 <- t(Ka2b)
    Kdb <- t(y.y2) %*% L %*% rM
    Kbd <- t(Kdb)
    Kbb <- t(rM) %*% D %*% rM
    
    K <- rbind(
      cbind(Ka1a1,Ka1p,Ka1t,Ka1a2,Ka1d,Ka1b),
      cbind(Kpa1 ,Kpp ,Kpt ,Kpa2 ,Kpd ,Kpb ),
      cbind(Kta1 ,Ktp ,Ktt ,Kta2 ,Ktd ,Ktb ),
      cbind(Ka2a1,Ka2p,Ka2t,Ka2a2,Ka2d,Ka2b),
      cbind(Kda1 ,Kdp ,Kdt ,Kda2 ,Kdd ,Kdb ),
      cbind(Kba1 ,Kbp ,Kbt ,Kba2 ,Kbd ,Kbb )
    )
    
    z$K <- K
    
    # #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    prec_prev[1:n] <- z$fitted_prec
    
    X_prev<- rbind(X,X_hat)
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha1 + X_prev[n+i,]%*%as.matrix(beta) +
        (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta))) +
        (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- exp(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
      
      prec_prev[n+i] <- exp(alpha2 + delta%*%(y_prev[n+i-1]/y_prev[n+i-2]))
    }
    
    z$serie <- y
    # z$bparma <- names_par
    z$forecast <- y_prev[(n+1):(n+h1)]
    z$forecast_prev <- prec_prev[(n+1):(n+h1)]
    
  }
  
  ################# RESIDUALS ##################
  
  res1 <- y-z$fitted
  vary <- z$fitted * (1+z$fitted)/z$fitted_prec
  
  z$resid1 <- (res1/sqrt(vary))[(m+1):n]
  
  a1 <- y * (1 + z$fitted_prec)
  b1 <- 2 + z$fitted_prec
  
  a2 <- z$fitted * (1 + z$fitted_prec)
  b2 <- 2 + z$fitted_prec
  
  # corrigir isso aqui
  l_tilde <- suppressWarnings(log(dbetapr(x = y,shape1 = a1,shape2 = b1,scale = 1.0)))
  l_hat <- suppressWarnings(log(dbetapr(x = y,shape1 = a2,shape2 = b2,scale = 1.0)))
  
  dt <- (l_tilde-l_hat)[(m+1):n]
  dt[which(dt<0)] <- 0
  
  z$resid2 <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))
  
  z$resid3 <- as.vector(qnorm(pbetapr(q = y[(m+1):n],shape1 = a2[(m+1):n],
                                      shape2 = b2[(m+1):n],scale = 1.0)))
  
  
  ################# SAIDAS ##################
  
  z$deviance <- 2*sum(dt)
  
  k <- sum(c(p1,q1),na.rm=T) # numero de parametros
  
  z$rank <- k
  
  Kchol <- tryCatch(chol(z$K), error = function(e) return("error"),
                    warning = function(o) return("error"))
  
  if(Kchol[1] == "error")
  {
    z$vcov <- try(chol2inv(z$K))
    warning("We have problems with information matrix inversion!")
    
  }else{
    vcov <- try(chol2inv(Kchol))
    z$vcov <- vcov
  }
  
  stderror <- sqrt(diag(z$vcov))
  z$stderror <- stderror
  
  z$zstat <- abs(z$coef/stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat))
  
  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])
  
  if(any(is.na(X)==F))
  {
    z$aic <- -2*z$loglik+2*(p1+q1+3+length(beta))
    z$bic <- -2*z$loglik+log(n)*(p1+q1+3+length(beta))
  }else{
    z$aic <- -2*z$loglik+2*(p1+q1+3)
    z$bic <- -2*z$loglik+log(n)*(p1+q1+3)
  }
  
  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),
                              round(z$pvalues,4))
  colnames(model_presentation) <- c("Estimate","Std. Error","z value","Pr(>|z|)")
  
  z$model <- model_presentation
  
  if(diag == 0){
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)),quote=F)
  }
  
  return(z)
  
}