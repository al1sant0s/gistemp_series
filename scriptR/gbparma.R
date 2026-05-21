setwd("C:/Users/Kleber Santos/Dropbox/PesquisaAlissonKleber/scriptR")

gbparma <- function(y, ar=NA, ma=NA, mu.link="log", prec.link="log", h1=6, X=NA, 
                    X_hat=NA,resid=3,diag=1){
  
  setwd("C:/Users/Kleber Santos/Dropbox/PesquisaAlissonKleber/scriptR")
  source("gbparma.fit.corrigido.r")
  
  if (min(y) <= 0)
    stop("OUT OF RANGE (0,Inf)!")
  
  if(is.ts(y)==T)
  {
    freq<-frequency(y)
  }else stop("data can be a time-series object")
  
  if(any(is.na(ar))==F) names_phi <- c(paste("phi",ar,sep=""))
  
  if(any(is.na(ma))==F) names_theta <- c(paste("theta",ma,sep=""))
  
  if(any(is.na(X))==F)
  {
    names_beta <- c(paste("beta",1:ncol( as.matrix(X) ),sep=""))
  }
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  linktemp1 <- substitute(mu.link)
  if (!is.character(linktemp1))
  {
    linktemp1 <- deparse(linktemp1)
    if (linktemp1 == "mu.link")
      linktemp1 <- eval(mu.link)
  }
  if (any(linktemp1 == c("log", "sqrt", "identity"))){
    stats <- make.link(linktemp1)
  } else stop(paste(linktemp1, "prec.link not available, available links are log, sqrt and identity"))
  
  mu.link <- structure(list(link = linktemp1,
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv,
                         mu.eta = stats$mu.eta,
                         diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  ))
  
  linkfun1 <- mu.link[[2]]
  linkinv1 <- mu.link[[3]]
  mu.eta1 <-  mu.link[[4]]
  diflink1 <- mu.link[[5]]
  
  linktemp2 <- substitute(prec.link)
  if (!is.character(linktemp2))
  {
    linktemp2 <- deparse(linktemp2)
    if (linktemp2 == "prec.link")
      linktemp2 <- eval(prec.link)
  }
  if (any(linktemp2 == c("log", "sqrt", "identity"))){
    stats <- make.link(linktemp2)
  } else stop(paste(linktemp2, "prec.link not available, available links are log, sqrt and identity"))
  
  prec.link <- structure(list(link = linktemp2,
                          linkfun = stats$linkfun,
                          linkinv = stats$linkinv,
                          mu.eta = stats$mu.eta,
                          diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  ))
  
  linkfun2 <- prec.link[[2]]
  linkinv2 <- prec.link[[3]]
  mu.eta2 <-  prec.link[[4]]
  diflink2 <- prec.link[[5]]
  
  fit <- gbparma.fit(y, ar, ma, mu.link, prec.link, names_phi, names_theta, 
                     names_beta, h1, X, X_hat, resid=resid, diag=diag)
  
  return(fit)
}

