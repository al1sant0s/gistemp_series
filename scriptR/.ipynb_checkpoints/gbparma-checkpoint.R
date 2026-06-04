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
  
  fit <- gbparma.fit(y, ar, ma, mu.link, prec.link, names_phi, names_theta, 
                     names_beta, h1, X, X_hat, resid=resid, diag=diag)
  
  return(fit)
}

