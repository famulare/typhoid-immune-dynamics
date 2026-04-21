library(Hmisc)

ticklab <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,1e2,2e2,5e2,1e3,
             2e3,5e3,1e4,2e4,5e4,1e5,2e5,5e5,1e6,2e6,5e6,1e7,2e7,5e7,1e8,
             2e8,5e8,1e9,2e9,5e9,1e10);
tickpos <- log10(ticklab);

bt <- function(y0,y1,t1) log(y1/y0)/t1

ab <- function(t,y0,y1,t1,alpha,shape) {
  beta <- bt(y0,y1,t1);
  yt <- 0;
  if(t <= t1) yt <- y0*exp(beta*t);
  if(t > t1) yt <- (y1^(1-shape)-(1-shape)*alpha*(t-t1))^(1/(1-shape));
  return(yt);
}

sero <- function(n,tvec,y0,y1,t1,alpha,shape) {
  tmp <- rep(NA,length(tvec));
  for(k in 1:length(tvec)) {
    tmp[k] <- ab(tvec[k],y0[n],y1[n],t1[n],alpha[n],shape[n])
  }
  return(tmp);
}

qsero <- function(t,q,y0,y1,t1,alpha,shape) {
  nmc <- length(y0);
  tmp <- rep(NA,nmc);
  for(k in 1:nmc)
    tmp[k] <- ab(t,y0[k],y1[k],t1[k],alpha[k],shape[k]);
  if(length(q)==1) if(q=="mean") return(exp(mean(log(tmp))));
  return(quantile(tmp,q));
}

serocourse <- function(tvec,q,y0,y1,t1,alpha,shape) {
  n.pts <- length(tvec);
  tmp <- rep(NA,n.pts);
  for(k in 1:length(tvec)) {
    tmp[k] <- qsero(tvec[k],q,y0,y1,t1,alpha,shape);
  }
  return(tmp);
}

wdens <- function(w,y1,alpha,shape){
  rho <- 1/(shape-1);
  dens <- 1/(w*gamma(rho))*(w*rho/alpha)^rho*exp(-w*rho*y1^(-1/rho)/alpha);
  return(dens);
}

wdistquan <- function(wvec,qvec,y1,alpha,shape){
  densvec <- array(NA,dim=c(length(wvec),length(qvec)));
  for(k.w in 1:length(wvec)){
    densvec[k.w,] <- quantile(wdens(wvec[k.w],y1,alpha,shape),qvec,na.rm=TRUE);
  }
  return(densvec);
}

wlogdistquan <- function(logwvec,qvec,y1,alpha,shape){
  densvec <- array(NA,dim=c(length(logwvec),length(qvec)));
  for(k.w in 1:length(logwvec)){
    densvec[k.w,] <- 10^logwvec[k.w]*log(10)*
        quantile(wdens(10^logwvec[k.w],y1,alpha,shape),qvec,na.rm=TRUE);
  }
  return(densvec);
}
