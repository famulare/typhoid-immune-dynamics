model {
 for(subj in  1:nsubj){
  for(test in 1:ntest) {
   beta[subj,test] <- log(y1[subj,test]/y0[subj,test])/t1[subj,test]
   for(obs in 1:nsmpl[subj]){
     mu.logy[subj,obs,test] <- ifelse(step(t1[subj,test]-smpl.t[subj,obs]),
       log(y0[subj,test])+(beta[subj,test]*smpl.t[subj,obs]),
       1/(1-shape[subj,test])*log(y1[subj,test]^(1-shape[subj,test])-
        (1-shape[subj,test])*alpha[subj,test]*(smpl.t[subj,obs]-t1[subj,test])))
     logy[subj,obs,test] ~ dnorm(mu.logy[subj,obs,test],prec.logy[test])
   }
   y0[subj,test]    <- exp(par[subj,test,1])
   y1[subj,test]    <- y0[subj,test]+exp(par[subj,test,2])
   t1[subj,test]    <- exp(par[subj,test,3])
   alpha[subj,test] <- exp(par[subj,test,4])
   shape[subj,test] <- exp(par[subj,test,5])+1
   par[subj,test,1:ndim] ~ dmnorm(mu.par[test,],prec.par[test,,])
  }
 }
 for(test in 1:ntest) {
  mu.par[test,1:ndim] ~ dmnorm(mu.hyp[test,],prec.hyp[test,,])
  prec.par[test,1:ndim,1:ndim] ~ dwish(omega[test,,],wishdf[test])
  prec.logy[test] ~ dgamma(prec.logy.hyp[test,1],prec.logy.hyp[test,2])
 }
}
