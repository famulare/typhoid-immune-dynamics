library(tidyverse)

datapath <-"your path here"


datafile <- "TypoidCaseData_github_09.30.21.csv"

#load data
raw <- read.csv(paste(datapath,datafile,sep=""),head=TRUE) 

#filter for country/age-specific estimates
#  filter(age >=10 & age<16) %>%
#  filter(Country == "Nepal") 

 # filter(bldculres == "typhi") %>%
 # filter(bldculres == "paratyphi") %>%
  
 # filter(Hospitalized =="Yes") %>% 
#  filter(Hospitalized =="No") %>% 
  
  #filter young ages 
  #filter(age<5)
# filter(age >=5 & age<16) 
 #filter(age >=16)


id <- raw$index_id;
age <- raw$age;
# nsmpl <- raw$nVisits; # not same as nr. observed y
hlye_iga <- cbind(raw$HlyE_IgA_visit1,raw$HlyE_IgA_visit2,
                  raw$HlyE_IgA_visit3,raw$HlyE_IgA_visit4,
                  raw$HlyE_IgA_visit5,raw$HlyE_IgA_visit6,
                  raw$HlyE_IgA_visit7);
hlye_igg <- cbind(raw$HlyE_IgG_visit1,raw$HlyE_IgG_visit2,
                  raw$HlyE_IgG_visit3,raw$HlyE_IgG_visit4,
                  raw$HlyE_IgG_visit5,raw$HlyE_IgG_visit6,
                  raw$HlyE_IgG_visit7);
lps_iga  <- cbind(raw$LPS_IgA_visit1, raw$LPS_IgA_visit2,
                  raw$LPS_IgA_visit3, raw$LPS_IgA_visit4,
                  raw$LPS_IgA_visit5, raw$LPS_IgA_visit6,
                  raw$LPS_IgA_visit7);
lps_igg  <- cbind(raw$LPS_IgG_visit1, raw$LPS_IgG_visit2,
                  raw$LPS_IgG_visit3, raw$LPS_IgG_visit4,
                  raw$LPS_IgG_visit5, raw$LPS_IgG_visit6,
                  raw$LPS_IgG_visit7);
mp_iga   <- cbind(raw$MP_IgA_visit1,  raw$MP_IgA_visit2,
                  raw$MP_IgA_visit3,  raw$MP_IgA_visit4,
                  raw$MP_IgA_visit5,  raw$MP_IgA_visit6,
                  raw$MP_IgA_visit7);
mp_igg   <- cbind(raw$MP_IgG_visit1,  raw$MP_IgG_visit2,
                  raw$MP_IgG_visit3,  raw$MP_IgG_visit4,
                  raw$MP_IgG_visit5,  raw$MP_IgG_visit6,
                  raw$MP_IgG_visit7);
vi_igg   <- cbind(raw$Vi_IgG_visit1,  raw$Vi_IgG_visit2,
                  raw$Vi_IgG_visit3,  raw$Vi_IgG_visit4,
                  raw$Vi_IgG_visit5,  raw$Vi_IgG_visit6,
                  raw$Vi_IgG_visit7);
visit.t <- cbind(raw$TimeInDays_visit1,raw$TimeInDays_visit2,
                raw$TimeInDays_visit3,raw$TimeInDays_visit4,
                raw$TimeInDays_visit5,raw$TimeInDays_visit6,
                raw$TimeInDays_visit7);

test <- data.frame(visit.t) %>% 
  rowwise() %>%
  mutate(sum.na = sum(is.na(X1), is.na(X2), is.na(X3), is.na(X4), is.na(X5), is.na(X6)))

pr.nm <- c("y0","y1","t1","alpha","shape");
ab.nm <- rbind(c("HlyE","IgA"),c("HlyE","IgG"),
               c("LPS","IgA"),c("LPS","IgG"),
               c("MP","IgA"),c("MP","IgG"),c("Vi","IgG"));
nsubj <- length(id);
ntest <- 7;
maxsmpl <- 7; # maximum nr. samples per subject
nsmpl <- rep(NA,nsubj+1);
smpl.t <- array(NA,dim=c(nsubj+1,maxsmpl));
smpl.y <- array(NA,dim=c(nsubj+1,maxsmpl,ntest));
indx   <- array(NA,dim=c(nsubj+1,maxsmpl));


for(k.subj in 1:nsubj){
  exst <- sort(which(!is.na(visit.t[k.subj,])));
  nsmpl[k.subj] <- length(exst);
  indx[k.subj,1:nsmpl[k.subj]] <- exst;
  smpl.t[k.subj,1:nsmpl[k.subj]] <- visit.t[k.subj,exst];
  smpl.y[k.subj,1:nsmpl[k.subj],1] <- hlye_iga[k.subj,exst];
  smpl.y[k.subj,1:nsmpl[k.subj],2] <- hlye_igg[k.subj,exst];
  smpl.y[k.subj,1:nsmpl[k.subj],3] <- lps_iga[k.subj,exst];
  smpl.y[k.subj,1:nsmpl[k.subj],4] <- lps_igg[k.subj,exst];
  smpl.y[k.subj,1:nsmpl[k.subj],5] <- mp_iga[k.subj,exst];
  smpl.y[k.subj,1:nsmpl[k.subj],6] <- mp_igg[k.subj,exst];
  smpl.y[k.subj,1:nsmpl[k.subj],7] <- vi_igg[k.subj,exst];
}


nsmpl[nsubj+1] <- 3;
smpl.t[nsubj+1,] <- c(5,30,90,NA,NA,NA,NA);
age <- c(age,10); # just made this up
smpl.y[smpl.y==0] <- 0.01;  # remove y=0
logy <- log(smpl.y);

npar <- 5; # y0, y1, t1, alpha, shape
ndim <- npar; # size of cov mat
mu.hyp   <- array(NA,dim=c(ntest,ndim));
prec.hyp <- array(NA,dim=c(ntest,ndim,ndim));
omega    <- array(NA,dim=c(ntest,ndim,ndim));
wishdf   <- rep(NA,ntest);
prec.logy.hyp <- array(NA,dim=c(ntest,2));
#                          log(c(y0,  y1,    t1,  alpha, shape-1))
for(k.test in 1:ntest){
  mu.hyp[k.test,] <-        c( 1.0,   7.0,   1.0,  -4.0, -1.0);
  prec.hyp[k.test,,] <- diag(c(1.0,   0.00001,   1.0,   0.001,  1.0));
  omega[k.test,,] <-    diag(c(1.0,  50.0,   1.0,   10.0,  1.0));
  wishdf[k.test] <- 20;
  prec.logy.hyp[k.test,] <- c(4.0,1.0);
}

longdata <- list("smpl.t"=smpl.t, "logy"=logy,
                 "ntest"=ntest,"nsmpl"=nsmpl, "nsubj"=nsubj+1,"ndim"=ndim,
                 "mu.hyp"=mu.hyp, "prec.hyp"=prec.hyp,
                 "omega"=omega, "wishdf"=wishdf,
                 "prec.logy.hyp"=prec.logy.hyp);
