# typhoid-like vaccine efficacy model

# helper functions
logit = function(p,max=1){log(p/max/(1-p/max))}
invlogit=function(l,max=1){max/(1+exp(-l))}
# the logits need the max to account for the possibility that the immune correlate isn't very good
# or there is no sterilizing immunity even at very high correlate.
# as in, that even at very high correlate, there may be some people who are not immune, because either
# no one can be immune, or we're measuring not quite the right thing
# 
#   if you knew the max efficacy came from people failing on "vaccine take" 
#   (like with live polio vax, where some people don't get infected), then you
#.  should model VE as a mixture, with some people having no response 
#.  and some having a correlate response that you model this way.
#.  
#.  from the TCV phase II immunogenicity studies I skimmed, it looks like both could be going on.
#.    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8298254/ if you subtract control seroconversion
#.    from vax seroconversion (table 3), you get 80%. So this would look like the max VE comes from vaccine failure
#.    But this one, also from burkina https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7762715/, is more like
#.    90% take (table 3 again). Intervals allow 80% but it looks a bit higher.
#.    and nepal https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8551681/ is in the mid-90% take
#
#  In what follows, for simplicity, I'm gonna assume 100% take and that all the action is in the logits.
#. we could revisit this.


# vaccine efficacy vs a correlate (either known, or we can assume one exists and treat it as a hidden/linking variable)
# logit-log model. General across pathogens
VE_vs_correlate = function(correlate=1,
                           logit_min_ve=0,
                           delta_logit_ve=1,
                           ve_max=1,
                           trans='logit'){
  
  if(is.data.frame(correlate)){ # to make nice data frames under composition
    times=correlate$time
    correlate=correlate$correlate
  }
  
  logit_ve = logit_min_ve + delta_logit_ve*log10(correlate) # logit-linear model
  res=logit_ve
  
  if(trans=='logit'){
  }else if(trans=='linear'){
    ve = invlogit(logit_ve,max=ve_max)
    res=ve
  }
  
  if(exists('times')){ # to make nice data frames under composition
    return(data.frame(time=times,VE=res))
  }else{
    return(res)
  }
}
#   (log-log) is more commonly used in analyzing trials because people (rightly)
#   don't want to assume that VE must be greater than zero.
#   but once you know that you have a good measure of an effective vaccine,
#   I don't think it makes sense to bother with the possibility of negative VE
#   Anyway, this doesn't really matter. So change to logs if you prefer. Logit
#.  is just a better transform at VE close to 1.

# default 50% VE
VE_vs_correlate(trans='logit')
VE_vs_correlate(trans='linear')


# power law waning
# if you make alpha much greater than 1, this returns to an exponential with time constant mean_decay.
# for small alpha, this has a slower decay at long times
#
#. NOTE THAT THE SHAPE AT EARLY TIMES I DREW ON THE NAPKIN IS WRONG :-(
#. I sadly mis-remembered how the functional form works for the simple power law.
#. I could justify something fancier, but I don't think we need to. 
correlate_vs_time = function(time=seq(1,10*365,by=1)/365, # 10 years
                             correlate_0 = 10^3.4, 
                             mean_decay=2, # 2 years
                             alpha=1,
                             old_model=FALSE){ # shape parameter. smaller is less exponential
  
  if(!old_model){
    correlate_vs_time = kronecker(correlate_0 , (1+time/(mean_decay*alpha))^(-alpha))
  } else {
    correlate_vs_time = kronecker(correlate_0 , pmin(1,(time/mean_decay)^(-alpha)))
  }
  return(data.frame(time=time,correlate=correlate_vs_time))
}

plot(correlate_vs_time(alpha=100)) # exponential decay with the same short-time behavior
lines(correlate_vs_time(alpha=1)) # default power law decay

plot(log10(correlate_vs_time(alpha=100,time=seq(1*365,80*365,by=1)/365)$time),
     log10(correlate_vs_time(alpha=100,time=seq(1*365,80*365,by=1)/365)$correlate)) # exponential decay with the same short-time behavior
lines(log10(correlate_vs_time(alpha=1,time=seq(1*365,80*365,by=1)/365)$time),
      log10(correlate_vs_time(alpha=1,time=seq(1*365,80*365,by=1)/365)$correlate)) # default power law decay
# ganna have to refit polio data...

plot(correlate_vs_time(alpha=1),ylim=c(1,10^3.5)) # exponential decay with the same short-time behavior
lines(correlate_vs_time(alpha=0.8,mean_decay=1,old_model = TRUE)) # default power law decay
lines(correlate_vs_time(alpha=0.87,mean_decay=30/365,old_model = TRUE)) # default power law decay


## put them together.
## I'm just gonna make up a correlate that starts at 10^3.4, because 
## that's the median 28 day anti-VI IgG titer here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7762715/
## the rest of the parameters are made up to make VE look like what you told me, very roughly
VE_vs_time = VE_vs_correlate(correlate=correlate_vs_time(correlate_0 = 10^3.4,time=seq(1,10*365,by=1)/365,
                                                         alpha=1,mean_decay=2),
                             logit_min_ve=-29,delta_logit_ve = 10,ve_max=0.8,
                             trans='linear')
plot(VE_vs_time,ylim=c(0,1),xlab='years')
## this is 5 free parameters, vs the 3 you are currently using (amplitude, width of box, decay rate), right?
# the extra two are the alpha shape parameter that describes deviation from exponential / power law of tail,
# and the relationship between the correlate and the VE

# the reward for the two extra parameters, if you can fit them, is you get a better long-time waning model, 
# and you can generalize VE to populations with different correlate distributions than the clinical trial pop

## if you want to do one extra parameter, keep the box VE but try to fit a power law tail

## also, if we don't want to assume we know anything about the correlate, we could just assume
## there is a magic number that goes from 1 to 10 (0 to 1 in log space) and fit VE vs time directly



# sequential boosting and what does fast vs slow waning mean?
VE_vs_time_boost = VE_vs_correlate(correlate=correlate_vs_time(correlate_0 = 10^3.8,time=seq(1,10*365,by=1)/365,
                                                         alpha=1,mean_decay=2),
                             logit_min_ve=-29,delta_logit_ve = 10,ve_max=0.8,
                             trans='linear')
plot(VE_vs_time,ylim=c(0,1),xlab='years')
lines(VE_vs_time_boost)



corr_v_time = correlate_vs_time(correlate_0 = 10^3.4,time=seq(1,50*365,by=1)/365,
                                alpha=1,mean_decay=2)
corr_v_time_boost = correlate_vs_time(correlate_0 = 10^3.8,time=seq(1,50*365,by=1)/365,
                                                               alpha=1,mean_decay=2)
corr_v_time_boost_two = correlate_vs_time(correlate_0 = 10^4,time=seq(1,50*365,by=1)/365,
                                      alpha=1,mean_decay=2)

plot(corr_v_time,ylim=c(1,10^4),xlab='years')
lines(corr_v_time_boost)
lines(corr_v_time_boost_two)
abline(h=1000)


sequential_boosting = corr_v_time
sequential_boosting[(5*365+1):length(sequential_boosting$time),2] = 
  corr_v_time_boost[1:(45*365),2]

sequential_boosting[(18*365+1):length(sequential_boosting$time),2] = 
  corr_v_time_boost_two[1:(32*365),2]

plot(sequential_boosting$time,log10(sequential_boosting$correlate),xlab='years',ylab='log10(anti-VI titer)')
abline(h=3)


VE_vs_time_sequential = VE_vs_correlate(correlate=sequential_boosting,
                                   logit_min_ve=-29,delta_logit_ve = 10,ve_max=0.8,
                                   trans='linear')
plot(VE_vs_time_sequential,ylim=c(0,1),xlab='years')
abline(h=0.5)




