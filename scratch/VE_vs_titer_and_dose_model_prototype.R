#' ---
#' title: "VE vs IgG and bacilli dose prototype"
#' output:
#'   md_document
#' knit: (function(input, ...) {
#'   out <- rmarkdown::render(input, output_dir = "./docs/blog/posts", ...)
#'   return(out)})
#' ---
#'
#' 
#' # VE vs IgG and bacilli dose prototype
#' 
#' This is a full run through of a vaccine efficacy vs anti-Vi IgG concentration and salmonella typhi dose model. 
#' 
#' **Why?** Show all the pieces and figure out where the challenges sit.
#' 
#' <!-- more -->
#' 
#' 

#+ echo=TRUE, message=FALSE, results = 'hide'
# set up the environment
library(tidyverse)
library(scales)

#' ## Define the model

#+ echo=TRUE, message=FALSE, results = 'hide'

# dose response functions for various outcomes

  ## set up a parameters data.frame that we will build up as we write the model

  # probability of fever given dose (quailes strain, skim milk) 
  # default used in Typhoidsim, from https://qmrawiki.canr.msu.edu/experiments/salmonella-typhi recommended model
  # param_df=data.frame(value=c(n50_fever_given_dose=1.11e6, 
  #                             alpha_fever_given_dose=0.175,  
  #                             gamma_fever_given_dose=0.2), # fuzzy guestimate from  https://pmc.ncbi.nlm.nih.gov/articles/PMC5720597/
  #                     component = 'p_fever_given_dose')  |>
  #   

  # I think there is some unknown background history of typhoid fever in the dose response challenge study cohorts
  # trying to make a crappy guess at what an equivalent mean titer might look like if inmates saw typhoid in childhood # https://www.cdc.gov/mmwr/preview/mmwrhtml/mm4840a1.htm

  # assume titer =50 because that's half of waned level after 5 years in Bangladesh teens...
  
  # rescale variable to remove immune confound
  gamma=0.4
  # rescale alpha
  alpha_prime=50^gamma*0.175
  
  # rescale n50
  N50_prime=1.11e6*(2^(1/alpha_prime)-1)/(2^(1/0.175)-1)
  N50_prime
  
  # note that this rescaling actually gets the control of the Jin2017 study decently right if bicarbonate right before challenge is worth ~a factor of 10 in dose
  param_df=data.frame(value=c(n50_fever_given_dose=N50_prime, 
                              alpha_fever_given_dose=alpha_prime,  
                              gamma_fever_given_dose=gamma), # guesstimated
                      component = 'p_fever_given_dose')  |>
  
  # probability of infection given dose
  # not well measured, but relevant for transmission. 
  # guesstimated from table 2 of https://doi.org/10.1056/NEJM197009242831306  
  # and stool vax vs control https://pmc.ncbi.nlm.nih.gov/articles/PMC5720597/
    rbind(data.frame(value=c(n50_infection_given_dose=N50_prime/10, 
                             alpha_infection_given_dose=alpha_prime*2, 
                             gamma_infection_given_dose=gamma*0.2), 
                     component = 'p_infection_given_dose') )
    

p_outcome_given_dose = function(dose=1e4,  CoP=1, outcome = 'fever_given_dose', params=param_df){
  
  if(outcome == 'fever_given_dose'){
    
    p = 1 - (1+dose*(2^(1/params['alpha_fever_given_dose','value'])-1)/params['n50_fever_given_dose','value'])^(-params['alpha_fever_given_dose','value']/(CoP^params['gamma_fever_given_dose','value']))
    
  } else if (outcome == 'infection_given_dose'){
    
    p = 1 - (1+dose*(2^(1/params['alpha_infection_given_dose','value'])-1)/params['n50_infection_given_dose','value'])^(-params['alpha_infection_given_dose','value']/(CoP^params['gamma_infection_given_dose','value']))
    
  } else if (outcome == 'fever_given_infection'){
    
    p = (1 - (1+dose*(2^(1/params['alpha_fever_given_dose','value'])-1)/params['n50_fever_given_dose','value'])^(-params['alpha_fever_given_dose','value']/(CoP^params['gamma_fever_given_dose','value']))) /
      (1 - (1+dose*(2^(1/params['alpha_infection_given_dose','value'])-1)/params['n50_infection_given_dose','value'])^(-params['alpha_infection_given_dose','value']/(CoP^params['gamma_infection_given_dose','value'])))
    
  }

  return(p)
}

# please forgive this disgusting data structure! It just kinda... happened.
plot_dat = data.frame(dose = 10^seq(0,10,by=0.05),
                      CoP=1) |>
  mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_dose'),
         outcome = 'fever given dose') |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=1) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=1) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=10) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=10) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=10) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=50) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=50) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=50) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |> 
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=200) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=200) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=200) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |> 
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=2000) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=2000) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=2000) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |> 
  mutate(CoP=factor(CoP)) |>
  mutate(outcome = factor(outcome, levels =c('infection given dose','fever given dose','fever given infection')))

ggplot(plot_dat ) +
  geom_line(aes(x=dose,y=p,group=CoP,color=CoP)) +
  facet_wrap('outcome') +
  theme_bw() +
  scale_x_continuous(trans='log10', breaks=10^seq(1,10,by=2),minor_breaks = NULL,labels = scales::trans_format("log10", math_format(10^.x)) ) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
  ylab('probability of outcome') +
  xlab('dose [bacilli]') +
  labs(color='Anti-Vi IgG\n[EU/ml]' )
ggsave('scratch/figures/dose_response.png',units='in',width=7, height=2.5)

# incubation period vs dose
# this is limited data on this, but it's a thing https://doi.org/10.1056/NEJM197009242831306
incubation_period_distribution_vs_dose = function(median,log_sd,params=param_df){
  # TO DO
}

#' ## Vaccine efficacy
#' 
#+ echo=TRUE, message=FALSE, results = 'hide'
vaccine_efficacy = function(dose=1e4,  CoP=1, outcome = 'fever_given_dose', params=param_df, CoP_control=1){
  
  VE = 1 - p_outcome_given_dose(dose,CoP=CoP,outcome = outcome)/p_outcome_given_dose(dose,CoP=CoP_control,outcome = outcome)
  return(VE)
}

# please forgive this disgusting data structure! It just kinda... happened.
plot_dat = data.frame(dose = 10^seq(0,10,by=0.05),
                      CoP=50) |>
  mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
         outcome = 'fever given dose') |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=200) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=2000) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=10) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
             CoP=50) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=200) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=2000) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=10) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=50) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=200) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=2000) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=10) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |> 
  mutate(CoP=factor(CoP)) |>
  mutate(outcome = factor(outcome, levels =c('infection given dose','fever given dose','fever given infection')))


ggplot(plot_dat) +
  geom_line(aes(x=dose,y=ve,group=CoP,color=CoP)) +
  facet_grid('~outcome') +
  theme_bw() +
  scale_x_continuous(trans='log10', breaks=10^seq(1,10,by=2),minor_breaks = NULL,labels = scales::trans_format("log10", math_format(10^.x)) ) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
  ylab('VE') +
  xlab('dose [bacilli]') +
  labs(color='Anti-Vi IgG\n[EU/ml]' )
ggsave('scratch/figures/vaccine_efficacy_vs_dose_CoP.png',units='in',width=7, height=2.5)

#' This kinda looks like it works!

#' # Waning model!
#' 

# time series function for waning
# units in days
t_peak=21
t = seq(0,1e3,by=1)

param_df = param_df |>
  rbind(data.frame(value=c(T_decay=9.6,k=1.04,
                           beta_C_max_age_mean=-0.055, beta_T_decay_age_mean=0.056, beta_k_age_mean=-0.037,
                           C_min=1,T_rise=1.5,t_start=3,t_peak=28), 
                   component = 'anti_Vi_IgG_elisa') )

titer_vs_time = function(t,C_max, age_years,params=param_df,C_pre=1){
  
  Cma=C_max*exp(params['beta_C_max_age_mean','value'] * age_years)
  Tda = params['T_decay','value'] *exp(params['beta_T_decay_age_mean','value']  * age_years)
  ka = params['k','value']*exp(params['beta_k_age_mean','value']  * age_years)

  # power law decay
  titer = params['C_min','value']  + (Cma-params['C_min','value'])*(1+(t-params['t_peak','value'] )/(ka*Tda))^(-ka)
  
  # simple logistic interpolator rise (this is just for continuity/realism. plays no role in the model)
  titer[t<params['t_peak','value'] ] =
    C_pre + (Cma-C_pre)*
    (1/(1+exp(-(t[t<params['t_peak','value'] ]-params['t_start','value'] *5)/params['T_rise','value'])) - 1/(1+exp(-(0-params['t_start','value']*5)/params['T_rise','value'])))/
    (1/(1+exp(-(params['t_peak','value'] -params['t_start','value']*5)/params['T_rise','value'])) - 1/(1+exp(-(0-params['t_start','value']*5)/params['T_rise','value'])))
  
  return(titer)
}

plot_dat = data.frame(t=seq(28,20*365,by=1),
                      dose = 1e3,
                      age_years=8,
                      dotname='Blantyre.Malawi',
                      vaccine='Typbar-TCV') |>
  mutate(CoP=titer_vs_time(t=t,C_max=4307,age_years=age_years)) |>
  mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
         outcome = 'fever given dose') |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 1e3,
                   age_years=1.4,
                   dotname='Blantyre.Malawi',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=4307,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 5e4,
                   age_years=8,
                   dotname='Mirpur.Dhaka.Bangladesh',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=4307,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 5e4,
                   age_years=1.4,
                   dotname='Mirpur.Dhaka.Bangladesh',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=4307,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                     dose = 1e3,
                     age_years=8,
                     dotname='Blantyre.Malawi',
                     vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=4307,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 1e3,
                   age_years=1.4,
                   dotname='Blantyre.Malawi',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=4307,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 5e4,
                   age_years=8,
                   dotname='Mirpur.Dhaka.Bangladesh',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=4307,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 5e4,
                   age_years=1.4,
                   dotname='Mirpur.Dhaka.Bangladesh',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=4307,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  mutate(year=t/365) |>
  mutate(year_groups = if_else(year<=2,'0-2',if_else(year<=5,'3-5','>5'))) |>
  group_by(year_groups,outcome,dotname,age_years,vaccine) |>
  mutate(mean=if_else(year<=5,mean(ve),NA)) |>
  ungroup()

ggplot(plot_dat) +
  geom_line(aes(x=year,y=CoP,group=interaction(dotname,age_years,vaccine),color=factor(age_years))) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  ylab('VE') +
  xlab('years post-TCV') 


ggplot(plot_dat) +
  geom_line(aes(x=year,y=ve,group=interaction(dotname,age_years,vaccine),color=dotname,linetype=factor(age_years,levels=c(8,1.4)))) +
  geom_line(aes(x=year,y=mean,group=interaction(dotname,age_years,vaccine),color=dotname,linetype=factor(age_years,levels=c(8,1.4)))) + 
  facet_grid('~outcome') +
  theme_bw() +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
  scale_linetype_discrete(labels=c('5-15','0.75-2')) +
  ylab('VE') +
  xlab('years post-TCV') +
  labs(linetype='age')



#' Boosting model
#' 
#' Parameters taken from boosting section of titer_model_fit_scratch.R
#' 
fold_rise_model = function(CoP_pre,mu_0=3.3,CoP_max=10^3.51, CoP_min=1){
  fold_rise = pmax(1,10^(mu_0*(1-(log10(CoP_pre)-log10(CoP_min))/(log10(CoP_max)-log10(CoP_min)))))
  return(fold_rise)
}


boost_dat = data.frame(pre_vax_elisa = 10^seq(0,3.5,by=0.1)) |>
  mutate(fold_rise = fold_rise_model(CoP_pre = pre_vax_elisa)) |>
  mutate(post_vax_elisa = fold_rise*pre_vax_elisa)
ggplot(boost_dat) +
  geom_line(aes(x=pre_vax_elisa,y=fold_rise)) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  xlab('pre-vaccination IgG [EU/ml]') + ylab('fold-rise')

ggplot(boost_dat) +
  geom_line(aes(x=pre_vax_elisa,y=post_vax_elisa)) +
  theme_bw() +
  scale_y_continuous(trans='log10',limits=c(1e0,10^3.6)) +
  scale_x_continuous(trans='log10') +
  xlab('pre-vaccination IgG [EU/ml]') +
  ylab('post-vaccination IgG [EU/ml]') 

#' What we see is the little bit of data on boosting are consistent with, for Typbar-TCV (and Vi-rEPA, another conjugate vaccine),
#' is that the immunogenicity appears to be maxed out such that titers cannot reliably get (much) higher than achieved by first vaccination.
#' Boosting thus restores Vi titers, but doesn't seem to be able to significantly raise them beyond what one dose achieves.
#' This is in contrast with the less immunogenic Vi-polysaccharide vaccine, where boosting may lead to increases in peak titers since the first dose 
#' doesn't max titer out...
#' 
#' Also, we almost certainly need to think about age. Some of the odd behavior may be not including age properly...

# let's add a boost to the waning model above
# boost at 10 years
boost_year = 10

plot_dat = data.frame(t=seq(28,20*365,by=1),
                      dose = 1e3,
                      age_years=8,
                      dotname='Blantyre.Malawi',
                      vaccine='Typbar-TCV') |>
  mutate(CoP=titer_vs_time(t=t,C_max=2400,age_years=age_years)) |>
  mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
         outcome = 'fever given dose') |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 1e3,
                   age_years=1.4,
                   dotname='Blantyre.Malawi',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=2400,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 5e4,
                   age_years=8,
                   dotname='Mirpur.Dhaka.Bangladesh',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=2400,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 5e4,
                   age_years=1.4,
                   dotname='Mirpur.Dhaka.Bangladesh',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=2400,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 1e3,
                   age_years=8,
                   dotname='Blantyre.Malawi',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=2400,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 1e3,
                   age_years=1.4,
                   dotname='Blantyre.Malawi',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=2400,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 5e4,
                   age_years=8,
                   dotname='Mirpur.Dhaka.Bangladesh',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3220,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(28,20*365,by=1),
                   dose = 1e3,
                   age_years=1.4,
                   dotname='Mirpur.Dhaka.Bangladesh',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3220,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  mutate(year=t/365) |>
  group_by(ceiling(year),outcome,dotname,age_years,vaccine) |>
  mutate(mean=if_else(year<=5,mean(ve),NA)) |> 
  ungroup() 

# I should've called it a day two hours ago but here we are
vaccines=unique(plot_dat$vaccine)
ages=unique(plot_dat$age_years)
dotnames=unique(plot_dat$dotname)
outcomes = unique(plot_dat$outcome)
for (k in 1:length(vaccines)){
  for (n in 1:length(ages)){
    for (m in 1:length(outcomes)){
      for (q in 1:length(dotnames)){
        idx = plot_dat$vaccine[k]==vaccines[k] & plot_dat$age_years==ages[n] & plot_dat$outcome==outcomes[m] & plot_dat$dotname==dotnames[q] & plot_dat$t >= boost_year*365
        plot_dat$CoP[idx] = titer_vs_time(t=plot_dat$t[idx]-boost_year*365,
                                          C_max=fold_rise_model(CoP_pre = plot_dat$CoP[idx & plot_dat$t==boost_year*365])*plot_dat$CoP[idx & plot_dat$t==boost_year*365],
                                          age_years=ages[n]+boost_year,
                                          C_pre=plot_dat$CoP[idx & plot_dat$t==boost_year*365])
      }
    }
  }
}

plot_dat = plot_dat |> group_by(dotname, outcome, age_years,vaccine) |>
  mutate(ve = vaccine_efficacy(dose=unique(dose),CoP=CoP,outcome=gsub(' ','_',unique(outcome)))) |>
  mutate(age_label_coarse=if_else(age_years==1.4,'toddler','children')) |>
  mutate(age_label_coarse=factor(age_label_coarse,levels=c('toddler','children'))) |>
  mutate(mean = if_else(dotname=='Mirpur.Dhaka.Bangladesh' & year<=2,mean(ve[dotname=='Mirpur.Dhaka.Bangladesh' & year<=2 ]),mean)) |>
  mutate(mean = if_else(dotname=='Mirpur.Dhaka.Bangladesh' & year>2 & year <=5,mean(ve[dotname=='Mirpur.Dhaka.Bangladesh' & year>2 & year <=5]),mean))


ggplot(plot_dat) +
  geom_line(aes(x=year,y=ve,group=interaction(dotname,age_years,vaccine),color=dotname)) +
  geom_line(aes(x=year,y=mean,group=interaction(dotname,age_years,vaccine),color=dotname)) + 
  facet_grid('outcome~age_label_coarse') +
  theme_bw() +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
  ylab('VE') +
  xlab('years post-TCV')

plot_dat_2 = plot_efficacy_df |> 
  filter(dotname %in% c('Blantyre.Malawi','Mirpur.Dhaka.Bangladesh')) |>
  filter(age_label_coarse %in% c('toddler','children')) |>
  filter(estimator == 'IRR_poisson_regression') |>
  filter(endpoint_duration == 'interval') |>
  filter(age_cohort == 'total') |>
  mutate(outcome='fever given dose') 

ggplot() +
  geom_line(data=plot_dat,aes(x=year,y=ve,group=interaction(dotname,age_years,vaccine),color=dotname),size=1) +
  # geom_line(data=plot_dat,aes(x=year,y=mean,group=interaction(dotname,age_years,vaccine),color=dotname),size=1) +
  geom_segment(data=plot_dat_2,aes(x=timepoint_months_authors_mean/12,y=lower/100,yend=upper/100,color=dotname),size=0.5,alpha=0.5) +
  geom_point(data=plot_dat_2,aes(x=timepoint_months_authors_mean/12,y=estimate/100,color=dotname)) +
  facet_grid('outcome~age_label_coarse') +
  theme_bw() +
  scale_y_continuous(limits=c(-0.15,1),breaks=round(seq(-0.1,1,by=0.1),1)) +
  ylab('VE') +
  xlab('years post-TCV') +
  xlim(c(0,20))
ggsave('scratch/figures/model_efficacy_vs_time_location.png',units='in',width=7, height=5)

plot_dat = plot_dat |> group_by(age_years,vaccine,outcome) |>
  mutate(rr = (1-ve)/(1-max(ve[year==1])))
ggplot() +
  geom_line(data=plot_dat,aes(x=year,y=rr,group=interaction(dotname,age_years,vaccine),color=dotname),size=1) +
  # geom_line(data=plot_dat,aes(x=year,y=mean,group=interaction(dotname,age_years,vaccine),color=dotname),size=1) +
  # geom_segment(data=plot_dat_2,aes(x=timepoint_months_authors_mean/12,y=1-lower/100,yend=1-upper/100,color=dotname),size=0.5,alpha=0.5) +
  # geom_point(data=plot_dat_2,aes(x=timepoint_months_authors_mean/12,y=1-estimate/100,color=dotname)) +
  facet_grid('outcome~age_label_coarse',scale='free_y') +
  theme_bw() +
  # scale_y_continuous(breaks=round(seq(1,10,by=0.5),1)) +
  scale_x_continuous(breaks=seq(0,20,by=5)) +
  ylab('risk relative to 1 year post-vax') +
  xlab('years post-TCV') +
  xlim(c(0,20))
ggsave('scratch/figures/model_RR_vs_time_location.png',units='in',width=7, height=5)


# /* back matter for exporting as a blog post

# source('./docs/docs_helper_functions.R')
# render_blog_post(input_file = './scratch/empirical_bayes_random_effects_via_optimization.R',
#                  categories_list = c('Calibration','Self-study'),
#                  date_created = '2025-04-10')

# */

