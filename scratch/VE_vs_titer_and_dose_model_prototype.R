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
library(stats4)

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
  gamma=0.3 #0.25
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
                   CoP=563) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=563) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=563) |>
          mutate(p = p_outcome_given_dose(dose,CoP=CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |> 
  mutate(CoP=factor(CoP))

ggplot(plot_dat ) +
  geom_line(aes(x=dose,y=p,group=outcome,color=outcome)) +
  facet_wrap('CoP') +
  theme_bw() +
  scale_x_continuous(trans='log10', breaks=10^seq(1,10,by=2),minor_breaks = NULL) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
  ylab('probability of outcome') +
  xlab('dose [bacilli]')


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
                   CoP=563) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=3300) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |> 
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
             CoP=50) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=563) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=3300) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=50) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=563) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |>
  rbind(data.frame(dose = 10^seq(0,10,by=0.05),
                   CoP=3300) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_infection'),
                 outcome = 'fever given infection')) |>
  mutate(CoP=factor(CoP))

ggplot(plot_dat) +
  geom_line(aes(x=dose,y=ve,group=CoP,color=CoP)) +
  facet_grid('~outcome') +
  theme_bw() +
  scale_x_continuous(trans='log10', breaks=10^seq(1,10,by=2),minor_breaks = NULL) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
  ylab('VE') +
  xlab('dose [bacilli]')

#' This kinda looks like it works!

#' # Waning model!
#' 

# time series function for waning
# units in days
t_peak=21
t = seq(0,1e3,by=1)

param_df = param_df |>
  rbind(data.frame(value=c(T_decay=9.3,k=0.98,
                           beta_C_max_age_mean=-0.029, beta_T_decay_age_mean=0.15, beta_k_age_mean=-0.023,
                           C_min=1,T_rise=1.5,t_start=3,t_peak=28), 
                   component = 'anti_Vi_IgG_elisa') )

titer_vs_time = function(t,C_max, age_years,params=param_df){
  
  Cma=C_max*(1+params['beta_C_max_age_mean','value'] * age_years)
  Tda = params['T_decay','value'] *(1+params['beta_T_decay_age_mean','value']  * age_years)
  ka = params['k','value']*(1+params['beta_k_age_mean','value']  * age_years)
  
  # power law decay
  titer = (1+(t-params['t_peak','value'] )/(ka*Tda))^(-ka)
  
  # simple logistic interpolator rise (this is just for continuity/realism. plays no role in the model)
  titer[t<params['t_peak','value'] ] =
    (1/(1+exp(-(t[t<params['t_peak','value'] ]-params['t_start','value'] *5)/params['T_rise','value'])) - 1/(1+exp(-(0-params['t_start','value']*5)/params['T_rise','value'])))/
    (1/(1+exp(-(params['t_peak','value'] -params['t_start','value']*5)/params['T_rise','value'])) - 1/(1+exp(-(0-params['t_start','value']*5)/params['T_rise','value'])))
  
  # scale dimensions
  titer = params['C_min','value'] + (Cma-params['C_min','value'])*titer
  
  return(titer)
}

plot_dat = data.frame(t=seq(1,20*365,by=1),
                      dose = 5e1,
                      age_years=8,
                      location='Malawi-like',
                      vaccine='Typbar-TCV') |>
  mutate(CoP=titer_vs_time(t=t,C_max=3650,age_years=age_years)) |>
  mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
         outcome = 'fever given dose') |>
  rbind(data.frame(t=seq(1,20*365,by=1),
                   dose = 1e3,
                   age_years=1.4,
                   location='Malawi-like',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3650,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(1,20*365,by=1),
                   dose = 5e4,
                   age_years=8,
                   location='Dhaka-like',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3650,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(1,20*365,by=1),
                   dose = 1e6,
                   age_years=1.4,
                   location='Dhaka-like',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3650,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'fever_given_dose'),
                 outcome = 'fever given dose')) |>
  rbind(data.frame(t=seq(1,20*365,by=1),
                     dose = 5e1,
                     age_years=8,
                     location='Malawi-like',
                     vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3650,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(1,20*365,by=1),
                   dose = 1e3,
                   age_years=1.4,
                   location='Malawi-like',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3650,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(1,20*365,by=1),
                   dose = 5e4,
                   age_years=8,
                   location='Dhaka-like',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3650,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  rbind(data.frame(t=seq(1,20*365,by=1),
                   dose = 1e6,
                   age_years=1.4,
                   location='Dhaka-like',
                   vaccine='Typbar-TCV') |>
          mutate(CoP=titer_vs_time(t=t,C_max=3650,age_years=age_years)) |>
          mutate(ve = vaccine_efficacy(dose,CoP,outcome = 'infection_given_dose'),
                 outcome = 'infection given dose')) |>
  mutate(year=t/365) |>
  mutate(year_groups = if_else(year<=2,'0-2',if_else(year<=5,'3-5','>5'))) |>
  group_by(year_groups,outcome,location,age_years,vaccine) |>
  mutate(mean=if_else(year<=5,mean(ve),NA)) |>
  ungroup()

ggplot(plot_dat) +
  geom_line(aes(x=year,y=ve,group=interaction(location,age_years,vaccine),color=location)) +
  geom_line(aes(x=year,y=mean,group=interaction(location,age_years,vaccine),color=location)) + 
  facet_grid('~outcome') +
  theme_bw() +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1)) +
  ylab('VE') +
  xlab('dose [bacilli]')


# /* back matter for exporting as a blog post

# source('./docs/docs_helper_functions.R')
# render_blog_post(input_file = './scratch/empirical_bayes_random_effects_via_optimization.R',
#                  categories_list = c('Calibration','Self-study'),
#                  date_created = '2025-04-10')

# */