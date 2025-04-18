# incidence model scratch

library(tidyverse)

# define cohort

N=1000
year_max=60
dt=365/12

population = expand.grid(id=1:N,
                         month=seq(0,year_max*12,by=1),
                         titer=1,
                         exposed=0,
                         infected=0,
                         fever=0,
                         exposure_month=0,
                         infection_month=0,
                         fever_month=0) |>
  mutate(year=month/12,
         age_group=cut(year,breaks=c(0,2,5,10,15,year_max),include.lowest = TRUE)) |>
  mutate(age_range = round(as.numeric(as.character(fct_recode( age_group, `2`='[0,2]',`3`='(2,5]',`5`='(5,10]',`5.1`='(10,15]',`45`='(15,60]')))))


# model equations
param_df = data.frame(value=c(T_decay=11,k=1.16,
                           beta_T_decay_age=0.057, beta_k_age=-0.060,
                           C_min=1,T_rise=1.5,t_start=3,t_peak=30.4), 
                   component = 'anti_Vi_IgG_elisa') 

titer_vs_time = function(t,C_max=10^3.5, age_years,params=param_df,C_pre=1){
  
  Cma=C_max
  Tda = params['T_decay','value'] *exp(params['beta_T_decay_age','value']  * age_years)
  ka = params['k','value']*exp(params['beta_k_age','value']  * age_years)
  
  # power law decay
  titer = params['C_min','value']  + (Cma-params['C_min','value'])*(1+(t-params['t_peak','value'] )/(ka*Tda))^(-ka)
  
  # simple logistic interpolator rise (this is just for continuity/realism. plays no role in the model)
  titer[t<params['t_peak','value'] ] =
    C_pre + (Cma-C_pre)*
    (1/(1+exp(-(t[t<params['t_peak','value'] ]-params['t_start','value'] *5)/params['T_rise','value'])) - 1/(1+exp(-(0-params['t_start','value']*5)/params['T_rise','value'])))/
    (1/(1+exp(-(params['t_peak','value'] -params['t_start','value']*5)/params['T_rise','value'])) - 1/(1+exp(-(0-params['t_start','value']*5)/params['T_rise','value'])))
  
  return(titer)
}

# natural immunity defaults
fold_rise_model = function(CoP_pre,mu_0=1.25,CoP_max=10^3.5, CoP_min=1){
  fold_rise = 10^(mu_0*(1-(log10(CoP_pre)-log10(CoP_min))/(log10(CoP_max)-log10(CoP_min))))
  return(fold_rise)
}

# dose response

# rescale variable to remove immune confound
gamma=0.4
# rescale alpha
alpha_prime=50^gamma*0.175

# rescale n50
N50_prime=1.11e6*(2^(1/alpha_prime)-1)/(2^(1/0.175)-1)
N50_prime

# note that this rescaling actually gets the control of the Jin2017 study decently right if bicarbonate right before challenge is worth ~a factor of 10 in dose
param_df= param_df |>
  rbind(data.frame(value=c(n50_fever_given_dose=N50_prime, 
                            alpha_fever_given_dose=alpha_prime,  
                            gamma_fever_given_dose=gamma), # guesstimated
                    component = 'p_fever_given_dose'))  |>
  
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


# run model

# very high
exposure_dose = 5e4
exposure_rate=1/(12*4) # per month

# high
exposure_dose = 5e3
exposure_rate=1/(12*4) # per month

# medium
exposure_dose = 5e2
exposure_rate=1/(12*15) # per month

exposure_rate_multiplier = c(rep(0.1,13),rep(0.5,12),rep(1,12*13),rep(1,year_max*12-24-12*13))

# expose
for (id in (1:N)){
  idx = population$id == id
  population$exposed[idx] = rpois(length(exposure_rate_multiplier),exposure_rate*exposure_rate_multiplier)
  population$exposure_month[idx & population$exposed>0] = population$month[idx & population$exposed>0]
}

# infect and titer
for (id in (1:N)){
  idx = population$id == id
  
  tmp_month=population$month[idx]
  tmp_titer =population$titer[idx]
  tmp_exposed = population$exposed[idx]
  tmp_infected = population$infected[idx]
  tmp_fever =population$fever[idx]
  
  for (exposure in which(tmp_exposed>0)){
    p_once = p_outcome_given_dose(dose=exposure_dose,CoP=tmp_titer[exposure],outcome='infection_given_dose')
    p_inf = 1-(1-p_once)^tmp_exposed[exposure]
    
    if(runif(1)<=p_inf){
      tmp_infected[exposure]=1
      tmp_titer[exposure:length(tmp_titer)]= titer_vs_time(t=(dt*(tmp_month-exposure+1))[(exposure):length(tmp_month)],
                                                           age_years = min(exposure/12,15),
                                                           C_pre = tmp_titer[exposure],
                                                           C_max=tmp_titer[exposure]*fold_rise_model(CoP_pre = tmp_titer[exposure]))
      tmp_titer[exposure:length(tmp_titer)] = round(tmp_titer[exposure:length(tmp_titer)],2)
    }
  }
  
  # fever
  infections = which(tmp_infected==1)
  p_fever = p_outcome_given_dose(dose=exposure_dose,CoP=tmp_titer[infections],outcome='fever_given_infection')
  fever_idx=runif(n=length(p_fever))<p_fever
  tmp_fever[infections[fever_idx]]=1
  
  population$titer[idx] = tmp_titer
  population$infected[idx] = tmp_infected
  population$fever[idx] = tmp_fever
  
  population$infection_month[idx & population$infected>0] = population$month[idx & population$infected>0]
  population$fever_month[idx & population$fever>0] = population$month[idx & population$fever>0]
  
}

ggplot(population |> filter(id<=10)) +
  geom_point(aes(x=exposure_month/12,y=id)) + xlim(c(1/12,70))

ggplot(population |> filter(id<=10)) +
  geom_point(aes(x=infection_month/12,y=id)) + xlim(c(1/12,70))

ggplot(population |> filter(id<=10)) +
  geom_point(aes(x=fever_month/12,y=id)) + xlim(c(1/12,70))

ggplot(population |> filter(id<=10)) +
  geom_line(aes(x=year,y=titer,color=id)) + xlim(c(1,70)) +
  facet_wrap('id')


# incidence targets
incidence_fever_targets = c(medium = 53,high=214,very_high=1255)
incidence_fever_targets

# incidence
incidence_vs_age =
  population |> group_by(age_group) |>
    summarize(incidence_fever = sum(fever==1)/(N*unique(age_range))*1e5,
              incidence_infection = sum(infected==1)/(N*unique(age_range))*1e5,
              symptomatic_fraction = sum(fever==1,na.rm=TRUE)/sum(infected==1,na.rm=TRUE))

ggplot(incidence_vs_age) +
  geom_bar(aes(x=age_group,y=incidence_fever),stat='identity') +
  ggtitle(paste('medium: dose = ',exposure_dose,', average years b/w exposures = ',1/(12*exposure_rate),sep=''))
ggsave('scratch/figures/cohort_model_incidence_by_age_medium.png',units='in',width=6,height=4)

# ggplot(incidence_vs_age) +
#   geom_bar(aes(x=age_group,y=incidence_infection),stat='identity')
# 
# 
# ggplot(incidence_vs_age) +
#   geom_bar(aes(x=age_group,y=symptomatic_fraction),stat='identity')
