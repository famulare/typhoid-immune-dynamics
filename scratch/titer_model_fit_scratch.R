# getting titer model going

library(tidyverse)
library(stats4)

# time series function for waning
# units in days
t_peak=21
t = seq(0,1e3,by=1)

titer_vs_time = function(t,T_decay=100,alpha=1,C_max=10,
                 C_min=1,T_rise=1.5,t_start=3,t_peak=28){
  # power law decay
  titer = (1+(t-t_peak)/(alpha*T_decay))^(-alpha)
  
  # simple logistic interpolator rise (this is just for continuity/realism. plays no role in the model)
  titer[t<t_peak] =
    (1/(1+exp(-(t[t<t_peak]-t_start*5)/T_rise)) - 1/(1+exp(-(0-t_start*5)/T_rise)))/
    (1/(1+exp(-(t_peak-t_start*5)/T_rise)) - 1/(1+exp(-(0-t_start*5)/T_rise)))
  
  # scale dimensions
  titer = C_min + (C_max-C_min)*titer
  
  return(titer)
}

plot_dat = data.frame(t=t,titer=titer_vs_time(t=t))
ggplot(plot_dat) +
  geom_line(aes(x=t,y=titer)) +
  theme_bw() +
  scale_y_continuous(limits=c(0,10),breaks=seq(0,10,by=1))



# bring in data!

# immunogenicitiy data
d_imm = readxl::read_excel('data/typhoid_vaccine_study_data.xlsx',sheet = 'immunogenicity') |>
  pivot_longer(
    cols = matches("^(treatment|control)_day"),    # select columns starting with "treatment" or "control" and "day"
    names_to = c("cohort", "day", "measure"),       # create three new columns from parts of the name
    names_pattern = "(treatment|control)_day(\\d+)_(.*)",  # regex to extract the cohort, day, and measure
    values_to = "value",                             # the data goes into a new column named 'value'
    values_transform = list(value = as.character)
  ) |>
  mutate(day = as.numeric(day)) |>
  pivot_wider(
    names_from = measure,
    values_from = value
  ) |>
  drop_na(N) |>
  mutate(across(c("age_min_years","age_max_years","N", "seropositive_percent", "fold_rise_median","fold_rise_lower",
                  "fold_rise_upper","seroconversion_percent"), as.numeric)) |>
  mutate(elisa_U_per_ml_median_numeric = if_else(elisa_U_per_ml_median=='<7.5',3.75,if_else(elisa_U_per_ml_median=='<3.0',1.5,as.numeric(elisa_U_per_ml_median)))) |>
  mutate(elisa_U_per_ml_lower_numeric = if_else(elisa_U_per_ml_lower=='<7.5',3.75,if_else(elisa_U_per_ml_lower=='<3.0',1.5,as.numeric(elisa_U_per_ml_lower)))) |>
  mutate(elisa_U_per_ml_upper_numeric = if_else(elisa_U_per_ml_upper=='<7.5',3.75,if_else(elisa_U_per_ml_upper=='<3.0',1.5,as.numeric(elisa_U_per_ml_upper)))) |>
  mutate(age_mean = 1/2*(age_min_years + age_max_years)) |>
  mutate(age_label = reorder(factor(paste(age_min_years,'-',age_max_years,sep='')),age_mean))

names(d_imm)


# need to rescale the Kossaczka1999 Vietnam data that uses a different elisa assay to match VaccZyme standard
# common reference is from Vi-polysaccharide vaccine

relevant_data = d_imm |> 
  filter(vaccine == 'Vi-polysaccharide') |>
  # comparable age cohorts
  filter(age_min_years==5) |>
  # comparable days
  filter(day ==0 | day == 42) |>
  # only comparable antigen left is IgG
  filter(antigen=='IgG') |>
  select(study,antigen, age_min_years,age_max_years,age_mean,day,elisa_U_per_ml_median_numeric)
relevant_data

# scaling multiplier thankfully looks homogenous from day 0 and day 42!
relevant_data |> group_by(age_min_years,day) |>
  reframe(ratio = elisa_U_per_ml_median_numeric/elisa_U_per_ml_median_numeric[study=='Kossaczka1999'])

Kossaczka_multiplier = (25.2 + 21.5)/2
Kossaczka_multiplier 

# rescale all Kossaczka data, including IgA (which is a leap!)

names(d_imm)
d_imm = d_imm |> 
  mutate(elisa_U_per_ml_median_numeric = if_else(study=='Kossaczka1999', Kossaczka_multiplier*elisa_U_per_ml_median_numeric,elisa_U_per_ml_median_numeric),
         elisa_U_per_ml_lower_numeric = if_else(study=='Kossaczka1999', Kossaczka_multiplier*elisa_U_per_ml_lower_numeric,elisa_U_per_ml_lower_numeric),
         elisa_U_per_ml_upper_numeric = if_else(study=='Kossaczka1999', Kossaczka_multiplier*elisa_U_per_ml_upper_numeric,elisa_U_per_ml_upper_numeric))

# plot total individual VE for bacteriologically-confirmed typhoid fever vs time
plot_dat = d_imm |>
  filter(cohort == 'treatment') 

ggplot(plot_dat, aes(shape=dotname,color=age_label,fill=age_label,group=age_label)) +
  geom_ribbon(aes(x=day,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric),color=NA,alpha=0.1) +
  geom_line(aes(x=day,y=elisa_U_per_ml_median_numeric)) +
  geom_point(aes(x=day,y=elisa_U_per_ml_median_numeric)) +
  facet_grid('antigen ~ vaccine') +
  theme_bw() +
  scale_y_continuous(trans='log10')


# fit to a single study group
fit_dat = d_imm |>
  filter(cohort == 'treatment') |>
  # filter(vaccine == 'Typbar-TCV') |> 
  filter(antigen == 'IgG') |>
  # filter(age_label %in% c('5-9','5-10')) |>
  filter(vaccine_schedule == 'single_dose') |>
  filter(day>0) |>
  filter(study !='Jin2017') |> filter(!(study=='Kossaczka1999' & age_min_years ==18 )) #excluding outliers for now...

fit_dat = fit_dat |> arrange(age_mean)

log_lsq_likelihood = function(C_max_TCV,T_decay,alpha,
                              beta_C_max_age_mean=0,beta_T_decay_age_mean=0,beta_alpha_age_mean=0,
                              Delta_C_max_Vips=0,Delta_C_max_VirEPA=0){
  
  titer=rep(NA,nrow(fit_dat))
  Tda = T_decay*exp(beta_T_decay_age_mean * fit_dat$age_mean)
  alpha_a = alpha*exp(beta_alpha_age_mean * fit_dat$age_mean)
  Cma=C_max_TCV*exp(beta_C_max_age_mean * fit_dat$age_mean)
 
   for (row in 1:length(titer)){
     if (fit_dat$vaccine[row] == 'Vi-polysaccharide'){
       Cma[row] = Cma[row]*exp(Delta_C_max_Vips)
     } else if (fit_dat$vaccine[row] == 'Vi-rEPA2'){
       Cma[row] = Cma[row]*exp(Delta_C_max_VirEPA)
     }
     
    titer[row]=titer_vs_time(t=fit_dat$day[row],C_max=Cma[row],T_decay=Tda[row],alpha=alpha_a[row])
  }
  
  # quick and dirty mean to mean log-lsq for now
    mlogL = sum((log(fit_dat$elisa_U_per_ml_median_numeric)-log(titer))^2/(log(fit_dat$elisa_U_per_ml_upper_numeric)-log(fit_dat$elisa_U_per_ml_lower_numeric))^2)
  
  return(mlogL)
}


mod = mle(log_lsq_likelihood,
          start=list(C_max_TCV=4000,T_decay=100,alpha=1,beta_C_max_age_mean=0,beta_T_decay_age_mean=0,beta_alpha_age_mean=0,Delta_C_max_Vips=-0,Delta_C_max_VirEPA=0),
          fixed=list(),
          control=list(parscale=c(1000,10,0.05,0.001,0.001,0.001,1,1)))
summary(mod) 

# fit
fitted= expand.grid(day=seq(0,2000,by=10),
                    age_mean=unique(fit_dat$age_mean),
                    antigen='IgG',
                    vaccine = unique(fit_dat$vaccine),
                    titer=NA) 
age_means=unique(fitted$age_mean)
vaccines=unique(fit_dat$vaccine)
Delta_coeff_C = c(0, coef(mod)['Delta_C_max_VirEPA'],coef(mod)['Delta_C_max_Vips'])
for (k in 1:length(age_means)){
  for (n in 1:length(vaccines)){
    idx = fitted$age_mean==age_means[k] & fitted$vaccine==vaccines[n]
    fitted$titer[idx] = titer_vs_time( t=seq(0,2000,by=10),
                                C_max=coef(mod)[1]*exp(coef(mod)[4]*age_means[k])*exp(Delta_coeff_C[n]),
                                T_decay=coef(mod)[2]*exp(coef(mod)[5]*age_means[k]),
                                alpha=coef(mod)[3]*exp(coef(mod)[6]*age_means[k]))
  }
}
fitted = fitted |> filter(interaction(age_mean,vaccine) %in% unique(interaction(fit_dat$age_mean,fit_dat$vaccine)))

ggplot() +
  geom_ribbon(data=fit_dat,aes(x=day,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric,
                                fill=age_label,group=age_label),alpha=0.1) +
  geom_line(data=fit_dat,aes(x=day,y=elisa_U_per_ml_median_numeric,
                              color=age_label,group=age_label)) +
  geom_point(data=fit_dat,aes(x=day,y=elisa_U_per_ml_median_numeric,
                               shape=dotname,color=age_label,group=age_label)) +
  geom_line(data=fitted,aes(x=day,y=titer,group=interaction(age_mean,vaccine)),color='black') +
  facet_wrap('vaccine ~ age_mean') +
  theme_bw() +
  scale_y_continuous(trans='log10')



## boosting

names(d_imm)

boosting_data = d_imm |> filter(study %in% c('Kossaczka1999','Mohan2015')) |>
  filter(antigen=='IgG') |>
  filter( (study == 'Kossaczka1999' & age_min_years == 2) | 
            (study == 'Mohan2015' & age_min_years == 2 & age_max_years == 45) | 
            (study == 'Mohan2015' & age_min_years == 0.5 & age_max_years == 2)) |>
  mutate(group_label = interaction(vaccine,age_label, vaccine_schedule))
boosting_data


ggplot(boosting_data) +
  geom_ribbon(aes(x=day,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric,
                               fill=vaccine_schedule,group=group_label),alpha=0.1) +
  geom_line(aes(x=day,y=elisa_U_per_ml_median_numeric,
                             color=vaccine_schedule,group=group_label)) +
  geom_point(aes(x=day,y=elisa_U_per_ml_median_numeric,
                              shape=study,color=vaccine_schedule,group=group_label)) +
  facet_wrap('vaccine ~ age_mean') +
  theme_bw() +
  scale_y_continuous(trans='log10')


# slice it into fold-rise
# first, exclude the time points that aren't immediately pre or post vaccination
boosting_data = boosting_data |> 
  filter(day <100) |>
  filter(!(study == 'Kossaczka1999' & vaccine_schedule == 'single_dose' & day==70))

# all remaining data is a pre-post pair or pre-(post1=pre2)-post2 triplet
ggplot(boosting_data) +
  geom_ribbon(aes(x=day,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric,
                  fill=vaccine_schedule,group=group_label),alpha=0.1) +
  geom_line(aes(x=day,y=elisa_U_per_ml_median_numeric,
                color=vaccine_schedule,group=group_label)) +
  geom_point(aes(x=day,y=elisa_U_per_ml_median_numeric,
                 shape=study,color=vaccine_schedule,group=group_label)) +
  facet_wrap('vaccine ~ age_mean') +
  theme_bw() +
  scale_y_continuous(trans='log10')


# let's duplicate and shift that triplet so everything is a pre-post pair
tmp = boosting_data |> filter(study == 'Kossaczka1999' & vaccine_schedule == 'two_dose_6_week') 
tmp$vaccine_schedule[1:2]='single_dose'
tmp[4,] = tmp[3,]
tmp[3,] = tmp[2,]
tmp[3,'day']=0
tmp[4,'day']=70-42
tmp[3:4,'vaccine_schedule']='booster'

boosting_data = boosting_data |>
  filter(!(study == 'Kossaczka1999' & vaccine_schedule == 'two_dose_6_week')) |>
  rbind(tmp) |>
  mutate(group_label = interaction(study,vaccine,age_label, vaccine_schedule))


# all remaining data is a pre-post pair!
ggplot(boosting_data) +
  geom_ribbon(aes(x=day,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric,
                  fill=vaccine_schedule,group=group_label),alpha=0.1) +
  geom_line(aes(x=day,y=elisa_U_per_ml_median_numeric,
                color=vaccine_schedule,group=group_label)) +
  geom_point(aes(x=day,y=elisa_U_per_ml_median_numeric,
                 shape=study,color=vaccine_schedule,group=group_label)) +
  facet_wrap('vaccine ~ age_mean') +
  theme_bw() +
  scale_y_continuous(trans='log10')

# convert to fold rise
boosting_data = boosting_data |>
  select(study, vaccine,age_label, vaccine_schedule,antigen, age_min_years,age_max_years,age_mean,age_label, 
         day,elisa_U_per_ml_median_numeric,group_label) 

boosting_data$fold_rise[boosting_data$day>0] = 
  boosting_data$elisa_U_per_ml_median_numeric[boosting_data$day>0]/boosting_data$elisa_U_per_ml_median_numeric[boosting_data$day==0]

boosting_data$pre_vax_elisa[boosting_data$day>0] = boosting_data$elisa_U_per_ml_median_numeric[boosting_data$day==0] 

boosting_data = boosting_data |> filter(day>0)
  
# fold-rise vs pre-challenge correlate. It's a joy when the theory holds for a new system!
# https://famulare.github.io/2024/03/18/Hypothesis-why-do-neutralizing-antibody-titers-max-out.html
ggplot(boosting_data) +
  geom_point(aes(x=pre_vax_elisa,y=fold_rise,
                 shape=age_label,color=vaccine,group=group_label)) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')

# adjust fold-rise with waning model to get true peak fold rise estimates
# probably need to think about age!

boosting_data$fold_rise_adjusted = boosting_data$fold_rise

for (row in 1:nrow(boosting_data)){
  boosting_data$fold_rise_adjusted[row] = boosting_data$fold_rise_adjusted[row]*
    titer_vs_time(t=28,
                  C_max=coef(mod)[1]*exp(coef(mod)[4]*age_means[k])*exp(Delta_coeff_C[n]),
                  T_decay=coef(mod)[2]*exp(coef(mod)[5]*age_means[k]),
                  alpha=coef(mod)[3]*exp(coef(mod)[6]*age_means[k])) /
    titer_vs_time(t=boosting_data$day[row],
                  C_max=coef(mod)[1]*exp(coef(mod)[4]*age_means[k])*exp(Delta_coeff_C[n]),
                  T_decay=coef(mod)[2]*exp(coef(mod)[5]*age_means[k]),
                  alpha=coef(mod)[3]*exp(coef(mod)[6]*age_means[k]))
}

plot(boosting_data$fold_rise,boosting_data$fold_rise_adjusted-boosting_data$fold_rise)

# boosting model equations and fitting
fold_rise_model = function(CoP_pre,mu_0,CoP_max, CoP_min=1){
  fold_rise = 10^(mu_0*(1-(log10(CoP_pre)-log10(CoP_min))/(log10(CoP_max)-log10(CoP_min))))
  return(fold_rise)
}


fold_rise_mlogL = function(mu_TCV,mu_Vips,mu_VirEPA,CoP_max){
  mlogL=0
  for (row in 1:nrow(boosting_data)){
    if (boosting_data$vaccine[row]=='Typbar-TCV'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_TCV,CoP_max=CoP_max)
    } else if (boosting_data$vaccine[row]=='Vi-polysaccharide'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_Vips,CoP_max=CoP_max)
    } else if (boosting_data$vaccine[row]=='Vi-rEPA2'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_VirEPA,CoP_max=CoP_max)
    }
    fold_rise = fold_rise
  mlogL = mlogL + (log(boosting_data$fold_rise_adjusted[row])-log(fold_rise))^2
  }
  return(mlogL)
}

mod = mle(fold_rise_mlogL,
          start=list(mu_TCV=3,mu_Vips=2,mu_VirEPA=3,CoP_max=2000),
          control=list(parscale=c(0.1,0.1,0.1,10)))
summary(mod) 

vaccines = unique(boosting_data$vaccine)
fitted_boosting = expand.grid(vaccine=vaccines,
                              pre_vax_elisa = 10^seq(0,4,by=0.01))
for (row in 1:nrow(fitted_boosting)){
  if (fitted_boosting$vaccine[row]=='Typbar-TCV'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=coef(mod)[1],CoP_max=coef(mod)[4])
  } else if (fitted_boosting$vaccine[row]=='Vi-polysaccharide'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=coef(mod)[2],CoP_max=coef(mod)[4])
  } else if (fitted_boosting$vaccine[row]=='Vi-rEPA2'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=coef(mod)[3],CoP_max=coef(mod)[4])
  }
}

# yay! this works, except that there are many observed values above the so-called maximum. 
ggplot() +
  geom_line(data=fitted_boosting,aes(x=pre_vax_elisa,y=fold_rise,
                                     color=vaccine,group=vaccine)) +
  geom_point(data=boosting_data,aes(x=pre_vax_elisa,y=fold_rise_adjusted,
                 shape=age_label,color=vaccine,group=group_label)) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')


# let's see what happens if we just force the intercept
max(d_imm$elisa_U_per_ml_median_numeric)
max(d_imm$elisa_U_per_ml_upper_numeric)

mod = mle(fold_rise_mlogL,
          start=list(mu_TCV=3,mu_Vips=2,mu_VirEPA=3),
          fixed=list(CoP_max=20e3),
          control=list(parscale=c(0.1,0.1,0.1)))
summary(mod) 

vaccines = unique(boosting_data$vaccine)
fitted_boosting = expand.grid(vaccine=vaccines,
                              pre_vax_elisa = 10^seq(0,4.3,by=0.01))
for (row in 1:nrow(fitted_boosting)){
  if (fitted_boosting$vaccine[row]=='Typbar-TCV'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=coef(mod)[1],CoP_max=20e3)
  } else if (fitted_boosting$vaccine[row]=='Vi-polysaccharide'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=coef(mod)[2],CoP_max=20e3)
  } else if (fitted_boosting$vaccine[row]=='Vi-rEPA2'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=coef(mod)[3],CoP_max=20e3)
  }
}

# yay! this works, except that there are many observed values above the so-called maximum. 
ggplot() +
  geom_line(data=fitted_boosting,aes(x=pre_vax_elisa,y=fold_rise,
                                     color=vaccine,group=vaccine)) +
  geom_point(data=boosting_data,aes(x=pre_vax_elisa,y=fold_rise_adjusted,
                                    shape=age_label,color=vaccine,group=group_label)) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')

# hmm... it is interesting that if, maybe, you happen to get to a higher titer than the boost max, that we wouldn't expect a boost again...
# Ah I see what's happening! It looks like the conjugate vaccines hit the maximum possible antibody level in one dose. 
# And so each re-dose just brings you back up.
# 
# It's also interesting that the polysaccharide data looks like it has a shared slope with the conjugates and not a shared intercept.
# My expectation is the intercept should be a human host property, not a vaccine property. But perhaps this is wrong. Not much to say confidently
# either way with just 2 points, but I'm sure there's more Vi-polysaccharide vaccine data out there to add...

