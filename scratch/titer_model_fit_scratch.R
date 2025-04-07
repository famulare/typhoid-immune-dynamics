# getting titer model going

library(tidyverse)
library(stat4)

# time series function for waning
# units in days
t_peak=21
t = seq(0,1e3,by=1)

titer_vs_time = function(t,T_decay=100,alpha=1,C_max=10,
                 C_min=1,T_rise=1.5,t_start=2,t_peak=21){
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
  filter(vaccine == 'Typbar-TCV') |> 
  filter(antigen == 'IgG') |>
  # filter(age_label %in% c('5-9','5-10')) |>
  filter(vaccine_schedule == 'single_dose') |>
  filter(day>0)

log_lsq_likelihood = function(C_max,T_decay,alpha,beta_C_max_age_mean=0,beta_T_decay_age_mean=0,beta_alpha_age_mean=0){
  
  titer=rep(NA,nrow(fit_dat))
  Cma=C_max*(1+beta_C_max_age_mean * (fit_dat$age_mean-fit_dat$age_mean[1]))
  Tda = T_decay*(1+beta_T_decay_age_mean * (fit_dat$age_mean-fit_dat$age_mean[1]))
  alpha_a = alpha*(1+beta_alpha_age_mean * (fit_dat$age_mean-fit_dat$age_mean[1]))
 
   for (row in 1:length(titer)){
    titer[row]=titer_vs_time(t=fit_dat$day[row],C_max=Cma[row],T_decay=Tda[row],alpha=alpha_a[row])
  }
  
  # quick and dirty mean to mean log-lsq for now
    mlogL = sum((log(fit_dat$elisa_U_per_ml_median_numeric)-log(titer))^2/(log(fit_dat$elisa_U_per_ml_upper_numeric)-log(fit_dat$elisa_U_per_ml_lower_numeric))^2)
  
  return(mlogL)
}


mod = mle(log_lsq_likelihood,
          start=list(C_max=4000,T_decay=100,alpha=1,beta_C_max_age_mean=0,beta_T_decay_age_mean=0,beta_alpha_age_mean=0),
          fixed=list(),
          control=list(parscale=c(1000,10,0.05,0.01,0.01,0.01)))
summary(mod) 

coef(mod)[1]*(1+coef(mod)[4]*(fit_dat$age_mean-fit_dat$age_mean[1]))
coef(mod)[2]*(1+coef(mod)[5]*(fit_dat$age_mean-fit_dat$age_mean[1]))
coef(mod)[3]*(1+coef(mod)[6]*(fit_dat$age_mean-fit_dat$age_mean[1]))
view(fit_dat)

# fit
fitted= expand.grid(day=seq(0,2000,by=10),
                    age_mean=unique(fit_dat$age_mean),
                    antigen='IgG',
                    vaccine = unique(fit_dat$vaccine),
                    titer=NA) 
age_means=unique(fitted$age_mean)
for (k in 1:length(age_means)){
  idx = fitted$age_mean==age_means[k]
  fitted$titer[idx] = titer_vs_time( t=seq(0,2000,by=10),
                              C_max=coef(mod)[1]*(1+coef(mod)[4]*(age_means[k]-fit_dat$age_mean[1])),
                              T_decay=coef(mod)[2]*(1+coef(mod)[5]*(age_means[k]-fit_dat$age_mean[1])),
                              alpha=coef(mod)[3]*(1+coef(mod)[6]*(age_means[k]-fit_dat$age_mean[1])))
}


ggplot() +
  geom_ribbon(data=fit_dat,aes(x=day,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric,
                                fill=age_label,group=age_label),alpha=0.1) +
  geom_line(data=fit_dat,aes(x=day,y=elisa_U_per_ml_median_numeric,
                              color=age_label,group=age_label)) +
  geom_point(data=fit_dat,aes(x=day,y=elisa_U_per_ml_median_numeric,
                               shape=dotname,color=age_label,group=age_label)) +
  geom_line(data=fitted,aes(x=day,y=titer,group=age_mean),color='black') +
  facet_wrap('age_mean') +
  theme_bw() +
  scale_y_continuous(trans='log10') 

