# getting titer model going

library(tidyverse)
library(stats4)

# time series function for waning
# units in days
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
  mutate(age_label = reorder(factor(paste(age_min_years,'-',age_max_years,sep='')),age_mean)) |>
  mutate(age_label_coarse=cut(age_mean,breaks=c(0,1,2.5,7.5,15,40),
                              labels=c('infant','toddler','younger child','children','adult'),ordered_result=TRUE)) |>
  mutate(assay = factor(assay)) |>
  mutate(assay_simple=fct_recode(assay,VaccZyme = 'anti_Vi_IgG_VaccZyme',
                                 VaccZyme = 'anti_Vi_IgA_VaccZyme',
                                 Jackson_ImmunoResearch = 'anti_Vi_IgG_Jackson_ImmunoResearch_Laboratories',
                                 George_Carlone_CDC='anti_Vi_IgA_HP6107_George_Carlone_CDC')) |>
  mutate(assay_simple=fct_relevel(assay_simple,c('VaccZyme','Jackson_ImmunoResearch','George_Carlone_CDC')))

sort(unique(d_imm$age_mean))
d_imm$assay_simple

d_imm$age_label_coarse

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

d_imm |> 
  filter(study == 'Kossaczka1999') |>
  filter(antigen =='IgG') |>
  select(study,antigen,age_min_years, age_max_years, age_mean,day, elisa_U_per_ml_median_numeric) |>
  print(n=35)


# rescale all Kossaczka data, including IgA (which is a leap!)

names(d_imm)
d_imm = d_imm |> 
  mutate(elisa_U_per_ml_median_numeric = if_else(study=='Kossaczka1999', Kossaczka_multiplier*elisa_U_per_ml_median_numeric,elisa_U_per_ml_median_numeric),
         elisa_U_per_ml_lower_numeric = if_else(study=='Kossaczka1999', Kossaczka_multiplier*elisa_U_per_ml_lower_numeric,elisa_U_per_ml_lower_numeric),
         elisa_U_per_ml_upper_numeric = if_else(study=='Kossaczka1999', Kossaczka_multiplier*elisa_U_per_ml_upper_numeric,elisa_U_per_ml_upper_numeric))

# linear proportional scaling seems unlikely, with the very high 18-35 day 42 titer...
d_imm |> 
  filter(study == 'Kossaczka1999') |>
  filter(antigen =='IgG') |>
  select(study,day,vaccine,age_min_years, age_max_years, age_mean,day, elisa_U_per_ml_median,elisa_U_per_ml_median_numeric) |>
  print(n=35)


# plot total individual VE for bacteriologically-confirmed typhoid fever vs time
plot_dat = d_imm |>
  filter(cohort == 'treatment') 

ggplot(plot_dat, aes(shape=dotname,color=age_label_coarse,fill=age_label_coarse,group=age_label,linetype=assay_simple)) +
  geom_ribbon(aes(x=day/365,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric),color=NA,alpha=0.1) +
  geom_line(aes(x=day/365,y=elisa_U_per_ml_median_numeric)) +
  geom_point(aes(x=day/365,y=elisa_U_per_ml_median_numeric)) +
  facet_grid('antigen ~ vaccine') +
  theme_bw() +
  scale_y_continuous(trans='log10',limits=c(0.4,12e3)) +
  scale_color_brewer(palette='Paired') +
  scale_fill_brewer(palette='Paired') +
  xlab('years since vaccination') +
  ylab('anti-Vi IgG [EU/ml]')
ggsave('scratch/figures/titer_data.png',units='in',width=9, height=5)

ggplot(plot_dat, aes(color=dotname,shape=age_label_coarse,fill=dotname,group=interaction(dotname,age_label),linetype=assay_simple)) +
  geom_ribbon(aes(x=day/365,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric),color=NA,alpha=0.1) +
  geom_line(aes(x=day/365,y=elisa_U_per_ml_median_numeric)) +
  geom_point(aes(x=day/365,y=elisa_U_per_ml_median_numeric)) +
  facet_grid('antigen ~ vaccine') +
  theme_bw() +
  scale_y_continuous(trans='log10',limits=c(0.4,1.2e4)) +
  xlab('years since vaccination') +
  ylab('anti-Vi IgG [EU/ml]')
ggsave('scratch/figures/titer_data_location.png',units='in',width=9, height=5)



# fit 
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
                              beta_C_max_age=0,beta_T_decay_age=0,beta_alpha_age=0,
                              Delta_C_max_Vips=0,Delta_C_max_VirEPA=0,tp=30.4,g=1){
  
  # use age_mean
  titer=rep(NA,nrow(fit_dat))
  # if(age_method=='age_mean'){
  #   Tda = T_decay*exp(beta_T_decay_age * fit_dat$age_mean)
  #   alpha_a = alpha*exp(beta_alpha_age * fit_dat$age_mean)
  #   Cma=C_max_TCV*exp(beta_C_max_age * fit_dat$age_mean)
  #  
  #    for (row in 1:length(titer)){
  #      if (fit_dat$vaccine[row] == 'Vi-polysaccharide'){
  #        Cma[row] = Cma[row]*exp(Delta_C_max_Vips)
  #      } else if (fit_dat$vaccine[row] == 'Vi-rEPA2'){
  #        Cma[row] = Cma[row]*exp(Delta_C_max_VirEPA)
  #      }
  #      
  #     titer[row]=titer_vs_time(t=fit_dat$day[row],C_max=Cma[row],T_decay=Tda[row],alpha=alpha_a[row])
  #   }
  # } else if (age_method=='endpoint'){
    Tda_lower = T_decay*1/(1+exp(beta_T_decay_age * fit_dat$age_min_years))
    alpha_a_lower = alpha*(exp(beta_alpha_age * (fit_dat$age_min_years^exp(g))))
    Cma_lower=C_max_TCV*exp(beta_C_max_age * fit_dat$age_min_years)
    
    Tda_upper = T_decay*1/(1+exp(beta_T_decay_age * fit_dat$age_max_years))
    alpha_a_upper = alpha*(exp(beta_alpha_age * (fit_dat$age_max_years^exp(g))))
    Cma_upper=C_max_TCV*exp(beta_C_max_age * fit_dat$age_max_years)
    
    for (row in 1:length(titer)){
      if (fit_dat$vaccine[row] == 'Vi-polysaccharide'){
        Cma_lower[row] = Cma_lower[row]*exp(Delta_C_max_Vips)
        Cma_upper[row] = Cma_upper[row]*exp(Delta_C_max_Vips)
      } else if (fit_dat$vaccine[row] == 'Vi-rEPA2'){
        Cma_lower[row] = Cma_lower[row]*exp(Delta_C_max_VirEPA)
        Cma_upper[row] = Cma_upper[row]*exp(Delta_C_max_VirEPA)
      }
      
      titer[row]=1/2*(titer_vs_time(t=fit_dat$day[row],C_max=Cma_lower[row],T_decay=Tda_lower[row],alpha=alpha_a_lower[row],t_peak=tp) +
                        titer_vs_time(t=fit_dat$day[row],C_max=Cma_upper[row],T_decay=Tda_upper[row],alpha=alpha_a_upper[row],t_peak=tp))
    }
  # }
  
  # quick and dirty mean to mean log-lsq for now
    mlogL = sum((log(fit_dat$elisa_U_per_ml_median_numeric)-log(titer))^2/(log(fit_dat$elisa_U_per_ml_upper_numeric)-log(fit_dat$elisa_U_per_ml_lower_numeric))^2)
  
  return(mlogL)
}


mod = mle(log_lsq_likelihood,
          start=list(C_max_TCV=4000,T_decay=10,alpha=1,Delta_C_max_Vips=-0,Delta_C_max_VirEPA=0,
                     beta_C_max_age=0,beta_T_decay_age=0,beta_alpha_age=0,g=0),
          fixed=list(tp=30.4),
          control=list(parscale=c(1000,1,0.05,1,1,0.1,0.001,0.001,0.001))) 
summary(mod)
data.frame(summary(mod)@coef) |> mutate(z = Estimate / `Std..Error`)

cov2cor(vcov(mod))

model_coeffs=coef(mod)

# fit
fitted= expand.grid(day=seq(0,2000,by=10),
                    age_label=unique(fit_dat$age_label),
                    antigen='IgG',
                    vaccine = unique(fit_dat$vaccine),
                    titer=NA) |>
  mutate(age_min_years = as.numeric(sub('-.*','',age_label))) |>
  mutate(age_max_years = as.numeric(sub('.*-','',age_label)))
age_labels=unique(fit_dat$age_label)
vaccines=unique(fit_dat$vaccine)
Delta_coeff_C = c(0,model_coeffs['Delta_C_max_VirEPA'],model_coeffs['Delta_C_max_Vips'])
for (k in 1:length(age_labels)){
  for (n in 1:length(vaccines)){
    idx = fitted$age_label==age_labels[k] & fitted$vaccine==vaccines[n]
    # fitted$titer[idx] = titer_vs_time( t=seq(0,2000,by=10),
    #                             C_max=model_coeffs['C_max_TCV']*exp(model_coeffs['beta_C_max_age']*age_means[k])*exp(Delta_coeff_C[n]),
    #                             T_decay=model_coeffs['T_decay']*exp(model_coeffs['beta_T_decay_age']*age_means[k]),
    #                             alpha=model_coeffs['alpha']*exp(model_coeffs['beta_alpha_age']*age_means[k]))
    age_min=unique(fitted$age_min_years[idx])
    age_max=unique(fitted$age_max_years[idx])
    fitted$titer[idx] = 1/2*(titer_vs_time( t=seq(0,2000,by=10),
                                       C_max=model_coeffs['C_max_TCV']*exp(model_coeffs['beta_C_max_age']*age_min)*exp(Delta_coeff_C[n]),
                                       T_decay=model_coeffs['T_decay']*1/(1+exp(model_coeffs['beta_T_decay_age']*age_min)),
                                       alpha=model_coeffs['alpha']*(exp(model_coeffs['beta_alpha_age']*(age_min^exp(model_coeffs['g'])))),
                                       t_peak=model_coeffs['tp']) +
                               titer_vs_time( t=seq(0,2000,by=10),
                                              C_max=model_coeffs['C_max_TCV']*exp(model_coeffs['beta_C_max_age']*age_max)*exp(Delta_coeff_C[n]),
                                              T_decay=model_coeffs['T_decay']*1/(1+exp(model_coeffs['beta_T_decay_age']*age_max)),
                                              alpha=model_coeffs['alpha']*(exp(model_coeffs['beta_alpha_age']*(age_max^exp(model_coeffs['g'])))),
                                              t_peak=model_coeffs['tp']))
  }
}
fitted = fitted |> filter(as.character(interaction(age_label,vaccine)) %in% as.character(unique(interaction(fit_dat$age_label,fit_dat$vaccine)))) |>
  mutate(age_label=factor(age_label,levels=levels(d_imm$age_label)))

ggplot() +
  geom_line(data=fitted,aes(x=day/365,y=titer,group=interaction(age_label,vaccine)),color='black') +
  geom_point(data=fit_dat,aes(x=day/365,y=elisa_U_per_ml_median_numeric,
  shape=dotname,color=age_label,group=age_label)) +
  geom_ribbon(data=fit_dat,aes(x=day/365,ymin=elisa_U_per_ml_lower_numeric,ymax=elisa_U_per_ml_upper_numeric,
                                fill=age_label,group=age_label),alpha=0.3) +
  geom_line(data=fit_dat,aes(x=day/365,y=elisa_U_per_ml_median_numeric,
  color=age_label,group=age_label)) +
    facet_wrap('vaccine ~ age_label',ncol=8) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  guides(color='none',fill='none') +
  xlab('years since vaccination') +
  ylab('anti-Vi IgG [EU/ml]')
ggsave('scratch/figures/fitted_titers.png',units='in',width=9, height=5)


model_coeffs['T_decay']*1/(1+exp(model_coeffs['beta_T_decay_age']*seq(0,75,by=1)))
model_coeffs['alpha']*(exp(model_coeffs['beta_alpha_age']*((seq(0,75,by=1)^exp(model_coeffs['g'])))))

## boosting

names(d_imm)

boosting_data = d_imm |> filter(study %in% c('Kossaczka1999','Mohan2015')) |>
  filter(antigen=='IgG') |>
  filter( (study == 'Kossaczka1999' & age_min_years == 2) | 
            (study == 'Kossaczka1999' & age_min_years == 18) | 
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
  facet_wrap('vaccine ~ age_label') +
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
                color=vaccine_schedule,group=group_label,linetype=assay_simple)) +
  geom_point(aes(x=day,y=elisa_U_per_ml_median_numeric,
                 shape=study,color=vaccine_schedule,group=group_label)) +
  facet_wrap('vaccine ~ age_mean') +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  xlab('pre-vaccination IgG [EU/ml]') + ylab('post-vaccination IgG [EU/ml]')
ggsave('scratch/figures/pre_post_titers.png',units='in',width=7, height=5)


# convert to fold rise
boosting_data = boosting_data |>
  select(study, vaccine,age_label, vaccine_schedule,antigen, age_min_years,age_max_years,age_mean,age_label, 
         day,elisa_U_per_ml_median_numeric,group_label) 

boosting_data$fold_rise[boosting_data$day>0] = 
  boosting_data$elisa_U_per_ml_median_numeric[boosting_data$day>0]/boosting_data$elisa_U_per_ml_median_numeric[boosting_data$day==0]

boosting_data$pre_vax_elisa[boosting_data$day>0] = boosting_data$elisa_U_per_ml_median_numeric[boosting_data$day==0] 

boosting_data = boosting_data |> filter(day>0)
  
# shoehorn in some natural infection pseudodata, from a mix of https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(22)00114-8/fulltext
# and explore_titer_censoring.R (which are both roughly compatible)

boosting_data = boosting_data |> 
  rbind(setNames(data.frame(matrix(nrow=3,ncol=ncol(boosting_data))),names(boosting_data)) |>
          mutate(study = 'natural',
                 vaccine = 'infection',
                 vaccine_schedule='none',
                 age_label = c('<5','5-15','16+'),
                 pre_vax_elisa = c(0.6,165,202), # units seem comparable but it is a different assay...
                 fold_rise = c(325/0.6,400/165,650/202)) |>
          mutate(group_label = interaction(study,vaccine,age_label, vaccine_schedule)))

# filter out natural infection point starting from limit of detection?
boosting_data = boosting_data |> filter(pre_vax_elisa>1)

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
                  C_max=model_coeffs[1]*exp(model_coeffs[4]*age_means[k])*exp(Delta_coeff_C[n]),
                  T_decay=model_coeffs[2]*exp(model_coeffs[5]*age_means[k]),
                  alpha=model_coeffs[3]*exp(model_coeffs[6]*age_means[k])) /
    titer_vs_time(t=boosting_data$day[row],
                  C_max=model_coeffs[1]*exp(model_coeffs[4]*age_means[k])*exp(Delta_coeff_C[n]),
                  T_decay=model_coeffs[2]*exp(model_coeffs[5]*age_means[k]),
                  alpha=model_coeffs[3]*exp(model_coeffs[6]*age_means[k]))
}
boosting_data$fold_rise_adjusted[is.na(boosting_data$fold_rise_adjusted)]=boosting_data$fold_rise[is.na(boosting_data$fold_rise_adjusted)]

plot(boosting_data$fold_rise,boosting_data$fold_rise_adjusted-boosting_data$fold_rise)

# boosting model equations and fitting
fold_rise_model = function(CoP_pre,mu_0,CoP_max, CoP_min=1){
  fold_rise = 10^(mu_0*(1-(log10(CoP_pre)-log10(CoP_min))/(log10(CoP_max)-log10(CoP_min))))
  return(fold_rise)
}


fold_rise_mlogL = function(mu_TCV,mu_Vips,mu_VirEPA,mu_inf,CoP_max){
  mlogL=0
  for (row in 1:nrow(boosting_data)){
    if (boosting_data$vaccine[row]=='Typbar-TCV'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_TCV,CoP_max=CoP_max)
    } else if (boosting_data$vaccine[row]=='Vi-polysaccharide'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_Vips,CoP_max=CoP_max)
    } else if (boosting_data$vaccine[row]=='Vi-rEPA2'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_VirEPA,CoP_max=CoP_max)
    } else if (boosting_data$vaccine[row]=='infection'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_inf,CoP_max=CoP_max)
    }
    fold_rise = fold_rise
  mlogL = mlogL + (log(boosting_data$fold_rise_adjusted[row])-log(fold_rise))^2
  }
  return(mlogL)
}
mod = mle(fold_rise_mlogL,
          start=list(mu_TCV=3,mu_Vips=2,mu_VirEPA=3,mu_inf=1,CoP_max=10^3.5),
          control=list(parscale=c(0.01,0.01,0.01,0.01,1)))
summary(mod) 
model_coeffs = coef(mod)


vaccines = unique(boosting_data$vaccine)
fitted_boosting = expand.grid(vaccine=vaccines,
                              pre_vax_elisa = 10^seq(0,4,by=0.01))
for (row in 1:nrow(fitted_boosting)){
  if (fitted_boosting$vaccine[row]=='Typbar-TCV'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=model_coeffs[1],CoP_max=model_coeffs[5])
  } else if (fitted_boosting$vaccine[row]=='Vi-polysaccharide'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=model_coeffs[2],CoP_max=model_coeffs[5])
  } else if (fitted_boosting$vaccine[row]=='Vi-rEPA2'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=model_coeffs[3],CoP_max=model_coeffs[5])
  } else if (fitted_boosting$vaccine[row]=='infection'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=model_coeffs[4],CoP_max=model_coeffs[5])
  }
}

# yay! this works, except that there are many observed values above the so-called maximum. 
ggplot() +
  geom_line(data=fitted_boosting,aes(x=pre_vax_elisa,y=fold_rise,
                                     color=vaccine,group=vaccine)) +
  geom_point(data=boosting_data,aes(x=pre_vax_elisa,y=fold_rise_adjusted,
                 shape=age_label,color=vaccine,group=group_label)) +
  theme_bw() +
  scale_y_continuous(trans='log10',limits=c(1,10^3.51),breaks=round(10^seq(0,3.5,by=0.5),0),minor_breaks = NULL) +
  scale_x_continuous(trans='log10',limits=c(1,10^3.51),breaks=round(10^seq(0,3.5,by=0.5),0),minor_breaks = NULL) +
  xlab('pre-vaccination IgG [EU/ml]') + ylab('fold-rise')

ggsave('scratch/figures/boosting_model_foldrise_data.png',units='in',width=5, height=3)

ggplot() +
  geom_line(data=fitted_boosting,aes(x=pre_vax_elisa,y=fold_rise*pre_vax_elisa,
                                     color=vaccine,group=vaccine)) +
  geom_point(data=boosting_data,aes(x=pre_vax_elisa,y=fold_rise_adjusted*pre_vax_elisa,
                                    shape=age_label,color=vaccine,group=group_label)) +
  theme_bw() +
  scale_y_continuous(trans='log10',limits=c(1,10^4.2),breaks=round(10^seq(0,4,by=0.5),0),minor_breaks = NULL) +
  scale_x_continuous(trans='log10',limits=c(1,10^3.51),breaks=round(10^seq(0,3.5,by=0.5),0),minor_breaks = NULL) +
  xlab('pre-vaccination IgG [EU/ml]') + ylab('post-vaccination IgG [EU/ml]')
ggsave('scratch/figures/boosting_model_pre_post.png',units='in',width=5, height=3)


# hmmm wondering if maybe Vi isn't actually a bad correlate of immunity for natural infection, but that
# natural infection isn't all that immunogenic, the human responses are quite heterogeneous, and
# the elisa assays aren't very good....


max(d_imm$elisa_U_per_ml_median_numeric)
max(d_imm$elisa_U_per_ml_upper_numeric)



# actually hmmm the one really high titer is an outlier, which is revealing that linear proportional scaling of the old assay to new is probably wrong
# otherwise, this all works to a factor of 2, which is within expected precision.
max_titer_df =d_imm |> 
  filter(antigen == 'IgG') |>
  filter(day<50 & day >0) |>
  filter(vaccine != 'Vi_polysaccharide') |>
  group_by(age_mean,age_label,study,day,assay,vaccine) |>
  summarize(max_EUpml_median = max(elisa_U_per_ml_median_numeric),
                                          max_EUpml_upper = max(elisa_U_per_ml_upper_numeric))
print(max_titer_df,n=30)

max_titer_df |>
  ggplot() +
    geom_point(aes(x=age_mean,y=max_EUpml_median,color=assay,shape=vaccine)) +
    geom_hline(yintercept =model_coeffs[4],linetype='dashed') +
    scale_y_continuous(trans='log2',limits=2^c(7,14),breaks=2^seq(7,14),minor_breaks = NULL)


# hmm... it is interesting that if, maybe, you happen to get to a higher titer than the boost max, that we wouldn't expect a boost again...
# Ah I see what's happening! It looks like the conjugate vaccines hit the maximum possible antibody level in one dose. 
# And so each re-dose just brings you back up.
# 

# so let's refit dropping Vi-rEPA and infection since I probably don't understand that titer data
boosting_data = boosting_data |> filter(vaccine %in% c('Typbar-TCV','Vi-polysaccharide'))
fold_rise_mlogL = function(mu_TCV,mu_Vips,CoP_max){
  mlogL=0
  for (row in 1:nrow(boosting_data)){
    if (boosting_data$vaccine[row]=='Typbar-TCV'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_TCV,CoP_max=CoP_max)
    } else if (boosting_data$vaccine[row]=='Vi-polysaccharide'){
      fold_rise = fold_rise_model(CoP_pre = boosting_data$pre_vax_elisa[row],mu_0=mu_Vips,CoP_max=CoP_max)
    }
    fold_rise = fold_rise
    mlogL = mlogL + (log(boosting_data$fold_rise_adjusted[row])-log(fold_rise))^2
  }
  return(mlogL)
}

mod = mle(fold_rise_mlogL,
          start=list(mu_TCV=3,mu_Vips=2,CoP_max=2000),
          control=list(parscale=c(0.1,0.1,10)))
summary(mod) 
model_coeffs = coef(mod)

vaccines = unique(boosting_data$vaccine)
fitted_boosting = expand.grid(vaccine=vaccines,
                              pre_vax_elisa = 10^seq(0,4,by=0.01))
for (row in 1:nrow(fitted_boosting)){
  if (fitted_boosting$vaccine[row]=='Typbar-TCV'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=model_coeffs[1],CoP_max=model_coeffs[3])
  } else if (fitted_boosting$vaccine[row]=='Vi-polysaccharide'){
    fitted_boosting$fold_rise[row] = fold_rise_model(CoP_pre = fitted_boosting$pre_vax_elisa[row],mu_0=model_coeffs[2],CoP_max=model_coeffs[3])
  } 
}

ggplot() +
  geom_line(data=fitted_boosting,aes(x=pre_vax_elisa,y=fold_rise,
                                     color=vaccine,group=vaccine)) +
  geom_point(data=boosting_data,aes(x=pre_vax_elisa,y=fold_rise_adjusted,
                                    shape=age_label,color=vaccine,group=group_label)) +
  theme_bw() +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')

max_titer_df |>
  ggplot() +
  geom_point(aes(x=age_mean,y=max_EUpml_median,color=assay,shape=vaccine)) +
  geom_hline(yintercept =model_coeffs[3],linetype='dashed') +
  scale_y_continuous(trans='log2',limits=2^c(7,14),breaks=2^seq(7,14),minor_breaks = NULL)



# It's also interesting that the polysaccharide data looks like it has a shared slope with the conjugates and not a shared intercept.
# My expectation is the intercept should be a human host property, not a vaccine property. But perhaps this is wrong. Not much to say confidently
# either way with just 2 points, but I'm sure there's more Vi-polysaccharide vaccine data out there to add...


