# incidence model scratch

library(tidyverse)
library(patchwork)
library(ggridges)

# age binning function
age_width = function(age_group,age_max){
  round(as.numeric(as.character(fct_recode( age_group,
    `2`='[0,2]',`3`='(2,5]',`5`='(5,10]',`5.1`='(10,15]',!!as.character(age_max-15):=paste('(15,',age_max,']',sep='')))))
}

##### indvidual-level model functions
# titer response function
titer_vs_time = function(t,age_years, CoP_max=10^3.5, CoP_pre=1,
                         T_decay=11, alpha=1.16,
                         beta_T_decay_age=0.057, beta_alpha_age=-0.060,
                         CoP_min=1,T_rise=1.5,t_start=3,t_peak=30.4){
  
  Tda = T_decay * exp(beta_T_decay_age * age_years)
  ka = alpha * exp(beta_alpha_age * age_years)
  
  # power law decay
  titer = CoP_min + (CoP_max-CoP_min)*(1+(t-t_peak)/(ka*Tda))^(-ka)
  
  # simple logistic interpolator rise (this is just for continuity/realism. plays no role in the model)
  titer[t<t_peak ] =
    CoP_pre + (CoP_max-CoP_pre)*
    (1/(1+exp(-(t[t<t_peak] - t_start*5)/T_rise)) - 1/(1+exp(-(0-t_start*5)/T_rise)))/
    (1/(1+exp(-(t_peak - t_start*5)/T_rise)) - 1/(1+exp(-(0-t_start*5)/T_rise)))
  
  return(titer)
}

# fold-rise model. defaults set to natural immunity defaults
fold_rise_model = function(CoP_pre,
                           mu_0=1.25*2,  # playing with hand-adjusting this, which was based on another model, to get better titer distributions
                           CoP_max=10^3.5, CoP_min=1){
  fold_rise = 10^(mu_0*(1-(log10(CoP_pre)-log10(CoP_min))/(log10(CoP_max)-log10(CoP_min))))
  return(fold_rise)
}

# dose response
p_outcome_given_dose = function(dose=1e4,  CoP_pre=1, outcome = 'fever_given_dose', 
                                n50_fever_given_dose=27800, alpha_fever_given_dose=0.84, gamma_fever_given_dose=0.4,
                                n50_infection_given_dose=27800/10,alpha_infection_given_dose = 0.84*2, 
                                gamma_infection_given_dose=0.4/2. # tweaking from previous hand guess to better match comment from Kyra about their estimate of adult immunity
                                ){
  
  if(outcome == 'fever_given_dose'){
    
    p = 1 - (1+dose*(2^(1/alpha_fever_given_dose)-1)/n50_fever_given_dose)^(-alpha_fever_given_dose/(CoP_pre^gamma_fever_given_dose))
    
  } else if (outcome == 'infection_given_dose'){
    
    p = 1 - (1+dose*(2^(1/alpha_infection_given_dose)-1)/n50_infection_given_dose)^(-alpha_infection_given_dose/(CoP_pre^gamma_infection_given_dose))
    
  } else if (outcome == 'fever_given_infection'){
    
    p = (1 - (1+dose*(2^(1/alpha_fever_given_dose)-1)/n50_fever_given_dose)^(-alpha_fever_given_dose/(CoP_pre^gamma_fever_given_dose))) /
      (1 - (1+dose*(2^(1/alpha_infection_given_dose)-1)/n50_infection_given_dose)^(-alpha_infection_given_dose/(CoP_pre^gamma_infection_given_dose)))
    
  }
  
  return(p)
}

# protective efficacy vs CoP_pre
protective_efficacy = function(dose=1e4,  CoP_pre=1, outcome = 'fever_given_dose', CoP_control=1){
  VE = 1 - p_outcome_given_dose(dose,CoP_pre=CoP_pre,outcome = outcome)/p_outcome_given_dose(dose,CoP_pre=CoP_control,outcome = outcome)
  return(VE)
}

expand.grid(dose = 10^seq(0,8,by=0.1),
            CoP_pre = round(10^seq(0,3.5,by=0.5)),
            outcome=factor(c('infection_given_dose','fever_given_dose'),
                           levels=c('infection_given_dose','fever_given_dose','fever_given_infection'))) |>
  group_by(outcome,CoP_pre,dose) |>
  mutate(protective_efficacy = protective_efficacy(dose=dose,CoP_pre=CoP_pre,outcome = outcome)) |>
  mutate(CoP_pre = factor(CoP_pre)) |>
  ggplot() +
  geom_line(aes(x=dose,y=protective_efficacy,group=CoP_pre,color=CoP_pre)) +
  facet_grid('~outcome') +
  theme_bw() +
  ylim(c(0,1)) +
  scale_x_continuous(trans='log10', breaks=10^seq(0,10,by=2),minor_breaks = NULL,labels = scales::trans_format("log10", math_format(10^.x)) ) 
ggsave('scratch/figures/cohort_model_vaccine_efficacy_vs_dose_CoP_pre.png',units='in',width=7, height=3)



### define titer observation model following https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8578 and https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(22)00114-8/fulltext
observe_titer = function(titer,
                         lod = 7, # limit of detection following Quadri2021 for VaccZyme IgG
                         sd_assay=0.3,  # sd from coefficient of variation in log10 space near log10(titer)=O(1) from
                                        # table s3 here https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(22)00114-8/fulltext
                                        # and validated by st control group day 0 here https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0008783
                                        # from looking at many post-vax high titer, I think this is actually an sd w/r/t scaling mean up.
                                        # yes, it should be additive in log space per cited eq 13 https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8578
                         sd_bio = 1/2, # supposedly additive in titer in absolute titer per eq 10 of https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8578
                                        # but I really doubt that it is additive as anitbodies go up
                                        # divide by two because reported value ~1 for IgG is 95% upper limit
                         bio_noise_model = 'additive', 
                         CoP_max=10^3.5){ # needed for bounded_multiplicative model

  if (bio_noise_model == 'additive'){
    measurement = exp(rnorm(length(titer),mean=log(titer),sd=sd_assay))
    measurement = pmax(lod, rnorm(length(titer),mean=measurement,sd=sd_bio))
  } else if (bio_noise_model == 'multiplicative'){ # assume sd_bio should really be a cv and thus log-additive, around lod~1 for their assay
    measurement = pmax(lod, exp(rnorm(length(titer),mean=log(titer),sd=sd_assay+sd_bio)))
  } else if (bio_noise_model == 'bounded_multiplicative'){ # noise decays toward max titer
    sd_bio = sd_bio * (1-log(titer)/log(CoP_max)) # assume proportion off-target antibodies declines as specific immunity goes up
    measurement = pmax(lod, exp(rnorm(length(titer),mean=log(titer),sd=sd_assay+sd_bio)))
  }
  
  is_positive = measurement > lod
  return(data.frame(titer_observed=measurement,is_positive=is_positive))
}

expand.grid(titer=10^seq(0,3.5,by=0.5),
            replicate=1:100) |>
  mutate(observe_titer(titer=titer,bio_noise_model = 'additive')) |>
  ggplot() +
  geom_density_ridges(aes(x=log10(titer_observed),y=factor(log10(titer)),color=is_positive),jittered_points=TRUE)

expand.grid(titer=10^seq(0,3.5,by=0.5),
            replicate=1:100) |>
  mutate(observe_titer(titer=titer,bio_noise_model = 'multiplicative')) |>
  ggplot() +
  geom_density_ridges(aes(x=log10(titer_observed),y=factor(log10(titer)),color=is_positive),jittered_points=TRUE)

expand.grid(titer=10^seq(0,3.5,by=0.5),
            replicate=1:100) |>
  mutate(observe_titer(titer=titer,bio_noise_model = 'bounded_multiplicative')) |>
  ggplot() +
  geom_density_ridges(aes(x=log10(titer_observed),y=factor(log10(titer)),color=is_positive),jittered_points=TRUE)


##### wrap the cohort model in a function
cohort_model = function(exposure_dose,exposure_rate,
                        N=1000,age_max=75,
                        titer_dt=365/12, # monthly timesteps so I don't have to worry about blocking reinfection during current infection
                        max_titer_id = 1000 # saving every titer above ~1000 gets really slow
                        ){
  
  simulation_months = seq(0,age_max*12-1,by=1)
  
  # titer tracker
  titer_df = expand.grid(id=1:max_titer_id, 
                         month=simulation_months,
                         titer=1) |>
    mutate(age_years=month/12,
           age_group=cut(age_years,breaks=c(0,2,5,10,15,age_max),include.lowest = TRUE)) |>
    mutate(age_width = age_width(age_group,age_max))
  
  # disease events tracker
  events_list =  replicate(N, list(exposure_month = NULL, exposure_count=NULL,
                                   infected = NULL, 
                                   fever = NULL), 
                           simplify = FALSE)
  
  # age-dependent exposure rate, to account for kids under 2y having less exposure to food and sewage
  exposure_rate_multiplier = c(rep(0.1,13),rep(0.5,12),rep(1,length(simulation_months)-25))
  
  # expose
  for (id in (1:N)){
    tmp_exposed = rpois(length(exposure_rate_multiplier),exposure_rate*exposure_rate_multiplier)
    events_list[[id]]$exposure_month = which(tmp_exposed>0)
    events_list[[id]]$exposure_count = tmp_exposed[events_list[[id]]$exposure_month]
    
    # declare needed fields
    events_list[[id]]$infected = rep(0, length(events_list[[id]]$exposure_month))
    events_list[[id]]$fever = rep(0, length(events_list[[id]]$exposure_month))
  }
  
  # infect and titer
  for (id in (1:N)){
    
    tmp_titer =rep(1,length(simulation_months))
    
    if (length(events_list[[id]]$exposure_month) >0){
      for (k in 1:length(events_list[[id]]$exposure_month)){
        
        exposure_month = events_list[[id]]$exposure_month[k]
        
        p_once = p_outcome_given_dose(dose=exposure_dose,CoP_pre=tmp_titer[exposure_month],outcome='infection_given_dose')
        p_inf = 1-(1-p_once)^events_list[[id]]$exposure_count[k]
        
        if(runif(1)>p_inf){
          events_list[[id]]$infected[k] = 0
          events_list[[id]]$fever[k]    = 0
        } else {
          
          # update infection list
          events_list[[id]]$infected[k]=1
          
          # update fever list 
          p_fever = p_outcome_given_dose(dose=exposure_dose,CoP_pre=tmp_titer[exposure_month],outcome='fever_given_infection')
          if (runif(n=1)<=p_fever){
            events_list[[id]]$fever[k]=1
          } else {
            events_list[[id]]$fever[k]=0
          }
          
          tmp_titer[exposure_month:length(tmp_titer)]= titer_vs_time(t=(titer_dt*(simulation_months-exposure_month+1))[(exposure_month):length(simulation_months)],
                                                               age_years = exposure_month/12,
                                                               CoP_pre = tmp_titer[exposure_month],
                                                               CoP_max=tmp_titer[exposure_month]*fold_rise_model(CoP_pre = tmp_titer[exposure_month]))
          tmp_titer[exposure_month:length(tmp_titer)] = round(tmp_titer[exposure_month:length(tmp_titer)],2)
        }
        
      }
      
      if(id <= max_titer_id){
        idx = titer_df$id == id
        titer_df$titer[idx] = tmp_titer
      }
    }
  }
  
  # make the events_list into a nice data frame
  not_empty_idx = which(!sapply(events_list, function(x){ is_empty(x$exposure_month)}))
  events_df = tibble(id = not_empty_idx, data = events_list[not_empty_idx]) |>
    unnest_wider(data) |>
    unnest_longer(c(exposure_month,exposure_count,infected,fever), keep_empty = TRUE) |>
    mutate(exposure_age=exposure_month/12,
           age_group=cut(exposure_age,breaks=c(0,2,5,10,15,age_max),include.lowest = TRUE)) |>
    mutate(age_width=age_width(age_group,age_max)) |>
    group_by(id) |>
    mutate(infection_age = if_else(infected==1,exposure_age,NA)) |>
    mutate(fever_age = if_else(fever==1,exposure_age,NA)) |>
    mutate(infected = factor(infected),
           fever = factor(fever)) |>
    mutate(outcome = interaction(infected,fever)) |>
    mutate(outcome = fct_recode(outcome,exposed='0.0',infected='1.0',fever='1.1')) |>
    mutate(outcome = factor(outcome,levels=c('exposed', 'infected', 'fever')))

  # combine some events with titer for plotting
  titer_df = titer_df |>
    left_join(events_df |> select(id,infection_age,infected,fever) |> 
                # mutate(fever = if_else(is.na(fever),0,fever)) |>
                drop_na(infected),by=join_by(id ==id, age_years == infection_age))

  return(list(titer_df=titer_df,
              events_df=events_df,
              config=data.frame(exposure_dose,exposure_rate,N,age_max,titer_dt,max_titer_id)))
  
}

###### run the model for the three different reference incidence levels

# define setting ecology: exposure rate and dose

output=list()

# N_cohort=2e4 # a lot faster for playing
N_cohort=1e6 # made huge to get good stats at lower incidence

# medium
output[['medium']] = cohort_model(exposure_dose = 5e2,
                                  exposure_rate=1/(12*40), # per month
                                  N=N_cohort)

# high
output[['high']] = cohort_model(exposure_dose = 5e3,
                                exposure_rate=1/(12*20), # per month
                                N=N_cohort/5)

# very high
output[['very_high']] = cohort_model(exposure_dose = 5e4,
                                     exposure_rate=1/(12*20), # per month
                                     N=N_cohort/10)

# plots 

# incidence targets
incidence_fever_targets = c(medium = 53,high=214,very_high=1255)

# calculate fever and infection incidence
for (k in 1:length(output)){
  
  N = output[[k]]$config$N
  
  output[[k]]$incidence_vs_age =
    output[[k]]$events_df |> group_by(age_group) |>
    summarize(incidence_fever = sum(fever==1)/(N*unique(age_width))*1e5,
              incidence_infection = sum(infected==1)/(N*unique(age_width))*1e5,
              symptomatic_fraction = sum(fever==1,na.rm=TRUE)/sum(infected==1,na.rm=TRUE),
              age_width = unique(age_width)) |>
    mutate(incidence_fever_overall = sum(incidence_fever*age_width/sum(age_width)),
           incidence_infection_overall = sum(incidence_infection*age_width/sum(age_width)),
           incidence_fever_target = incidence_fever_targets[names(output)[k]])
}

# calculate seroincidence
for (k in 1:length(output)){
  
  N = output[[k]]$config$max_titer_id
  
  output[[k]]$titer_df = output[[k]]$titer_df[,1:8] |>
    group_by(id) |>
    mutate(fold_rise = c(0,diff(titer))) |>
    mutate(seroconversion = fold_rise>0) |>
    mutate(observe_titer(titer,bio_noise_model = 'additive')) |>
    rename(titer_observed_additive = titer_observed, is_positive_additive=is_positive) |>
    mutate(fold_rise_additive = c(0,diff(titer_observed_additive))) |>
    mutate(seroconversion_observed_additive = fold_rise_additive>3*sqrt(2)) |>
    mutate(observe_titer(titer,bio_noise_model = 'bounded_multiplicative')) |>
    mutate(fold_rise_bounded_multiplicative = c(0,diff(titer_observed))) |>
    mutate(seroconversion_observed_bounded_multiplicative = fold_rise_bounded_multiplicative>3*sqrt(2)) 
# 
#   ggplot(output[[k]]$titer_df |> filter(id<=20)) +
#     geom_point(aes(x=age_years,y=fold_rise,color=seroconversion)) +
#     facet_wrap('id')
#   
#   ggplot(output[[k]]$titer_df |> filter(id<=20)) +
#     geom_point(aes(x=age_years,y=fold_rise_additive,color=seroconversion_observed_additive)) +
#     facet_wrap('id')
#   
#   ggplot(output[[k]]$titer_df |> filter(id<=20)) +
#     geom_point(aes(x=age_years,y=fold_rise_bounded_multiplicative,color=seroconversion_observed_bounded_multiplicative)) +
#     facet_wrap('id')
  
  output[[k]]$incidence_vs_age= output[[k]]$incidence_vs_age |>
    left_join(
      output[[k]]$titer_df |>
        group_by(age_group) |>
        summarize(incidence_sero = sum(seroconversion==TRUE)/(N*unique(age_width))*1e5,
                  incidence_sero_additive = sum(seroconversion_observed_additive==TRUE)/(N*unique(age_width))*1e5,
                  incidence_sero_bounded_multiplicative = sum(seroconversion_observed_bounded_multiplicative==TRUE)/(N*unique(age_width))*1e5)) |>
    mutate(ratio_additive_sero = incidence_sero_additive/incidence_sero,
           ratio_bounded_multiplicative_sero = incidence_sero_bounded_multiplicative/incidence_sero,
           ratio_additive_fever = incidence_sero_additive/incidence_fever,
           ratio_bounded_multiplicative_fever = incidence_sero_bounded_multiplicative/incidence_fever) |>
    mutate(symptomatic_fraction_additive_sero = symptomatic_fraction / ratio_additive_sero,
           symptomatic_fraction_bounded_multiplicative_sero = symptomatic_fraction / ratio_bounded_multiplicative_sero) |>
    mutate(incidence_sero_overall = sum(incidence_sero*age_width/sum(age_width)),
           incidence_sero_additive_overall = sum(incidence_sero_additive*age_width/sum(age_width)),
           incidence_sero_bounded_multiplicative_overall = sum(incidence_sero_bounded_multiplicative*age_width/sum(age_width)),
           symptomatic_fraction_overall = sum(symptomatic_fraction*age_width/sum(age_width)),
           symptomatic_fraction_overall_additive_sero = sum(symptomatic_fraction_additive_sero*age_width/sum(age_width)),
           symptomatic_fraction_overall_bounded_multiplicative_sero = sum(symptomatic_fraction_bounded_multiplicative_sero*age_width/sum(age_width)))
     
}

# I THINK the additive observation model is equivalent to the sero-epi incidence assumptions 
# used by these people https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(22)00114-8/fulltext#fig3
# if seroincidence is then based on the observed titers without accounting for the noise, 
# it incidence will be over-estimated. https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8578 Factor of 10 not impossible.
# But I can't believe this is really equivalent to what anyone is doing. The daily differencing assumed
# above has to be not equivalent to how anyone does this, I hope...

##### make lots of plots

p_incidence_fever=list()
p_incidence_infection=list()
p_symptomatic_fraction=list()
p_incidence_sero_additive=list()
p_symptomatic_fraction_additive_sero=list()
for (k in 1:length(output)){
  p_incidence_fever[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_fever),stat='identity') +
    geom_hline(aes(yintercept=incidence_fever_overall[1]),linetype='solid') +
    geom_hline(aes(yintercept=incidence_fever_target[1]),linetype='dashed') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of fever per 100k') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
  
  p_incidence_infection[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_infection),stat='identity') +
    geom_hline(aes(yintercept=incidence_infection_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of infection per 100k') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
  
  p_symptomatic_fraction[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=symptomatic_fraction),stat='identity') +
    geom_hline(aes(yintercept=symptomatic_fraction_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('symptomatic fraction') +
    ylim(c(0,1)) +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
  
  p_incidence_sero_additive[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_sero_additive),stat='identity') +
    geom_hline(aes(yintercept=incidence_sero_additive_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of infection per 100k (additive sero model)') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
  
  p_symptomatic_fraction_additive_sero[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=symptomatic_fraction_additive_sero),stat='identity') +
    geom_hline(aes(yintercept=symptomatic_fraction_overall_additive_sero[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('symptomatic fraction (additive seroepi model)') +
    ylim(c(0,1)) +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
}

wrap_plots(p_incidence_fever) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_incidence_fever_by_age.png',units='in',width=7,height=4)

wrap_plots(c(p_incidence_infection,p_incidence_sero_additive),ncol=3) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_incidence_infection_by_age.png',units='in',width=7,height=4)

wrap_plots(c(p_symptomatic_fraction,p_symptomatic_fraction_additive_sero),ncol=3) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_symptomatic_fraction.png',units='in',width=7,height=4)


wrap_plots(c(p_symptomatic_fraction[3],p_symptomatic_fraction_additive_sero[3]),ncol=2)  + 
  plot_layout(guides = "collect",axes='collect') & 
  scale_y_continuous(minor_breaks = seq(0,1,by=0.05),limits=c(0,0.65)) 
ggsave('scratch/figures/cohort_model_symptomatic_fraction.png',units='in',width=7,height=4)


p_exposure=list()
for (k in 1:length(output)){
  p_exposure[[k]] = ggplot(output[[k]]$events_df |> filter(id<=40)) +
                      geom_point(aes(x=exposure_age,color=outcome,y=id)) +
                      theme_bw() +
                      xlab('age') +
      labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
              subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                    1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) +
    scale_color_discrete(drop = FALSE)
}
wrap_plots(p_exposure) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_individual_level_exposure_examples.png',units='in',width=7,height=6)

p_titer_examples=list()
for (k in 1:length(output)){
  p_titer_examples[[k]] = ggplot() +
    geom_line(data=output[[k]]$titer_df |> filter(id<=20),aes(x=age_years,y=titer)) + 
    geom_point(data=output[[k]]$titer_df |> filter(id<=20) |> filter(!is.na(infected)),
               aes(x=age_years,y=titer,color=fever)) +
    scale_y_continuous(trans='log10',limits=c(1,10^3.5)) +
    facet_wrap('id') +
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))  +
    scale_color_discrete(drop = FALSE)
}
wrap_plots(p_titer_examples) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_individual_level_titer_examples.png',units='in',width=12,height=6)


p_titer_summary=list()
p_titer_density=list()
p_protective_efficacy_summary=list()
p_titer_observed_density_additive_noise=list()
p_titer_observed_density_bounded_multiplicative_noise=list()
for (k in 1:length(output)){
  
  exposure_dose = output[[k]]$config$exposure_dose
  
  tmp_titer_summary =output[[k]]$titer_df |>
    mutate(protective_efficacy = protective_efficacy(dose = exposure_dose,CoP_pre = titer, outcome='infection_given_dose')) |>
    group_by(age_years) |>
    summarize(titer_gmt = exp(mean(log(titer))),
           titer_median = median(titer),
           titer_upper_iqr = quantile(titer,probs=0.75),
           titer_lower_iqr = quantile(titer,probs=0.25),
           titer_upper_95 = quantile(titer,probs=0.975),
           titer_lower_95 = quantile(titer,probs=0.025),
           protective_efficacy_median = median(protective_efficacy),
           protective_efficacy_upper_iqr = quantile(protective_efficacy,probs=0.75),
           protective_efficacy_lower_iqr = quantile(protective_efficacy,probs=0.25),
           protective_efficacy_upper_95 = quantile(protective_efficacy,probs=0.975),
           protective_efficacy_lower_95 = quantile(protective_efficacy,probs=0.025))
  
  p_titer_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_ribbon(aes(x=age_years,ymin=titer_lower_95,ymax=titer_upper_95),alpha=0.2)+
    geom_ribbon(aes(x=age_years,ymin=titer_lower_iqr,ymax=titer_upper_iqr),alpha=0.2)+
    geom_line(aes(x=age_years,y=titer_median)) + 
    geom_line(aes(x=age_years,y=titer_gmt),linetype='dashed') + 
    scale_y_continuous(limits=c(1,10^3.5)) +
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) +
    ylab('anti-Vi IgG EU/ml')
  
  p_protective_efficacy_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_ribbon(aes(x=age_years,ymin=protective_efficacy_lower_95,ymax=protective_efficacy_upper_95),alpha=0.2)+
    geom_ribbon(aes(x=age_years,ymin=protective_efficacy_lower_iqr,ymax=protective_efficacy_upper_iqr),alpha=0.2)+
    geom_line(aes(x=age_years,y=protective_efficacy_median)) + 
    scale_y_continuous(limits=c(0,1)) +
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) +
    ylab('protective efficacy')
  
  
  sampled_titer_df = output[[k]]$titer_df |>
    filter((age_years>0.75 & age_years<=2 & id<=200) | # sample sizes like quadri2021-like sample sizes
             (age_years >2 & age_years<=5 & id<600) | 
             (age_years >5 & age_years<=10 & id<300) | 
             (age_years >10 & age_years<=15 & id<300) | 
             (age_years >15 & id<300)) |>
    group_by(age_group,id) |>
    slice_sample(n=1) |> 
    ungroup() |>
    mutate(observe_titer(titer,bio_noise_model = 'additive')) |>
    rename(titer_observed_additive = titer_observed, is_positive_additive=is_positive) |>
    mutate(observe_titer(titer,bio_noise_model = 'bounded_multiplicative')) 
  
  p_titer_density[[k]] = ggplot(sampled_titer_df |> filter(titer>7)) +
    geom_density_ridges(aes(x=titer,y=age_group),scale=0.9,jittered_points = TRUE,bandwidth = 0.1) + 
    scale_x_continuous(trans='log10',limits=c(7,10^4.5)) +
    theme_bw() + 
    xlab('titer, given above limit of detection') +
    ylab('') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) 
  
  p_titer_observed_density_additive_noise[[k]] = ggplot(sampled_titer_df |> filter(is_positive_additive)) +
    geom_density_ridges(aes(x=titer_observed_additive,y=age_group),scale=0.9,jittered_points = TRUE,bandwidth = 0.1) + 
    scale_x_continuous(trans='log10',limits=c(7,10^4.5)) +
    theme_bw() + 
    xlab('observed titer (additive model)') +
    ylab('') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) 
  
  p_titer_observed_density_bounded_multiplicative_noise[[k]] = ggplot(sampled_titer_df |> filter(is_positive)) +
    geom_density_ridges(aes(x=titer_observed,y=age_group),scale=0.9,jittered_points = TRUE,bandwidth = 0.1) + 
    scale_x_continuous(trans='log10',limits=c(7,10^4.5)) +
    theme_bw() + 
    xlab('observed titer (bounded multiplicative model)') +
    ylab('') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) 
}
wrap_plots(p_titer_summary) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_titer_vs_age_summary.png',units='in',width=7,height=4)

wrap_plots(p_protective_efficacy_summary) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_protective_efficacy_infection_vs_age_summary.png',units='in',width=7,height=4)

wrap_plots(c(p_titer_density,
             p_titer_observed_density_bounded_multiplicative_noise,
             p_titer_observed_density_additive_noise),
           ncol=3) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_titer_density.png',units='in',width=7,height=4)
# I have strong prior that the bounded_multiplicative biology noise model is the right one, 
# because the fraction of the titer that is truly pathogen specific should increase with higher true pathogen exposure




################
# for comparison to https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(22)00114-8/fulltext#fig3
# last age bin only goes to age 25

output_aiemjoy =output  

# calculate fever and infection incidence
for (k in 1:length(output_aiemjoy)){
  
  N = output_aiemjoy[[k]]$config$N
  
  output_aiemjoy[[k]]$events_df = output_aiemjoy[[k]]$events_df |>
    filter(exposure_age<=25) |>
    ungroup() |>
    mutate(age_group=cut(exposure_age,breaks=c(0,2,5,10,15,25),include.lowest = TRUE)) |>
    mutate(age_width = age_width(age_group,25))
  
  output_aiemjoy[[k]]$incidence_vs_age =
    output_aiemjoy[[k]]$events_df |> group_by(age_group) |>
    summarize(incidence_fever = sum(fever==1)/(N*unique(age_width))*1e5,
              incidence_infection = sum(infected==1)/(N*unique(age_width))*1e5,
              symptomatic_fraction = sum(fever==1,na.rm=TRUE)/sum(infected==1,na.rm=TRUE),
              age_width = unique(age_width)) |>
    mutate(incidence_fever_overall = sum(incidence_fever*age_width/sum(age_width)),
           incidence_infection_overall = sum(incidence_infection*age_width/sum(age_width)),
           incidence_fever_target = incidence_fever_targets[names(output_aiemjoy)[k]])
}

# calculate seroincidence
for (k in 1:length(output_aiemjoy)){
  
  N = output_aiemjoy[[k]]$config$max_titer_id
  
  output_aiemjoy[[k]]$titer_df = output_aiemjoy[[k]]$titer_df[,1:8] |>
    filter(age_years<=25) |>
    ungroup() |>
    mutate(age_group=cut(age_years,breaks=c(0,2,5,10,15,25),include.lowest = TRUE)) |>
    mutate(age_width = age_width(age_group,25)) |>
    group_by(id) |>
    mutate(fold_rise = c(0,diff(titer))) |>
    mutate(seroconversion = fold_rise>0) |>
    mutate(observe_titer(titer,bio_noise_model = 'additive')) |>
    rename(titer_observed_additive = titer_observed, is_positive_additive=is_positive) |>
    mutate(fold_rise_additive = c(0,diff(titer_observed_additive))) |>
    mutate(seroconversion_observed_additive = fold_rise_additive>3*sqrt(2)) |>
    mutate(observe_titer(titer,bio_noise_model = 'bounded_multiplicative')) |>
    mutate(fold_rise_bounded_multiplicative = c(0,diff(titer_observed))) |>
    mutate(seroconversion_observed_bounded_multiplicative = fold_rise_bounded_multiplicative>3*sqrt(2)) 

  output_aiemjoy[[k]]$incidence_vs_age= output_aiemjoy[[k]]$incidence_vs_age |>
    left_join(
      output_aiemjoy[[k]]$titer_df |>
        group_by(age_group) |>
        summarize(incidence_sero = sum(seroconversion==TRUE)/(N*unique(age_width))*1e5,
                  incidence_sero_additive = sum(seroconversion_observed_additive==TRUE)/(N*unique(age_width))*1e5,
                  incidence_sero_bounded_multiplicative = sum(seroconversion_observed_bounded_multiplicative==TRUE)/(N*unique(age_width))*1e5)) |>
    mutate(ratio_additive_sero = incidence_sero_additive/incidence_sero,
           ratio_bounded_multiplicative_sero = incidence_sero_bounded_multiplicative/incidence_sero,
           ratio_additive_fever = incidence_sero_additive/incidence_fever,
           ratio_bounded_multiplicative_fever = incidence_sero_bounded_multiplicative/incidence_fever) |>
    mutate(symptomatic_fraction_additive_sero = symptomatic_fraction / ratio_additive_sero,
           symptomatic_fraction_bounded_multiplicative_sero = symptomatic_fraction / ratio_bounded_multiplicative_sero) |>
    mutate(incidence_sero_overall = sum(incidence_sero*age_width/sum(age_width)),
           incidence_sero_additive_overall = sum(incidence_sero_additive*age_width/sum(age_width)),
           incidence_sero_bounded_multiplicative_overall = sum(incidence_sero_bounded_multiplicative*age_width/sum(age_width)),
           symptomatic_fraction_overall = sum(symptomatic_fraction*age_width/sum(age_width)),
           symptomatic_fraction_overall_additive_sero = sum(symptomatic_fraction_additive_sero*age_width/sum(age_width)),
           symptomatic_fraction_overall_bounded_multiplicative_sero = sum(symptomatic_fraction_bounded_multiplicative_sero*age_width/sum(age_width)))
  
}

# yeah, this can't be right the seroincidence by age grows exponentially the way I'm doing it, 
# but published estimates don't do that.
# so it's presumably a coincinence I get very close to the same all-ages average seroincidence...
##### make lots of plots

p_incidence_fever=list()
p_incidence_infection=list()
p_symptomatic_fraction=list()
p_incidence_sero_additive=list()
p_symptomatic_fraction_additive_sero=list()
for (k in 1:length(output_aiemjoy)){
  p_incidence_fever[[k]]=ggplot(output_aiemjoy[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_fever),stat='identity') +
    geom_hline(aes(yintercept=incidence_fever_overall[1]),linetype='solid') +
    geom_hline(aes(yintercept=incidence_fever_target[1]),linetype='dashed') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of fever per 100k') +
    labs(title = paste(sub('_',' ',names(output_aiemjoy)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output_aiemjoy[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output_aiemjoy[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
  
  p_incidence_infection[[k]]=ggplot(output_aiemjoy[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_infection),stat='identity') +
    geom_hline(aes(yintercept=incidence_infection_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of infection per 100k') +
    labs(title = paste(sub('_',' ',names(output_aiemjoy)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output_aiemjoy[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output_aiemjoy[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
  
  p_symptomatic_fraction[[k]]=ggplot(output_aiemjoy[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=symptomatic_fraction),stat='identity') +
    geom_hline(aes(yintercept=symptomatic_fraction_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('symptomatic fraction') +
    ylim(c(0,1)) +
    labs(title = paste(sub('_',' ',names(output_aiemjoy)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output_aiemjoy[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output_aiemjoy[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
  
  p_incidence_sero_additive[[k]]=ggplot(output_aiemjoy[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_sero_additive),stat='identity') +
    geom_hline(aes(yintercept=incidence_sero_additive_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of infection per 100k (additive sero model)') +
    labs(title = paste(sub('_',' ',names(output_aiemjoy)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output_aiemjoy[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output_aiemjoy[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
  
  p_symptomatic_fraction_additive_sero[[k]]=ggplot(output_aiemjoy[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=symptomatic_fraction_additive_sero),stat='identity') +
    geom_hline(aes(yintercept=symptomatic_fraction_overall_additive_sero[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('symptomatic fraction (additive seroepi model)') +
    ylim(c(0,1)) +
    labs(title = paste(sub('_',' ',names(output_aiemjoy)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output_aiemjoy[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output_aiemjoy[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
}

wrap_plots(p_incidence_fever) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/aiemjoy/cohort_model_incidence_fever_by_age.png',units='in',width=7,height=4)

wrap_plots(c(p_incidence_infection,p_incidence_sero_additive),ncol=3) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/aiemjoy/cohort_model_incidence_infection_by_age.png',units='in',width=7,height=4)

wrap_plots(c(p_symptomatic_fraction,p_symptomatic_fraction_additive_sero),ncol=3) + plot_layout(guides = "collect",axes='collect') & 
  scale_y_continuous(minor_breaks = seq(0,1,by=0.05),limits=c(0,0.65)) 
ggsave('scratch/figures/aiemjoy/cohort_model_symptomatic_fraction.png',units='in',width=7,height=4)


