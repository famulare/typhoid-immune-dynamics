# cohort_incidence_model_ve_study_scratch.R
#
# EXPLORATORY BRANCH - VACCINE EFFICACY STUDY SIMULATION
#
# This file extends the cohort incidence model with:
# - TCV vaccination simulation (vaccination_rate parameter, random vaccination age)
# - Separate fold-rise parameters for infection vs TCV (mu_0_inf, mu_0_tcv)
# - VE study analysis: matched treatment/control arms, VE by time interval
# - Different titer waning parameters than canonical version
#
# Key additions beyond canonical cohort model:
# - fold_rise_model() has source='infection' or 'tcv' options
# - cohort_model() has vaccination_rate parameter
# - VE calculation over 1-year intervals out to 20+ years
#
# For the canonical deployable model, see main branch.
#
# =============================================================================

library(tidyverse)
library(patchwork)
library(ggridges)
library(scales)

# age width function
age_width = function(age_group,age_max){
  round(as.numeric(as.character(fct_recode( age_group,
                                            `2`='[0,2]',`3`='(2,5]',`5`='(5,10]',`5.1`='(10,15]',!!as.character(age_max-15):=paste('(15,',age_max,']',sep='')))))
}

# =============================================================================
# INTRAHOST IMMUNITY MODEL FUNCTIONS
# =============================================================================

# Titer response function
# Different parameters than canonical: T_decay=14.7, alpha=1.18, different age formulas
titer_vs_time = function(t,age_years, CoP_peak=1000, CoP_pre=1,
                         T_decay=14.7, alpha=1.18,
                         beta_T_decay_age=-0.61, beta_alpha_age=-0.091,lambda_alpha_age=0.9,
                         CoP_min=1,T_rise=1.5,t_start=3,t_peak=30.4){

  Tda = T_decay * 1/(1+exp(beta_T_decay_age * age_years))
  ka = alpha * exp(beta_alpha_age * (age_years^lambda_alpha_age))

  # power law decay
  titer = CoP_min + (CoP_peak-CoP_min)*(1+(t-t_peak)/(ka*Tda))^(-ka)

  # simple logistic interpolator rise
  titer[t<t_peak ] =
    CoP_pre + (CoP_peak-CoP_pre)*
    (1/(1+exp(-(t[t<t_peak] - t_start*5)/T_rise)) - 1/(1+exp(-(0-t_start*5)/T_rise)))/
    (1/(1+exp(-(t_peak - t_start*5)/T_rise)) - 1/(1+exp(-(0-t_start*5)/T_rise)))

  return(titer)
}

# Demo plot
expand.grid(t=seq(0,10,by=1/24), age_years = factor(c(1,5,15,45))) |>
  mutate(titer = titer_vs_time(t=t*365,age=as.numeric(age_years))) |>
  ggplot() +
  geom_line(aes(x=t,y=titer,color=age_years)) +
  theme_bw() + scale_y_continuous(trans='log10') + xlab('years post response') + ylab('')

# Fold-rise model with separate parameters for infection vs TCV
# KEY DIFFERENCE: source parameter allows 'infection' or 'tcv'
fold_rise_model = function(CoP_pre,
                           mu_0_inf=1.5,
                           mu_0_tcv=3.16,
                           CoP_max=10^3.3, CoP_min=1,
                           sigma_0=(2.9/7.2)*1.5,
                           response='individual',
                           source='infection'){
  if (source=='infection'){
    mu = mu_0_inf
  } else if (source == 'tcv'){
    mu = mu_0_tcv
  }

  if (response=='individual'){
    mu = pmax(0,rnorm(length(CoP_pre),mean=mu,sd=sigma_0))
  }
  fold_rise = 10^(mu*(1-(log10(CoP_pre)-log10(CoP_min))/(log10(CoP_max)-log10(CoP_min))))
  return(fold_rise)
}

# Demo plot comparing infection vs TCV fold-rise
pl_df=expand.grid(CoP_pre=10^seq(0,3.3,by=0.1),
                  source = c('infection','tcv')) |>
  group_by(CoP_pre,source) |>
  mutate(fold_rise = fold_rise_model(CoP_pre=CoP_pre, response='median',source=source)) |>
  mutate(CoP_post = fold_rise * CoP_pre)
(ggplot(pl_df) +
  geom_line(aes(x=CoP_pre,y=fold_rise,color=source)) +
  theme_bw() + scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') +
  xlab('pre-challenge titer') + ylab('fold-rise')) +
  (ggplot(pl_df) +
     geom_line(aes(x=CoP_pre,y=CoP_post,color=source)) +
     theme_bw() + scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') +
     xlab('pre-challenge titer') + ylab('post-challenge titer'))

# Example: lifetime titer trajectory with infections
data.frame(t=seq(0,15,by=1/12), titer=1) |>
  mutate(titer = if_else(t<=2, titer,
                         titer_vs_time(t=(t-2)*365,age=2,
                                       CoP_pre=titer[t==2],
                                       CoP_peak = titer[t==2]*fold_rise_model(CoP_pre = titer[t==2])))) |>
  mutate(titer = if_else(t<=7, titer,
                         titer_vs_time(t=(t-7)*365,age=7,
                                       CoP_pre=titer[t==7],
                                       CoP_peak = titer[t==7]*fold_rise_model(CoP_pre = titer[t==7])))) |>
  ggplot() +
  geom_line(aes(x=t,y=titer)) +
  theme_bw() + scale_y_continuous(trans='log10') + xlab('age [years]') + ylab('Anti-Vi IgG [EU/ml]')


# =============================================================================
# DOSE-RESPONSE MODEL
# =============================================================================

p_outcome_given_dose = function(dose=1e4,  CoP_pre=1, outcome = 'fever_given_dose',
                                n50_fever_given_dose=27800, alpha_fever_given_dose=0.84,
                                gamma_fever_given_dose=0.4,
                                n50_infection_given_dose=27800/40,alpha_infection_given_dose = 0.84*2,
                                gamma_infection_given_dose=0.4/2
                                ){

  if(outcome == 'fever_given_dose'){
    p = 1 - (1+dose*(2^(1/alpha_fever_given_dose)-1)/n50_fever_given_dose)^(-alpha_fever_given_dose/(CoP_pre^gamma_fever_given_dose))
  } else if (outcome == 'infection_given_dose'){
    p = 1 - (1+dose*(2^(1/alpha_infection_given_dose)-1)/n50_infection_given_dose)^
      (-alpha_infection_given_dose/(CoP_pre^gamma_infection_given_dose))
  } else if (outcome == 'fever_given_infection'){
    p = (1 - (1+dose*(2^(1/alpha_fever_given_dose)-1)/n50_fever_given_dose)^
           (-alpha_fever_given_dose/(CoP_pre^gamma_fever_given_dose))) /
      (1 - (1+dose*(2^(1/alpha_infection_given_dose)-1)/n50_infection_given_dose)^
         (-alpha_infection_given_dose/(CoP_pre^gamma_infection_given_dose)))
  }
  return(p)
}

# Demo plot
plot_dat = expand.grid(dose = 10^seq(0,9,by=0.1),
            CoP_pre = c(50,round(10^seq(0,3.3,by=0.5))),
            outcome=factor(c('infection given dose','fever given dose','fever given infection'),
                           levels=c('infection given dose','fever given dose','fever given infection'))) |>
  group_by(outcome,CoP_pre,dose) |>
  mutate(probability = p_outcome_given_dose(dose=dose,CoP_pre=CoP_pre,outcome = gsub(' ','_',outcome))) |>
  mutate(CoP_pre = factor(CoP_pre))

ggplot() +
  geom_line(data = plot_dat |> filter(CoP_pre!=50),
            aes(x=dose,y=probability,group=CoP_pre,color=CoP_pre)) +
  geom_line(data = plot_dat |> filter(CoP_pre==50 & outcome == 'fever given dose'),
            aes(x=dose,y=probability,group=CoP_pre),color='black',linewidth=1) +
  facet_grid('~outcome') +
  theme_bw() +
  ylim(c(0,1)) +
  scale_x_continuous(trans='log10', breaks=10^seq(0,10,by=2),minor_breaks = NULL,
                     labels = trans_format("log10", math_format(10^.x)) ) +
  xlab('bacilli') +
  labs(color = "anti-Vi titer\n[EU/ml]")
ggsave('scratch/figures/cohort_model_susceptibility_vs_dose_CoP_pre.png',units='in',width=6, height=3)

# Protective efficacy function
protective_efficacy = function(dose=1e4,  CoP_pre=1, outcome = 'fever_given_dose', CoP_control=1){
  VE = 1 - p_outcome_given_dose(dose,CoP_pre=CoP_pre,outcome = outcome)/
    p_outcome_given_dose(dose,CoP_pre=CoP_control,outcome = outcome)
  return(VE)
}

# Demo plot
expand.grid(dose = 10^seq(0,8,by=0.1),
            CoP_pre = round(10^seq(0,3.5,by=0.5)),
            outcome=factor(c('infection given dose','fever given dose'),
                           levels=c('infection given dose','fever given dose','fever given infection'))) |>
  group_by(outcome,CoP_pre,dose) |>
  mutate(protective_efficacy = protective_efficacy(dose=dose,CoP_pre=CoP_pre,outcome = gsub(' ','_',outcome))) |>
  mutate(CoP_pre = factor(CoP_pre)) |>
  ggplot() +
  geom_line(aes(x=dose,y=protective_efficacy,group=CoP_pre,color=CoP_pre)) +
  facet_grid('~outcome') +
  theme_bw() +
  ylim(c(0,1)) +
  scale_x_continuous(trans='log10', breaks=10^seq(0,10,by=2),minor_breaks = NULL,
                     labels = scales::trans_format("log10", math_format(10^.x)) ) +
  xlab('bacilli') + ylab('protective efficacy')
ggsave('scratch/figures/cohort_model_vaccine_efficacy_vs_dose_CoP_pre.png',units='in',width=4.5, height=3)


# =============================================================================
# COHORT INCIDENCE MODEL WITH VACCINATION
# =============================================================================

# KEY DIFFERENCE: vaccination_rate parameter and vaccination logic
cohort_model = function(exposure_dose,exposure_rate,
                        N=1000,age_max=75,
                        titer_dt=365/12,
                        max_titer_id = 1000,
                        exposure_rate_multiplier = c(rep(0.2,13),rep(0.5,12),rep(0.8,36),rep(1,120),rep(0.5,75*12-25-36-120)),
                        vaccination_rate = 1.0
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

  # expose everyone
  for (id in (1:N)){
    tmp_exposed = rpois(length(exposure_rate_multiplier),exposure_rate*exposure_rate_multiplier)
    events_list[[id]]$exposure_month = which(tmp_exposed>0)
    events_list[[id]]$exposure_count = tmp_exposed[events_list[[id]]$exposure_month]
    events_list[[id]]$infected = rep(0, length(events_list[[id]]$exposure_month))
    events_list[[id]]$fever = rep(0, length(events_list[[id]]$exposure_month))
  }

  # vaccinate (random age >=9 months)
  for (id in (1:N)){
    events_list[[id]]$vaccinated = FALSE
    if (runif(n=1) <= vaccination_rate){
      events_list[[id]]$vaccination_month = sample(seq(9,age_max*12-1,by=1),1)
    } else {
      events_list[[id]]$vaccination_month = max(simulation_months)+1
    }
  }

  # infect and update titers (with vaccination)
  for (id in (1:N)){
    tmp_titer =rep(1,length(simulation_months))

    if (length(events_list[[id]]$exposure_month) >0){
      for (k in 1:(length(events_list[[id]]$exposure_month)+1)){

        exposure_month = events_list[[id]]$exposure_month[k]
        vaccination_month = events_list[[id]]$vaccination_month

        # vaccinate if it's time
        if (!events_list[[id]]$vaccinated & vaccination_month <= max(simulation_months) &
            ((exposure_month >= vaccination_month) | k == (1+length(events_list[[id]]$exposure_month)))){

          titer_pre = tmp_titer[vaccination_month]
          titer_post = titer_pre * fold_rise_model(CoP_pre = titer_pre,source = 'tcv')

          future_times = vaccination_month:length(simulation_months)
          future_times_from_vaccination = titer_dt*(simulation_months-vaccination_month+1)[future_times]

          tmp_titer[future_times] = titer_vs_time(t=future_times_from_vaccination,
                                                  age_years = vaccination_month/12,
                                                  CoP_pre = titer_pre,
                                                  CoP_peak= titer_post)

          tmp_titer[vaccination_month:length(tmp_titer)] = round(tmp_titer[vaccination_month:length(tmp_titer)],2)
          events_list[[id]]$vaccinated=TRUE
        }

        if (k <= length(events_list[[id]]$exposure_month)){
          # infection probability
          p_once = p_outcome_given_dose(dose=exposure_dose,CoP_pre=tmp_titer[exposure_month],
                                        outcome='infection_given_dose')
          p_inf = 1-(1-p_once)^events_list[[id]]$exposure_count[k]

          if(runif(1)>p_inf){
            events_list[[id]]$infected[k] = 0
            events_list[[id]]$fever[k]    = 0
          } else {
            events_list[[id]]$infected[k]=1

            p_fever = p_outcome_given_dose(dose=exposure_dose,CoP_pre=tmp_titer[exposure_month],
                                           outcome='fever_given_infection')

            if (runif(n=1)<=p_fever){
              events_list[[id]]$fever[k]=1
            } else {
              events_list[[id]]$fever[k]=0
            }

            # update titer from infection
              titer_pre = tmp_titer[exposure_month]
              titer_post = titer_pre * fold_rise_model(CoP_pre = titer_pre)

              future_times = exposure_month:length(simulation_months)
              future_times_from_new_infection = titer_dt*(simulation_months-exposure_month+1)[future_times]

              tmp_titer[future_times] = titer_vs_time(t=future_times_from_new_infection,
                                                      age_years = exposure_month/12,
                                                      CoP_pre = titer_pre,
                                                      CoP_peak= titer_post)

            tmp_titer[exposure_month:length(tmp_titer)] = round(tmp_titer[exposure_month:length(tmp_titer)],2)
          }
        }
      }

      if(id <= max_titer_id){
        idx = titer_df$id == id
        titer_df$titer[idx] = tmp_titer
      }
    }
  }

  # tidy up events data
    not_empty_idx = which(!sapply(events_list, function(x){ is_empty(x$exposure_month)}))

    events_df = tibble(id = not_empty_idx, data = events_list[not_empty_idx]) |>
      unnest_wider(data) |>
      unnest_longer(c(exposure_month,exposure_count,infected,fever,vaccination_month), keep_empty = TRUE) |>
      mutate(vaccination_month = if_else(vaccination_month>max(simulation_months),NaN,vaccination_month)) |>
      mutate(exposure_age=exposure_month/12,
             vaccination_age=vaccination_month/12,
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

  # join titers with events
  titer_df = titer_df |>
    left_join(events_df |> select(id,infection_age,infected,fever) |>
                drop_na(infected),by=join_by(id ==id, age_years == infection_age))
    titer_df = titer_df |>
      left_join(events_df |> select(id,vaccination_age) |> filter(!duplicated(id)) |> drop_na(vaccination_age),
                by=join_by(id ==id)) |>
      mutate(vaccination_age = if_else(abs(vaccination_age-age_years)>0.01,NA,vaccination_age))

  return(list(titer_df=titer_df,
              events_df=events_df,
              config=list(exposure_dose=exposure_dose,exposure_rate=exposure_rate,
                          N=N,age_max=age_max,titer_dt=titer_dt,max_titer_id=max_titer_id,
                          vaccination_rate=vaccination_rate,
                          exposure_rate_multiplier=exposure_rate_multiplier)))
}


# =============================================================================
# RUN SIMULATIONS
# =============================================================================

if (!file.exists('scratch/output_cache_VE_highdose.RData')){
  output=list()

  N_cohort=1e6

  # medium incidence
  output[['medium']] = cohort_model(exposure_dose = 2.9e2,
                                    exposure_rate=1/(12*16),
                                    N=N_cohort,
                                    exposure_rate_multiplier = c(rep(0.25,13),rep(0.25,12),rep(1,36),rep(1,60),rep(1,60),rep(1,75*12-25-36-120)))

  # high incidence
  output[['high']] = cohort_model(exposure_dose = 2.9e2,
                                  exposure_rate=1/(12*1.9),
                                  N=N_cohort/5,
                                  exposure_rate_multiplier = c(rep(0.25,13),rep(0.25,12),rep(1,36),rep(1,60),rep(1,60),rep(1,75*12-25-36-120)))

  # very high incidence
  output[['very_high']] = cohort_model(exposure_dose = 3.1e3,
                                       exposure_rate=1/(12*1.9),
                                       N=N_cohort/10,
                                       exposure_rate_multiplier = c(rep(0.25,13),rep(0.25,12),rep(1,36),rep(1,60),rep(1,60),rep(1,75*12-25-36-120)))

  save(output,N_cohort,file='scratch/output_cache_VE_highdose.RData')
} else {
  load(file='scratch/output_cache_VE_highdose.RData')
}


# =============================================================================
# CALCULATE INCIDENCE
# =============================================================================

incidence_fever_targets = c(medium = 53,high=214,very_high=1255)

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
           symptomatic_fraction_overall = sum(symptomatic_fraction*age_width/sum(age_width)),
           incidence_fever_target = incidence_fever_targets[names(output)[k]])
}


# =============================================================================
# INCIDENCE PLOTS
# =============================================================================

incidence_fever_by_age_targets = rep(list(data.frame(age_group=unique(output[[1]]$incidence_vs_age$age_group),
                                                     incidence_fever = NaN,
                                                     lower=NaN,
                                                     upper=NaN)),3)
incidence_fever_by_age_targets[[1]]$incidence_fever <- c(22,50,65,72,35)
incidence_fever_by_age_targets[[2]]$incidence_fever = c(160,520,380,330,120)
incidence_fever_by_age_targets[[3]]$incidence_fever = c(1300,4200,2750,1900,600)

p_incidence_fever=list()
p_incidence_infection=list()
p_symptomatic_fraction=list()
for (k in 1:length(output)){
  p_incidence_fever[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_fever),stat='identity') +
    geom_point(data=incidence_fever_by_age_targets[[k]],aes(x=age_group,y=incidence_fever),stat='identity',color='orangered') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of fever per 100k') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),axis.title = element_text(size=10))

  p_incidence_infection[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_infection),stat='identity') +
    geom_hline(aes(yintercept=incidence_infection_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of infection per 100k') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),axis.title = element_text(size=10))

  p_symptomatic_fraction[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=symptomatic_fraction),stat='identity') +
    geom_hline(aes(yintercept=symptomatic_fraction_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('symptomatic fraction') +
    ylim(c(0,0.15)) +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),axis.title = element_text(size=10))
}

wrap_plots(p_incidence_fever) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_incidence_fever_by_age_VE.png',units='in',width=8,height=3)

wrap_plots(p_incidence_infection) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_incidence_infection_by_age_VE.png',units='in',width=8,height=3)

wrap_plots(p_symptomatic_fraction) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_symptomatic_fraction_VE.png',units='in',width=8,height=3)


# =============================================================================
# INDIVIDUAL-LEVEL DIAGNOSTICS WITH VACCINATION
# =============================================================================

p_exposure=list()
for (k in 1:length(output)){
  p_exposure[[k]] = ggplot(output[[k]]$events_df |> filter(id<=40)) +
                      geom_point(aes(x=exposure_age,color=outcome,y=id)) +
                      geom_point(aes(x=vaccination_age,color='vaccinated',y=id)) +
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
ggsave('scratch/figures/cohort_model_individual_level_exposure_examples_VE.png',units='in',width=7,height=6)

p_titer_examples=list()
for (k in 1:length(output)){
  p_titer_examples[[k]] = ggplot() +
    geom_line(data=output[[k]]$titer_df |> filter(id<=20),aes(x=age_years,y=titer)) +
    geom_point(data=output[[k]]$titer_df |> filter(id<=20) |> filter(!is.na(infected)),
               aes(x=age_years,y=titer,color=fever)) +
    geom_point(data=output[[k]]$titer_df |> filter(id<=20) |> filter(!is.na(vaccination_age)),
               aes(x=age_years,y=titer,color='vaccinated')) +
    scale_y_continuous(trans='log10',limits=c(1,10^4)) +
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
ggsave('scratch/figures/cohort_model_individual_level_titer_examples_VE.png',units='in',width=12,height=6)


# =============================================================================
# POPULATION IMMUNITY PLOTS
# =============================================================================

p_titer_summary=list()
p_seropositive_summary=list()
p_titer_density=list()
p_protective_efficacy_infection_summary=list()
p_protective_efficacy_fever_summary=list()

for (k in 1:length(output)){

  exposure_dose = output[[k]]$config$exposure_dose

  tmp_titer_summary =output[[k]]$titer_df |>
    mutate(protective_efficacy_infection = protective_efficacy(dose = exposure_dose,CoP_pre = titer, outcome='infection_given_dose')) |>
    mutate(protective_efficacy_fever = protective_efficacy(dose = exposure_dose,CoP_pre = titer, outcome='fever_given_dose')) |>
    mutate(seropositive_VaccZyme = 10^(log10(titer)*(1 + 0.37*rnorm(length(titer))))>7) |>
    group_by(age_years) |>
    summarize(titer_gmt = exp(mean(log(titer))),
           titer_median = median(titer),
           titer_upper_iqr = quantile(titer,probs=0.75),
           titer_lower_iqr = quantile(titer,probs=0.25),
           titer_upper_95 = quantile(titer,probs=0.975),
           titer_lower_95 = quantile(titer,probs=0.025),
           fraction_seropositive_VaccZyme = mean(seropositive_VaccZyme),
           protective_efficacy_infection_median = median(protective_efficacy_infection),
           protective_efficacy_infection_upper_iqr = quantile(protective_efficacy_infection,probs=0.75),
           protective_efficacy_infection_lower_iqr = quantile(protective_efficacy_infection,probs=0.25),
           protective_efficacy_fever_median = median(protective_efficacy_fever),
           protective_efficacy_fever_upper_iqr = quantile(protective_efficacy_fever,probs=0.75),
           protective_efficacy_fever_lower_iqr = quantile(protective_efficacy_fever,probs=0.25))

  p_titer_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_ribbon(aes(x=age_years,ymin=log2(titer_lower_95),ymax=log2(titer_upper_95)),alpha=0.2)+
    geom_ribbon(aes(x=age_years,ymin=log2(titer_lower_iqr),ymax=log2(titer_upper_iqr)),alpha=0.2)+
    geom_line(aes(x=age_years,y=log2(titer_median))) +
    geom_line(aes(x=age_years,y=log2(titer_gmt)),linetype='dashed') +
    scale_y_continuous(breaks=seq(0,10,by=1),minor_breaks=NULL)+
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),axis.title = element_text(size=10)) +
    ylab('anti-Vi IgG EU/ml')

  p_seropositive_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_line(aes(x=age_years,y=fraction_seropositive_VaccZyme)) +
    scale_y_continuous(limits=c(0,1)) +
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),axis.title = element_text(size=10)) +
    ylab('fraction seropositive')

  p_protective_efficacy_fever_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_ribbon(aes(x=age_years,ymin=protective_efficacy_fever_lower_iqr,ymax=protective_efficacy_fever_upper_iqr),alpha=0.2)+
    geom_line(aes(x=age_years,y=protective_efficacy_fever_median)) +
    scale_y_continuous(limits=c(0,1)) +
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),axis.title = element_text(size=10)) +
    ylab('protective efficacy against fever')

  sampled_titer_df = output[[k]]$titer_df |>
    group_by(age_group,id) |>
    slice_sample(n=1)
  p_titer_density[[k]] = ggplot(sampled_titer_df) +
    geom_density_ridges(aes(x=titer,y=age_group),scale=0.9,jittered_points = TRUE,bandwidth = 0.1) +
    scale_x_continuous(trans='log10',limits=c(7,10^3.5)) +
    theme_bw() +
    xlab('titer, given above LOD') +
    ylab('') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),axis.title = element_text(size=10))
}

wrap_plots(wrap_plots(p_titer_summary),wrap_plots(p_seropositive_summary),nrow=2) +
  plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_titer_vs_age_summary_VE.png',units='in',width=7,height=5)

wrap_plots(p_protective_efficacy_fever_summary) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_protective_efficacy_vs_age_summary_VE.png',units='in',width=7,height=2.5)

wrap_plots(p_titer_density) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_titer_density_VE.png',units='in',width=7,height=4)


# =============================================================================
# VACCINE EFFICACY STUDY SIMULATION
# =============================================================================

set.seed(42)
for (k in 1:3){
  N_pairs=output[[k]]$config$N

  tmp_events =  output[[k]]$events_df |>
    select(id,exposure_age,infected,fever,vaccination_age,age_group) |>
    mutate(vaccinated = exposure_age > vaccination_age) |>
    ungroup()

  # add subjects with no exposures
  if (length(setdiff(1:output[[k]]$config$N,tmp_events$id))>0){
    zeros_df = data.frame(id = setdiff(1:output[[k]]$config$N,tmp_events$id),
                          exposure_age=output[[k]]$config$age_max+1,
                          infected=0,fever=0) |>
      mutate(tmp_events |> select(vaccination_age,age_group,vaccinated) |> slice_sample(n=n()))
    tmp_events = tmp_events |> rbind(zeros_df)
  }

  # treatment arm: vaccinated fever cases
  treatment_set = tmp_events |> filter(vaccinated & fever==1) |>
    mutate(enrollment_age = vaccination_age) |>
    mutate(age_group=cut(enrollment_age,breaks=c(0,2,5,10,15,75),include.lowest = TRUE)) |>
    mutate(years_at_risk = exposure_age-enrollment_age) |>
    filter(years_at_risk>0) |>
    group_by(id) |>
    slice_head(n=1) |> ungroup()

  # control arm: unvaccinated fever cases
  control_set = tmp_events |> filter(!vaccinated & fever==1) |>
    group_by(id) |>
    mutate(enrollment_age = sample(seq(0,unique(vaccination_age),by=1/12),size=1)) |>
    mutate(age_group=cut(enrollment_age,breaks=c(0,2,5,10,15,75),include.lowest = TRUE)) |>
    ungroup() |>
    mutate(years_at_risk = exposure_age-enrollment_age) |>
    filter(years_at_risk>0) |>
    group_by(id) |> slice_head(n=1) |>
    ungroup()

  tmp_pairs = rbind(treatment_set,control_set) |> select(age_group, vaccinated, enrollment_age,years_at_risk)

  # calculate fever by 1-year intervals
  tmp_pairs = tmp_pairs |>
    mutate(fever_0.1=if_else(years_at_risk<=1,1,0),
          fever_1.2=if_else(years_at_risk<=1,NaN, if_else(years_at_risk>1 & years_at_risk<=2,1,0)),
          fever_2.3=if_else(years_at_risk<=2,NaN, if_else(years_at_risk>2 & years_at_risk<=3,1,0)),
          fever_3.4=if_else(years_at_risk<=3,NaN, if_else(years_at_risk>3 & years_at_risk<=4,1,0)),
          fever_4.5=if_else(years_at_risk<=4,NaN, if_else(years_at_risk>4 & years_at_risk<=5,1,0)),
          fever_5.6=if_else(years_at_risk<=5,NaN, if_else(years_at_risk>5 & years_at_risk<=6,1,0)),
          fever_6.7=if_else(years_at_risk<=6,NaN, if_else(years_at_risk>6 & years_at_risk<=7,1,0)),
          fever_7.8=if_else(years_at_risk<=7,NaN, if_else(years_at_risk>7 & years_at_risk<=8,1,0)),
          fever_8.9=if_else(years_at_risk<=8,NaN, if_else(years_at_risk>8 & years_at_risk<=9,1,0)),
          fever_9.10=if_else(years_at_risk<=9,NaN, if_else(years_at_risk>9 & years_at_risk<=10,1,0)),
          fever_10.11=if_else(years_at_risk<=10,NaN, if_else(years_at_risk>10 & years_at_risk<=11,1,0)),
          fever_11.12=if_else(years_at_risk<=11,NaN, if_else(years_at_risk>11 & years_at_risk<=12,1,0)),
          fever_12.13=if_else(years_at_risk<=12,NaN, if_else(years_at_risk>12 & years_at_risk<=13,1,0)),
          fever_13.14=if_else(years_at_risk<=13,NaN, if_else(years_at_risk>13 & years_at_risk<=14,1,0)),
          fever_14.15=if_else(years_at_risk<=14,NaN, if_else(years_at_risk>14 & years_at_risk<=15,1,0)),
          fever_15.16=if_else(years_at_risk<=15,NaN, if_else(years_at_risk>15 & years_at_risk<=16,1,0)),
          fever_16.17=if_else(years_at_risk<=16,NaN, if_else(years_at_risk>16 & years_at_risk<=17,1,0)),
          fever_17.18=if_else(years_at_risk<=17,NaN, if_else(years_at_risk>17 & years_at_risk<=18,1,0)),
          fever_18.19=if_else(years_at_risk<=18,NaN, if_else(years_at_risk>18 & years_at_risk<=19,1,0)),
          fever_19.20=if_else(years_at_risk<=19,NaN, if_else(years_at_risk>19 & years_at_risk<=20,1,0)),
          fever_20.21=if_else(years_at_risk<=20,NaN, if_else(years_at_risk>20 & years_at_risk<=21,1,0))) |>
    pivot_longer(cols=contains('fever_'),names_to='interval',
                            names_prefix='fever_',values_to = 'fever') |>
    mutate(interval = factor(interval, levels = c('0.1','1.2','2.3','3.4','4.5','5.6',
                                                  '6.7','7.8','8.9','9.10','10.11','11.12',
                                                  '12.13','13.14','14.15','15.16',
                                                  '16.17','17.18','18.19','19.20','20.21')))

  tmp_pairs = tmp_pairs |>
    mutate(interval_width = 1,
           interval_startpoint = seq(0,20,by=1)[as.numeric(interval)]) |>
    ungroup() |>
    mutate(years_at_risk = pmin(years_at_risk-interval_startpoint,interval_width)) |>
    mutate(years_at_risk = if_else(years_at_risk<=0,NaN,years_at_risk))

  # matched pairing
  controls <- tmp_pairs %>%
    filter(!vaccinated) %>%
    filter(!is.na(years_at_risk)) |>
    group_by(enrollment_age) %>%
    mutate(match_id = row_number())

  vaccinees <- tmp_pairs %>%
    filter(vaccinated) %>%
    filter(!is.na(years_at_risk)) |>
    group_by(enrollment_age) %>%
    mutate(match_id = row_number())

  matched_pairs <- inner_join(
    controls,
    vaccinees,
    by   = c("enrollment_age", "match_id"),
    suffix = c("_ctrl", "_vac")
  ) %>%
    arrange(enrollment_age, match_id)

  matched_long <- matched_pairs %>%
    select(match_id,enrollment_age,ends_with("_ctrl")) %>%
    rename_with(~ sub("_ctrl$", "", .x)) %>%
    bind_rows(
      matched_pairs %>%
        select(match_id,enrollment_age,ends_with("_vac")) %>%
        rename_with(~ sub("_vac$",  "", .x))
    ) %>%
    arrange(enrollment_age,match_id,interval,vaccinated)

  # calculate VE
  output[[k]]$VE_df = matched_long |> group_by(age_group,interval,vaccinated) |>
    summarize(frac_fever = sum(fever==1,na.rm=TRUE)/sum(years_at_risk,na.rm=TRUE)) |>
    group_by(age_group,interval) |>
    summarize(VE_fever = 1-frac_fever[vaccinated]/frac_fever[!vaccinated]) |>
    mutate(interval_midpoint = c(seq(1,21,by=1)-0.5)[as.numeric(interval)],
           interval_startpoint = c(seq(0,20,by=1))[as.numeric(interval)])
}

# combine VE results
VE_df = rbind(data.frame(incidence=names(output[1]),output[[1]]$VE_df),
              data.frame(incidence=names(output[2]),output[[2]]$VE_df),
              data.frame(incidence=names(output[3]),output[[3]]$VE_df)) |>
  mutate(incidence=factor(incidence,levels=c('medium','high','very_high')))

# fit logistic model to VE
mod = lm(log((0.00+VE_fever)/(1.00-VE_fever))~incidence*log(interval_midpoint) + incidence*age_group,
         data=VE_df)
summary(mod)

pred_data = expand_grid(incidence=factor(c('medium','high','very_high'),
                                         levels=c('medium','high','very_high')),
                        age_group = unique(VE_df$age_group),
                        interval_midpoint = seq(1/12,21,by=1/12))
pred_data = pred_data |>
  mutate(log_ve = predict(mod, newdata=pred_data)) |>
  mutate(VE_fever = 1/(1+exp(-log_ve)))

# VE by age group
VE_df |>
  ggplot() +
  geom_point(aes(y=VE_fever,x=interval_midpoint,color=incidence),size=0.5) +
  geom_step(aes(y=VE_fever,x=interval_midpoint-0.5,color=incidence),alpha=0.3) +
  geom_line(data = pred_data,
            aes(y=VE_fever,x=interval_midpoint,color=incidence),linewidth=0.75) +
  theme_bw() +
  facet_grid('~age_group') +
  ylim(c(0,1)) +
  xlim(c(0,20))+
  ylab('vaccine efficacy\nagainst fever') +
  xlab('years since enrollment') +
  theme(plot.title=element_text(size=10),axis.title = element_text(size=10))
ggsave('scratch/figures/cohort_model_VE_fever_by_age.png',units='in',width=7,height=2.5)

# VE by incidence level
ggplot(VE_df) +
  geom_point(aes(y=VE_fever,x=interval_midpoint,color=age_group),size=0.5) +
  geom_step(aes(y=VE_fever,x=interval_midpoint-0.5,color=age_group),alpha=0.3) +
  geom_line(data = pred_data, aes(y=VE_fever,x=interval_midpoint,color=age_group),linewidth=0.75) +
  theme_bw() +
  facet_grid('~incidence') +
  ylim(c(0,1)) +
  xlim(c(0,20))+
  ylab('vaccine efficacy\nagainst fever') +
  xlab('interval midpoint [years]') +
  theme(plot.title=element_text(size=10),axis.title = element_text(size=10))
ggsave('scratch/figures/cohort_model_VE_fever_incidence_level.png',units='in',width=7,height=2.5)
