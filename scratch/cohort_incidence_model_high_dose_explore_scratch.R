# cohort_incidence_model_high_dose_explore_scratch.R
#
# EXPLORATORY BRANCH - NOT THE CANONICAL VERSION
#
# This file diverged from scratch/cohort_incidence_model_proof_of_concept.R
# (canonical version in main branch) during high-dose exploration work.
#
# Key differences from canonical version:
# - Modified dose-response parameters for high-dose exploration
# - Experimental changes to exposure rate multipliers
# - Various parameter tweaks that may have broken logical consistency
# - Stripped of blog post formatting (roxygen comments)
#
# For the working, deployable model, see the canonical version in main branch.
# This file is preserved for reference and potential cherry-picking of ideas.
#
# =============================================================================

library(tidyverse)
library(patchwork)
library(ggridges)
library(scales)

# age width function
# Re-extracts age_group widths from a cut() call
age_width = function(age_group,age_max){
  round(as.numeric(as.character(fct_recode( age_group,
                                            `2`='[0,2]',`3`='(2,5]',`5`='(5,10]',`5.1`='(10,15]',!!as.character(age_max-15):=paste('(15,',age_max,']',sep='')))))
}

# =============================================================================
# INTRAHOST IMMUNITY MODEL FUNCTIONS
# =============================================================================

# Titer response function
# Describes time dynamics of the correlate of protection (CoP).
# - Initial baseline CoP = 1 (below typical LOD of 7 EU/ml)
# - Peak at ~30 days post-exposure
# - Power-law decay with age-dependent parameters fit to titer data
titer_vs_time = function(t,age_years, CoP_peak=1000, CoP_pre=1,
                         T_decay=11, alpha=1.16,
                         beta_T_decay_age=0.057, beta_alpha_age=-0.060,
                         CoP_min=1,T_rise=1.5,t_start=3,t_peak=30.4){

  Tda = T_decay * exp(beta_T_decay_age * age_years)
  ka = alpha * exp(beta_alpha_age * age_years)

  # power law decay
  titer = CoP_min + (CoP_peak-CoP_min)*(1+(t-t_peak)/(ka*Tda))^(-ka)

  # simple logistic interpolator rise (just for continuity, no model role)
  titer[t<t_peak ] =
    CoP_pre + (CoP_peak-CoP_pre)*
    (1/(1+exp(-(t[t<t_peak] - t_start*5)/T_rise)) - 1/(1+exp(-(0-t_start*5)/T_rise)))/
    (1/(1+exp(-(t_peak - t_start*5)/T_rise)) - 1/(1+exp(-(0-t_start*5)/T_rise)))

  return(titer)
}

# Demo plot: titer waning curves by age
expand.grid(t=seq(0,10,by=1/24), age_years = factor(c(1,5,15,45))) |>
  mutate(titer = titer_vs_time(t=t*365,age=as.numeric(age_years))) |>
  ggplot() +
  geom_line(aes(x=t,y=titer,color=age_years)) +
  theme_bw() + scale_y_continuous(trans='log10') + xlab('years post response') + ylab('')

# Fold-rise model
# Describes CoP response after immunizing event given pre-challenge CoP.
# - CoP_max ~10^3.5 from analysis of pre-post vaccine responses
# - mu_0 parameter loosely informed by infection dynamics models
fold_rise_model = function(CoP_pre,
                           mu_0=1.25,
                           CoP_max=10^3.5, CoP_min=1,
                           sigma_0=0.5,
                           response='individual'){# 'median'
  if(response == 'median'){
    mu = mu_0
  } else if (response=='individual'){
    mu = pmax(0,rnorm(length(CoP_pre),mean=mu_0,sd=sigma_0))
  }
  fold_rise = 10^(mu*(1-(log10(CoP_pre)-log10(CoP_min))/(log10(CoP_max)-log10(CoP_min))))
  return(fold_rise)
}

# Demo plot: fold-rise model
pl_df=expand.grid(CoP_pre=10^seq(0,3.5,by=0.1)) |>
  mutate(fold_rise = fold_rise_model(CoP_pre=CoP_pre, response='median')) |>
  mutate(CoP_post = fold_rise * CoP_pre)
(ggplot(pl_df) +
  geom_line(aes(x=CoP_pre,y=fold_rise)) +
  theme_bw() + scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') +
  xlab('pre-challenge titer') + ylab('fold-rise')) +
  (ggplot(pl_df) +
     geom_line(aes(x=CoP_pre,y=CoP_post)) +
     theme_bw() + scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') +
     xlab('pre-challenge titer') + ylab('post-challenge titer'))

# Example: titer trajectory over lifetime with two infections
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

# Dose-response as a function of immunity and bacilli ingested
# Based on challenge study data; parameters hand-tuned
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

# Demo plot: dose-response curves
plot_dat = expand.grid(dose = 10^seq(0,9,by=0.1),
            CoP_pre = c(50,round(10^seq(0,3.5,by=0.5))),
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
  geom_text(
    data = data.frame(
      x = 10^2,
      y = 0.3,
      label = "Hornick/\nLevine",
      outcome = factor("fever given dose",levels=levels(plot_dat$outcome))
    ),
    aes(x = x, y = y, label = label),fontface = "bold",size=3
  ) +
  labs(color = "anti-Vi titer\n[EU/ml]")
ggsave('scratch/figures/cohort_model_susceptibility_vs_dose_CoP_pre.png',units='in',width=6, height=3)

# Protective efficacy function
# Relative risk reduction vs naive (CoP=1) control
protective_efficacy = function(dose=1e4,  CoP_pre=1, outcome = 'fever_given_dose', CoP_control=1){
  VE = 1 - p_outcome_given_dose(dose,CoP_pre=CoP_pre,outcome = outcome)/
    p_outcome_given_dose(dose,CoP_pre=CoP_control,outcome = outcome)
  return(VE)
}

# Demo plot: protective efficacy vs dose
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
# COHORT INCIDENCE MODEL
# =============================================================================

# Main simulation function
# Constant force-of-infection cohort model demonstrating immunity over lifetime
cohort_model = function(exposure_dose,exposure_rate,
                        N=1000,age_max=75,
                        titer_dt=365/12, # monthly timesteps
                        max_titer_id = 1000, # limit titer tracking for performance
                        exposure_rate_multiplier = c(rep(0.2,13),rep(0.5,12),rep(0.8,36),rep(1,120),rep(0.5,75*12-25-36-120))
){
  # exposure_dose: bacilli per exposure event
  # exposure_rate: Poisson rate for exposures (per month)
  # exposure_rate_multiplier: age-dependent exposure adjustment

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

  # infect and update titers
  for (id in (1:N)){
    tmp_titer =rep(1,length(simulation_months))

    if (length(events_list[[id]]$exposure_month) >0){
      for (k in 1:length(events_list[[id]]$exposure_month)){

        exposure_month = events_list[[id]]$exposure_month[k]

        # infection probability
        p_once = p_outcome_given_dose(dose=exposure_dose,CoP_pre=tmp_titer[exposure_month],
                                      outcome='infection_given_dose')
        p_inf = 1-(1-p_once)^events_list[[id]]$exposure_count[k]

        if(runif(1)>p_inf){
          events_list[[id]]$infected[k] = 0
          events_list[[id]]$fever[k]    = 0
        } else {
          events_list[[id]]$infected[k]=1

          # fever probability given infection
          p_fever = p_outcome_given_dose(dose=exposure_dose,CoP_pre=tmp_titer[exposure_month],
                                         outcome='fever_given_infection')

          if (runif(n=1)<=p_fever){
            events_list[[id]]$fever[k]=1
          } else {
            events_list[[id]]$fever[k]=0
          }

          # update titer curve from infection forward
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

      if(id <= max_titer_id){
        idx = titer_df$id == id
        titer_df$titer[idx] = tmp_titer
      }
    }
  }

  # convert events_list to tidy dataframe
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

  # join infection events with titers for plotting
  titer_df = titer_df |>
    left_join(events_df |> select(id,infection_age,infected,fever) |>
                drop_na(infected),by=join_by(id ==id, age_years == infection_age))

  return(list(titer_df=titer_df,
              events_df=events_df,
              config=data.frame(exposure_dose,exposure_rate,N,age_max,titer_dt)))

}


# =============================================================================
# RUN SIMULATIONS
# =============================================================================

# Three WHO archetype scenarios: medium, high, very high incidence
if (TRUE | !file.exists('scratch/output_cache.RData')){
  output=list()

  N_cohort=1e6 # large for good stats at lower incidence

  # medium incidence: ~50/100k fever, flat age distribution
  output[['medium']] = cohort_model(exposure_dose = 4e2,
                                    exposure_rate=1/(12*20),
                                    N=N_cohort,
                                    exposure_rate_multiplier = c(rep(0.1,13),rep(0.5,12),rep(0.7,36),rep(1,60),rep(1,60),rep(0.5,75*12-25-36-120)))

  # high incidence: ~200/100k fever, peak at young ages
  output[['high']] = cohort_model(exposure_dose = 4e2,
                                  exposure_rate=1/(12*3),
                                  N=N_cohort/5,
                                  exposure_rate_multiplier = c(rep(0.1,13),rep(0.5,12),rep(1,36),rep(1,60),rep(0.7,60),rep(0.5,75*12-25-36-120)))

  # very high incidence: ~1000/100k fever, peak at young ages
  output[['very_high']] = cohort_model(exposure_dose = 5e3,
                                       exposure_rate=1/(12*3),
                                       N=N_cohort/10,
                                       exposure_rate_multiplier = c(rep(0.1,13),rep(0.5,12),rep(1,36),rep(0.7,60),rep(0.5,60),rep(0.5,75*12-25-36-120)))

  save(output,N_cohort,file='scratch/output_cache.RData')
} else {
  load(file='scratch/output_cache.RData')
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

# exposure rate by age (for reference)
exposure_rate_by_age_targets = expand.grid(age_group=unique(output[[1]]$incidence_vs_age$age_group),
                                           setting = names(output)) |>
  mutate(age_group_numeric = as.numeric(age_group)-0.5) |>
  mutate(exposures_per_year = c(1/20*c(0.3,0.7,1,1,0.5),
                                1/3*c(0.3,1,1,0.7,0.5),
                                1/3*c(0.3,1,0.7,0.5,0.5))) |>
  mutate(bacilli_per_exposure = c(4e2*c(1,1,1,1,1),
                              4e2*c(1,1,1,1,1),
                              5e3*c(1,1,1,1,1))) |>
  mutate(bacilli_per_year = bacilli_per_exposure * exposures_per_year) |>
  rbind(data.frame(age_group=NA,
                   age_group_numeric=5.5,
                   exposures_per_year=c(1/20*0.5,1/3*0.5,1/3*0.5),
                   bacilli_per_exposure = c(4e2,4e2,5e3),
                   bacilli_per_year=c(4e2/20*0.5,4e2/3*0.5,5e3/3*0.5),
                   setting = names(output)))

ggplot(exposure_rate_by_age_targets) +
  geom_step(aes(x=age_group_numeric,y=bacilli_per_year,group=setting)) +
  theme_bw() +
  facet_grid('~setting') +
  scale_x_continuous(breaks=1:5,labels = as.character(unique(output[[1]]$incidence_vs_age$age_group))) +
  scale_y_continuous(trans='log10',breaks=c(1,2,5,10,20,50,100,200,500,1000,2000),
                     minor_breaks = NULL) +
  xlab('age [years]') +
  ylab('bacilli\ningested\nper exposure') +
  theme(strip.text = element_blank())
ggsave('scratch/figures/cohort_model_mean_bacilli_ingested_per_year_by_age.png',units='in',width=8, height=1.25)

ggplot(exposure_rate_by_age_targets) +
  geom_step(aes(x=age_group_numeric,y=bacilli_per_exposure,group=setting)) +
  theme_bw() +
  facet_grid('~setting') +
  scale_x_continuous(breaks=1:5,labels = as.character(unique(output[[1]]$incidence_vs_age$age_group))) +
  scale_y_continuous(trans='log10',breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000),
                     minor_breaks = NULL) +
  xlab('age [years]') +
  ylab('bacilli\ningested\nper exposure') +
  theme(strip.text = element_blank())
ggsave('scratch/figures/cohort_model_mean_bacilli_ingested_per_exposure_by_age.png',units='in',width=8, height=1.25)

ggplot(exposure_rate_by_age_targets) +
  geom_step(aes(x=age_group_numeric,y=exposures_per_year,group=setting)) +
  theme_bw() +
  facet_grid('~setting') +
  scale_x_continuous(breaks=1:5,labels = as.character(unique(output[[1]]$incidence_vs_age$age_group))) +
  scale_y_continuous(minor_breaks = NULL) +
  theme(strip.text = element_blank()) +
  xlab('age [years]') +
  ylab('mean\nexposures\nper year')
ggsave('scratch/figures/cohort_model_fitted_exposure_rates_by_age.png',units='in',width=8, height=1.25)


# =============================================================================
# PLOTTING: INCIDENCE BY AGE
# =============================================================================

# fever incidence targets from WHO TCV working group
incidence_fever_by_age_targets = rep(list(data.frame(age_group=unique(output[[1]]$incidence_vs_age$age_group),
                                                     incidence_fever = NaN,
                                                     lower=NaN,
                                                     upper=NaN)),3)
incidence_fever_by_age_targets[[1]]$incidence_fever <- c(20,53,71,72,35)
incidence_fever_by_age_targets[[1]]$lower <- c(15,40,58,57,30)
incidence_fever_by_age_targets[[1]]$upper <- c(27,66,83,84,40)
incidence_fever_by_age_targets[[2]]$incidence_fever = c(160,500,420,275,140)
incidence_fever_by_age_targets[[2]]$lower <- c(120,440,370,240,160)
incidence_fever_by_age_targets[[2]]$upper <- c(200,560,470,320,120)
incidence_fever_by_age_targets[[3]]$incidence_fever = c(1400,4350,2950,1750,600)
incidence_fever_by_age_targets[[3]]$lower <- c(1200,4000,2750,1600,550)
incidence_fever_by_age_targets[[3]]$upper <- c(1600,4700,3150,1900,650)

p_incidence_fever=list()
p_incidence_infection=list()
p_symptomatic_fraction=list()
for (k in 1:length(output)){
  p_incidence_fever[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_fever),stat='identity') +
    geom_segment(data=incidence_fever_by_age_targets[[k]],aes(x=age_group,y=lower,yend=upper,group=age_group),stat='identity',color='orangered') +
    geom_point(data=incidence_fever_by_age_targets[[k]],aes(x=age_group,y=incidence_fever),stat='identity',color='orangered') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of fever per 100k') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))

  p_incidence_infection[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=incidence_infection),stat='identity') +
    geom_hline(aes(yintercept=incidence_infection_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('annual incidence of infection per 100k') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))

  p_symptomatic_fraction[[k]]=ggplot(output[[k]]$incidence_vs_age) +
    geom_bar(aes(x=age_group,y=symptomatic_fraction),stat='identity') +
    geom_hline(aes(yintercept=symptomatic_fraction_overall[1]),linetype='solid') +
    theme_bw() +
    xlab('') +
    ylab('symptomatic fraction') +
    ylim(c(0,0.2)) +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
}

wrap_plots(p_incidence_fever) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_incidence_fever_by_age.png',units='in',width=8,height=3)

wrap_plots(p_incidence_infection) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_incidence_infection_by_age.png',units='in',width=8,height=3)

wrap_plots(p_symptomatic_fraction) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_symptomatic_fraction.png',units='in',width=8,height=3)


# =============================================================================
# PLOTTING: INDIVIDUAL-LEVEL DIAGNOSTICS
# =============================================================================

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


# =============================================================================
# PLOTTING: POPULATION IMMUNITY
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
           protective_efficacy_infection_upper_95 = quantile(protective_efficacy_infection,probs=0.975),
           protective_efficacy_infection_lower_95 = quantile(protective_efficacy_infection,probs=0.025),
           protective_efficacy_fever_median = median(protective_efficacy_fever),
           protective_efficacy_fever_upper_iqr = quantile(protective_efficacy_fever,probs=0.75),
           protective_efficacy_fever_lower_iqr = quantile(protective_efficacy_fever,probs=0.25),
           protective_efficacy_fever_upper_95 = quantile(protective_efficacy_fever,probs=0.975),
           protective_efficacy_fever_lower_95 = quantile(protective_efficacy_fever,probs=0.025))

  p_titer_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_ribbon(aes(x=age_years,ymin=log2(titer_lower_95),ymax=log2(titer_upper_95)),alpha=0.2)+
    geom_ribbon(aes(x=age_years,ymin=log2(titer_lower_iqr),ymax=log2(titer_upper_iqr)),alpha=0.2)+
    geom_line(aes(x=age_years,y=log2(titer_median))) +
    geom_line(aes(x=age_years,y=log2(titer_gmt)),linetype='dashed') +
    scale_y_continuous(breaks=seq(0,10,by=1),minor_breaks=NULL)+
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) +
    ylab('anti-Vi IgG EU/ml')

p_seropositive_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_line(aes(x=age_years,y=fraction_seropositive_VaccZyme)) +
    scale_y_continuous(limits=c(0,1)) +
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) +
    ylab('fraction seropositive\n[anti-Vi IgG EU/ml > 7]')

p_protective_efficacy_infection_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_ribbon(aes(x=age_years,ymin=protective_efficacy_infection_lower_95,ymax=protective_efficacy_infection_upper_95),alpha=0.2)+
    geom_ribbon(aes(x=age_years,ymin=protective_efficacy_infection_lower_iqr,ymax=protective_efficacy_infection_upper_iqr),alpha=0.2)+
    geom_line(aes(x=age_years,y=protective_efficacy_infection_median)) +
    scale_y_continuous(limits=c(0,1)) +
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) +
    ylab('protective efficacy against infection')

p_protective_efficacy_fever_summary[[k]] = ggplot(tmp_titer_summary) +
    geom_ribbon(aes(x=age_years,ymin=protective_efficacy_fever_lower_95,ymax=protective_efficacy_fever_upper_95),alpha=0.2)+
    geom_ribbon(aes(x=age_years,ymin=protective_efficacy_fever_lower_iqr,ymax=protective_efficacy_fever_upper_iqr),alpha=0.2)+
    geom_line(aes(x=age_years,y=protective_efficacy_fever_median)) +
    scale_y_continuous(limits=c(0,1)) +
    theme_bw() +
    xlab('age') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10)) +
    ylab('protective efficacy against fever')


  sampled_titer_df = output[[k]]$titer_df |>
    group_by(age_group,id) |>
    slice_sample(n=1)
  p_titer_density[[k]] = ggplot(sampled_titer_df) +
    geom_density_ridges(aes(x=titer,y=age_group),scale=0.9,jittered_points = TRUE,bandwidth = 0.1) +
    scale_x_continuous(trans='log10',limits=c(7,10^3.5)) +
    theme_bw() +
    xlab('titer, given above limit of detection') +
    ylab('') +
    labs(title = paste(sub('_',' ',names(output)[k]),' incidence',sep=''),
         subtitle=paste('dose = ',output[[k]]$config$exposure_dose,' bacilli\nmean years b/w exposures = ',
                        1/(12*output[[k]]$config$exposure_rate),sep='')) +
    theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
          axis.title = element_text(size=10))
}

wrap_plots(wrap_plots(p_titer_summary),wrap_plots(p_seropositive_summary),nrow=2) +
  plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_titer_vs_age_summary.png',units='in',width=7,height=5)

wrap_plots(wrap_plots(p_protective_efficacy_fever_summary),wrap_plots(p_protective_efficacy_infection_summary),nrow=2) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_protective_efficacy_vs_age_summary.png',units='in',width=7,height=4)

wrap_plots(p_titer_density) + plot_layout(guides = "collect",axes='collect')
ggsave('scratch/figures/cohort_model_titer_density.png',units='in',width=7,height=4)
