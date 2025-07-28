# fix dose response model to be properly monotonic

library(tidyverse)
library(stats4)

# current parameterization with the wrong conditional
p_outcome_given_dose_wrong = function(dose=1e4,  CoP_pre=1, outcome = 'fever_given_dose', 
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


# corrected parameterization with the right conditional
p_outcome_given_dose = function(dose=1e4,  CoP_pre=1, outcome = 'fever_given_dose', 
                                      n50_fever_given_infection=26250, alpha_fever_given_infection=0.79, 
                                      gamma_fever_given_infection=0.39,
                                      n50_infection_given_dose=27800/40,alpha_infection_given_dose = 0.84*2, 
                                      gamma_infection_given_dose=0.4/2 
){
  
  if(outcome == 'fever_given_infection'){
    
    p = 1 - (1+dose*(2^(1/alpha_fever_given_infection)-1)/n50_fever_given_infection)^(-alpha_fever_given_infection/(CoP_pre^gamma_fever_given_infection))
    
  } else if (outcome == 'infection_given_dose'){
    
    p = 1 - (1+dose*(2^(1/alpha_infection_given_dose)-1)/n50_infection_given_dose)^
      (-alpha_infection_given_dose/(CoP_pre^gamma_infection_given_dose))
    
  } else if (outcome == 'fever_given_dose'){
    
    p = (1 - (1+dose*(2^(1/alpha_infection_given_dose)-1)/n50_infection_given_dose)^
         (-alpha_infection_given_dose/(CoP_pre^gamma_infection_given_dose))) *
      (1 - (1+dose*(2^(1/alpha_fever_given_infection)-1)/n50_fever_given_infection)^(-alpha_fever_given_infection/(CoP_pre^gamma_fever_given_infection)))
  }
  
  return(p)
}


## descriptive plot before refitting
plot_dat = expand.grid(dose = 10^seq(0,9,by=0.1),
                       CoP_pre = c(50,round(10^seq(0,3.5,by=0.5))),
                       outcome=factor(c('infection given dose','fever given dose','fever given infection'),
                                      levels=c('infection given dose','fever given dose','fever given infection'))) |>
  group_by(outcome,CoP_pre,dose) |>
  mutate(probability_wrong = p_outcome_given_dose_wrong(dose=dose,CoP_pre=CoP_pre,outcome = gsub(' ','_',outcome))) |>
  mutate(probability = p_outcome_given_dose(dose=dose,CoP_pre=CoP_pre,outcome = gsub(' ','_',outcome),
                                            n50_fever_given_infection=27800, alpha_fever_given_infection=0.84, 
                                            gamma_fever_given_infection=0.4)) |>
  mutate(CoP_pre = factor(CoP_pre))

ggplot() +
  geom_line(data = plot_dat |> filter(CoP_pre!=50),
            aes(x=dose,y=probability_wrong,group=CoP_pre,color=CoP_pre),linetype='dashed') +
  geom_line(data = plot_dat |> filter(CoP_pre!=50),
            aes(x=dose,y=probability,group=CoP_pre,color=CoP_pre)) +
  geom_line(data = plot_dat |> filter(CoP_pre==50 & outcome == 'fever given dose'),
            aes(x=dose,y=probability,group=CoP_pre),color='black',linewidth=1) +
  facet_grid('~outcome') +
  theme_bw() +
  scale_x_continuous(trans='log10', breaks=10^seq(0,10,by=2),minor_breaks = NULL,
                     labels = trans_format("log10", math_format(10^.x)) ) +
  scale_y_continuous(trans='logit') +
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



# fit correct fever given dose form to incorrect form


mlogl = function(n50_fever_given_infection=27800, alpha_fever_given_infection=0.84, 
                         gamma_fever_given_infection=0.4){
  
  fit_dat = expand.grid(dose = 10^seq(0,9,by=0.1), # typical dosing # default used in Typhoidsim, from https://qmrawiki.canr.msu.edu/experiments/salmonella-typhi recommended model
                        CoP_pre = (10^seq(0,3.5,by=0.1)),
                        outcome=c('fever_given_dose')) |>
    group_by(outcome,CoP_pre,dose) |>
    mutate(probability_wrong = p_outcome_given_dose_wrong(dose=dose,CoP_pre=CoP_pre,outcome=outcome)) |>
    mutate(probability = p_outcome_given_dose(dose=dose,CoP_pre=CoP_pre,outcome=outcome,
                                              n50_fever_given_infection=n50_fever_given_infection, 
                                              alpha_fever_given_infection=alpha_fever_given_infection, 
                                              gamma_fever_given_infection=gamma_fever_given_infection)) 
  
  # mlogl = sum((log(fit_dat$probability/(1-fit_dat$probability)) - 
  #                 log(fit_dat$probability_wrong/(1-fit_dat$probability_wrong)))^2)
  mlogl = sum((fit_dat$probability - 
                 fit_dat$probability_wrong)^2)
  return(mlogl)
}
mlogl()

mod <- mle(mlogl, start = list(n50_fever_given_infection=27000, alpha_fever_given_infection=0.84, 
                                    gamma_fever_given_infection=0.4),
        method="L-BFGS-B",
        lower= c(0.1,0.1,0.1),
        upper=c(1e6,2,2),
        control=list(parscale = c(100,0.1,0.1)))
summary(mod)
round(coef(mod),2)

mlogl(n50_fever_given_infection = coef(mod)[1],
      alpha_fever_given_infection= coef(mod)[2],
      gamma_fever_given_infection= coef(mod)[3])

## descriptive plot after fit
plot_dat = expand.grid(dose = 10^seq(0,9,by=0.1),
                       CoP_pre = c(50,round(10^seq(0,3.5,by=0.5))),
                       outcome=factor(c('infection given dose','fever given dose','fever given infection'),
                                      levels=c('infection given dose','fever given dose','fever given infection'))) |>
  group_by(outcome,CoP_pre,dose) |>
  mutate(probability_wrong = p_outcome_given_dose_wrong(dose=dose,CoP_pre=CoP_pre,outcome = gsub(' ','_',outcome))) |>
  mutate(probability = p_outcome_given_dose(dose=dose,CoP_pre=CoP_pre,outcome = gsub(' ','_',outcome),
                                            n50_fever_given_infection = coef(mod)[1],
                                            alpha_fever_given_infection= coef(mod)[2],
                                            gamma_fever_given_infection= coef(mod)[3])) |>
  mutate(CoP_pre = factor(CoP_pre))

ggplot() +
  geom_line(data = plot_dat |> filter(CoP_pre!=50),
            aes(x=dose,y=probability_wrong,group=CoP_pre,color=CoP_pre),linetype='dashed') +
  geom_line(data = plot_dat |> filter(CoP_pre!=50),
            aes(x=dose,y=probability,group=CoP_pre,color=CoP_pre)) +
  geom_line(data = plot_dat |> filter(CoP_pre==50 & outcome == 'fever given dose'),
            aes(x=dose,y=probability,group=CoP_pre),color='black',linewidth=1) +
  facet_grid('~outcome') +
  theme_bw() +
  # ylim(c(0,1)) +
  scale_x_continuous(trans='log10', breaks=10^seq(0,10,by=2),minor_breaks = NULL,
                     labels = trans_format("log10", math_format(10^.x)) ) +
  # scale_y_continuous(trans='logit',breaks=c(0.01,0.1,0.2,0.5,0.7,0.9,0.99)) +
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



# protective efficacy vs CoP_pre
protective_efficacy = function(dose=1e4,  CoP_pre=1, outcome = 'fever_given_dose', CoP_control=1){
  VE = 1 - p_outcome_given_dose(dose,CoP_pre=CoP_pre,outcome = outcome)/
    p_outcome_given_dose(dose,CoP_pre=CoP_control,outcome = outcome)
  return(VE)
}

expand.grid(dose = 10^seq(0,8,by=0.1),
            CoP_pre = round(10^seq(0,3.5,by=0.5)),
            outcome=factor(c('infection given dose','fever given dose','fever given infection'),
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
