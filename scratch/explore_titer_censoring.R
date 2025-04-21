# explore true negative vs below limit of detection



library(tidyverse)
library(ggridges)
library(MASS)
library(stats4)


# individual-level titer data
d_titer = readxl::read_excel('data/typhoid_vaccine_study_data.xlsx',sheet = 'Qadri2021_pre_post_given_positi')  |>
  mutate(day=as.numeric(day),
         Vi_IgG_U_per_ml=as.numeric(Vi_IgG_U_per_ml),
         age_min_years = as.numeric(age_lower),
         age_max_years = as.numeric(age_upper)) |>
  mutate(age_mean=1/2*(age_min_years + age_max_years)) |>
  mutate(age_label = factor(paste(age_min_years,'-',age_max_years,sep=''))) |>
  mutate(age_label = reorder(age_label, age_mean)) |>
  mutate(day_factor=factor(day,levels=sort(unique(day)))) |>
  mutate(grouping_var = interaction(age_label,group)) |>
  uncount(N)

names(d_titer)


ggplot(d_titer,aes(x=day_factor,y=Vi_IgG_U_per_ml)) +
  facet_grid('~grouping_var') +
  geom_jitter(width=0.25) +
  theme_bw() +
  scale_y_continuous(trans='log10')

ggplot(d_titer,aes(y=day_factor,x=Vi_IgG_U_per_ml)) +
  facet_grid('group~age_label') +
  geom_density_ridges(jittered_points = TRUE, fill = NA, scale = 0.85,alpha=0.1) +
  theme_bw() +
  scale_x_continuous(trans='log10')

p1=ggplot(d_titer |> filter(Vi_IgG_U_per_ml>7),
       aes(y=day_factor,x=Vi_IgG_U_per_ml)) +
  facet_grid('group~age_label') +
  geom_density_ridges(jittered_points = TRUE, fill = NA, scale = 0.85, bandwidth=0.15) +
  theme_bw() +
  scale_x_continuous(trans='log10')
p1

p2=ggplot(d_titer |> filter(Vi_IgG_U_per_ml>7) |> filter(day ==0),
          aes(y=age_label,x=Vi_IgG_U_per_ml)) +
  geom_density_ridges(jittered_points = TRUE, scale = 0.85, bandwidth=0.1) +
  theme_bw() +
  scale_x_continuous(trans='log10',limits=c(7,10^3.5))
p2
ggsave('scratch/figures/quadri2021_titer_density_given_detection_day0.png',units='in',width=4,height=3)

# best single guess at a once-infected-ish group is 0.75-2 and 2-4  controls + day 0 treatments
d_titer = d_titer |> mutate(probably_infected_once_group = 
                              (age_min_years<5 & !(group=='treatment' & day==28)))

ggplot(d_titer |> filter(Vi_IgG_U_per_ml>7),
       aes(y=day_factor,x=Vi_IgG_U_per_ml,color=probably_infected_once_group)) +
  facet_grid('group~age_label') +
  geom_density_ridges(jittered_points = TRUE, fill = NA, scale = 0.85,alpha=0.1) +
  theme_bw() +
  scale_x_continuous(trans='log10')

ggplot(d_titer |> filter(probably_infected_once_group==TRUE) |> filter(Vi_IgG_U_per_ml>7),
       aes(y=probably_infected_once_group,x=Vi_IgG_U_per_ml))+
  geom_density_ridges(jittered_points = TRUE, fill = NA, scale = 0.85,alpha=0.1) +
  theme_bw() +
  scale_x_continuous(trans='log10')



# defined left-truncated normal distribution
dtnorm<- function(x, mean, sd, lower=0, upper=Inf) {
  tmp = dnorm(x, mean, sd)/(pnorm(upper, mean, sd)-pnorm(lower, mean, sd))
  tmp[x<lower]=0
  return(tmp)
}
ptnorm <- function(x, mean, sd, lower=0, upper=Inf) {
  (pnorm(x,mean,sd) - pnorm(lower,mean,sd)) / 
    (pnorm(upper,mean,sd) - pnorm(lower,mean,sd))
}

TLN_likelihood = function(x,mu,sd,lower=log10(7),upper=100){
  dtnorm(x=x, mean=mu, sd=sd,lower,upper)
}

log10_data=log10(d_titer$Vi_IgG_U_per_ml[
  d_titer$probably_infected_once_group & 
    d_titer$Vi_IgG_U_per_ml>7 &
    d_titer$Vi_IgG_U_per_ml<100])
mean(log10_data)
sd(log10_data)

mod = fitdistr( x=log10_data, 
              TLN_likelihood, 
              method="L-BFGS-B",
              start=list(mu=1.34,sd=0.36),
              lower=c(-Inf,0.3),upper=c(Inf,3),
              hessian=FALSE)
mod


fit_dat = data.frame(Vi_IgG_U_per_ml=10^seq(-1,3,by=0.1)) |>
  mutate(model_censored = dtnorm(x=log10(Vi_IgG_U_per_ml),mean=mod$estimate[1],sd=mod$estimate[2],lower=log10(7))) |>
  mutate(model = dnorm(x=log10(Vi_IgG_U_per_ml),mean=mod$estimate[1],sd=mod$estimate[2]))

ggplot()+
  geom_density_ridges(data = d_titer |> filter(probably_infected_once_group==TRUE) |> 
                        filter(Vi_IgG_U_per_ml>=7 & Vi_IgG_U_per_ml<300), 
                      aes(y=probably_infected_once_group,x=Vi_IgG_U_per_ml),
                      jittered_points = TRUE, fill = NA, scale = 1,alpha=0.1) +
  geom_line(data=fit_dat,aes(x=Vi_IgG_U_per_ml,y=1+model_censored*0.88),color='red') +
  geom_line(data=fit_dat,aes(x=Vi_IgG_U_per_ml,y=1+model*1.33),color='blue') +
  theme_bw() +
  scale_x_continuous(trans='log10')

p1

# probability positive but below detection
pnorm(q=log10(7),mean=mod$estimate[1],sd=mod$estimate[2], lower.tail=TRUE)

# looks like the assay misses about a third of the positives



# quick check of variance on treated vs controls, to get a feel for what's assay variability vs real human diversity

log10_data=log10(d_titer$Vi_IgG_U_per_ml[(d_titer$group=='treatment' & d_titer$day==28)])
log10_data=log10_data[log10_data>log10(300)]
mean(log10_data)
sd(log10_data)

tmp=function(x,mu,sd){TLN_likelihood(x,mu,sd, upper=Inf)}
tmp(x=log10_data,mu=3.5,sd=0.54)
mod2 = fitdistr( x=log10_data, 
                tmp, 
                method="L-BFGS-B",
                start=list(mu=3.5,sd=0.54),
                lower=c(-Inf,0.3),upper=c(Inf,3),
                hessian=FALSE)
mod2
10^3.5 # true median peak titer
mod

# so basically, the implied variance near the threshold of detection, ater accounting for censoring,
# and the variance after boosting, are the same. 
# To a first approximation, this implies to me that what we're seeing with the variance is dominated by
# assay variation and not biology.
# To a second, insofar as the variance does appear to get wider with age, that is showing the biology effect,
# as the life histories of people will affect individual-level antibody expression and boosting etc. 