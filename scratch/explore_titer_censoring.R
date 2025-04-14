# explore true negative vs below limit of detection



library(tidyverse)
library(ggridges)
library(MASS)
library(stats4)


# efficacy data
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

ggplot(d_titer,aes(y=day_factor,x=Vi_IgG_U_per_ml,color=age_label)) +
  facet_grid('day~group') +
  geom_density_ridges(jittered_points = TRUE, fill = NA, scale = 0.85,alpha=0.1) +
  theme_bw() +
  scale_x_continuous(trans='log10')


# best signle guess at a once-infected-ish group is 0.75-2 and 2-4  controls + day 0 treatments
d_titer = d_titer |> mutate(probably_infected_once_group = 
                              (age_min_years<5 & !(group=='treatment' & day==28)))

ggplot(d_titer,aes(y=day_factor,x=Vi_IgG_U_per_ml,color=probably_infected_once_group)) +
  facet_grid('group~age_label') +
  geom_density_ridges(jittered_points = TRUE, fill = NA, scale = 0.85,alpha=0.1) +
  theme_bw() +
  scale_x_continuous(trans='log10')

ggplot(d_titer |> filter(probably_infected_once_group==TRUE),
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


TLN_likelihood = function(x,mu,sd,lower=log10(7),upper=Inf){
  sum(dtnorm(x, mu, sd,lower,upper))
}

log10_data=log10(d_titer$Vi_IgG_U_per_ml[d_titer$probably_infected_once_group & d_titer$Vi_IgG_U_per_ml>7])
mean(log10_data)
sd(log10_data)
TLN_likelihood(x=log10_data,mu=1.44,sd=0.55)

mod = fitdistr( x=log10_data, 
              TLN_likelihood, 
              method="L-BFGS-B",
              start=list(mu=1.44,sd=0.54),
              lower=c(0,0.3),upper=c(5,3),
              hessian=TRUE)

TLN_likelihood = function(mu,sd){
  sum(dtnorm(log10_data, mu, sd,lower=log10(7),upper=Inf))
}

mod = mle(
                TLN_likelihood, 
                method="L-BFGS-B",
                start=list(mu=1.44,sd=0.54),
                lower=c(-Inf,0.01),upper=c(Inf,3))

summary(mod)

hist(log10_data)
