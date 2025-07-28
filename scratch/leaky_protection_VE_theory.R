# leaky protection VE theory

library(tidyverse)

# basic version, constant risks, 

# medium/high dose
p_infection_no_vax = 0.7
p_fever_given_infection_no_vax = 0.7

p_infection_vax = 0.3
p_fever_given_infection_vax = 0.3


# very high dose
p_infection_no_vax = 0.85
p_fever_given_infection_no_vax = 0.6

p_infection_vax = 0.6
p_fever_given_infection_vax = 0.2



# single exposure
VE_single_exposure = 1 - (p_infection_vax * p_fever_given_infection_vax)/
                            (p_infection_no_vax*p_fever_given_infection_no_vax)
VE_single_exposure


# poisson exposure: sum over geometric distribution with poisson weights
lambda = 1/7 # medium
lambda = 1/0.8 # high
RR_vax = 0
RR_no_vax = 0
for (k in 0:(lambda + 10*(sqrt(lambda)))){
  RR_vax = RR_vax + 
    # probability of a fever when infected upon an exposure
    (p_infection_vax * p_fever_given_infection_vax) *
    # probability of no fever on exposure k: either no infection or infected but no fever
    (1 - p_infection_vax*p_fever_given_infection_vax)^k  *
    # poisson probability
    exp(-lambda)*lambda^k/gamma(k+1)
  
  RR_no_vax = RR_no_vax + (p_infection_no_vax * p_fever_given_infection_no_vax)*
    (1 - p_infection_no_vax*p_fever_given_infection_no_vax)^k *
    exp(-lambda)*lambda^k/gamma(k+1)
  
}
VE_lambda_expsoures = 1- RR_vax/RR_no_vax
VE_lambda_expsoures
