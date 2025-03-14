# /* empirical_bayes_random_effects_via_optimization scratch */

#' ---
#' title: "Hand rolling empirical Bayes estimation of a hierarchical model to learn how it works"
#' output:
#'   md_document
#' ---
#'
#' # Introduction
#' 
#' This is a quick experiment to teach myself if I can just use full likelihood optimization
#' to get well-behaved empirical bayes estimate of both the trial-level random effects
#' and the metastudy group-level hyperparameters in a mixed model. 
#' 
#' **Why?** Because statisticians tell me you aren't supposed to use a full likelihood 
#' to estimate random effects and hyperparameters because that is biased for the group-level variance components,
#' but do I really care for model fitting metastudy applications? 
#' 
#' It matters to me because it's nice to use
#' simple optimizations when initially calibrating models on the kind of strange data types you get from papers when you don't have
#' the individual-level data. I want to better understand when I can take shortcuts and still get decent statistical properties
#' during model calibration.
#' 
#' **What?** Compare full maximum likelihood (ML) estimate where I jointly optimize over
#' the subject-level random effects and the group level hyperparameters, vs
#' the restricted maximum likelihood (REML) estimate where I first optimize only over the hyperparamaters 
#' and then estimate the random effects given the hyperparameters.  For the best tight description of the two algorithms that I've found,
#' see this documentation for 
#' [Estimating Parameters in Linear Mixed-Effects Models from The Mathworks](https://www.mathworks.com/help/stats/estimating-parameters-in-linear-mixed-effects-models.html). 
#' 
#' **What did I learn?** WHELP, it looks full likelihood optimization is just fine, at least for this example, which is on a relevant scale 
#' for a vaccine model building metastudy. (Of order 10 trial arms, with of order 100 subjects per arm, and binary outcomes.)
#' 
#' To learn more, and very quickly, about what I'm trying to do, here's a 
#' [useful chat log with ChatGPT4.5](https://chatgpt.com/share/e/67d487c3-f080-800f-b14f-ce399681d9bc) to go along with this. 
#'
#' # Setting up the study
#' 
#' First, we set up a study with 10 trial arms and 100 people measured per arm. 
#' The observations for each trial arm are zero or one for the outcome, and each arm
#' has a true probability of the outcome drawn from a beta distribution for the ensemble of studies.
#' Note that this isn't set up as real vaccine efficacy measurement. I'm just interested in looking at
#' cohort parameter estimation for now. VE would be a transform on that, if the arms were labeled by treatment vs control.


#+ echo=TRUE, message=FALSE
# set up the environment
library(tidyverse)
library(rmutil)
library(stats4)
library(knitr)

# seed
set.seed(100)

##simulate binomial draws from a beta random-effects model
# typical example of a metastudy, 10 trials, ~100 subjects per trial
n_trials = 10 # don't change! I hard-coded 10 later for dumb reasons I didn't feel like debugging.
n_subjects=100

# set up the measurement-level data frame
subjects = expand_grid(ID=1:n_trials,rep = 1:n_subjects)

# draw probability of positive for each ID from a beta distribution

# rmutil::betabinomial parameterization
m=0.3
s=8

# rbeta parameterization
alpha=m*s
beta=s-alpha

# add true trial-level outcome probality to the subject data
subjects = subjects |> left_join(
  data.frame(ID=1:n_trials,p=rbeta(n_trials,shape1=alpha,shape2=beta)))

# draw binomial samples for each subject
subjects = subjects |> cbind(
  data.frame(response=rbinom(n=nrow(subjects),size=1,prob=subjects$p)))

# collapse the measurement-level data into subject level data for fitting later
observed = subjects |> group_by(ID,p) |>
            summarize(n_pos=sum(response),
                      n_trials=n(),
                      p_hat = n_pos/n_trials)
#+ message=TRUE
observed |> kable()

#' # Parameter inference with the full likelihood

minus_log_lik = function(m=alpha/(alpha+beta),s=(alpha+beta),
                         p1=observed$p[1],p2=observed$p[2],p3=observed$p[3],p4=observed$p[4], # there's probably a nicer way to do this...
                         p5=observed$p[5],p6=observed$p[6],p7=observed$p[7],p8=observed$p[8],
                         p9=observed$p[9],p10=observed$p[10]){
  a = m*s
  b = s-a
  p =c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
  
  -sum(dbeta(x=p,shape1=a,shape2=b,log=TRUE) + dbinom(x=observed$n_pos,size=observed$n_trials,prob=p,log=TRUE))
} 
minus_log_lik()

minus_log_REML()/minus_log_lik() # not sure why these aren't identical... I guess it must not integrating over the p's is like overfitting. 

model_lik = stats4::mle(minus_log_lik, start = list(m=0.5,s=10,p1=0.5,p2=0.5,p3=0.5,p4=0.5,
                                                    p5=0.5,p6=0.5,p7=0.5,p8=0.5,p9=0.5,p10=0.5),
                        method='L-BFGS-B',
                        upper=c(0.99999,Inf,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999),
                        lower=c(0.00001,0,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001))
summary(model_lik)


full_lik_params = data.frame(true=c(m=m,s=s,observed$p), estimate =coef(model_lik),  se = sqrt(diag(vcov(model_lik )))) |>
  mutate( z = (true-estimate)^2/se^2,
          lower = estimate - 1.96*se,
          upper = estimate + 1.96*se) |>
  rownames_to_column(var='param') |>
  select(param, everything()) |> 
  mutate(param=factor(param,levels=c('m','s','p1','p2','p3','p4','p5','p6','p7','p8','p9','p10'))) |>
  mutate(level=factor(c('group_mean','group_dispersion',rep('subject',10)))) |>
  mutate(param_index = 1:12)

full_lik_params


## estimate shape parameters with random effects "p" integrated out
# REML

minus_log_REML = function(m=alpha/(alpha+beta),s=(alpha+beta)){
  -sum(rmutil::dbetabinom(y=observed$n_pos,size=observed$n_trials,m=m,s=s,log=TRUE))
}
minus_log_REML()

model_REML = stats4::mle(minus_log_REML, start = list(m=0.5,s=1),
            method='L-BFGS-B',
            upper=c(0.99999,Inf),
            lower=c(0.00001,0))

summary(model_REML)


minus_log_REML_p = function(p1=observed$p[1],p2=observed$p[2],p3=observed$p[3],p4=observed$p[4], # there's probably a nicer way to do this...
                            p5=observed$p[5],p6=observed$p[6],p7=observed$p[7],p8=observed$p[8],
                            p9=observed$p[9],p10=observed$p[10]){
  
  a = coef(model_REML)[1]*coef(model_REML)[2]
  b = coef(model_REML)[2]-a
  p =c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
  
  -sum(dbeta(x=p,shape1=a,shape2=b,log=TRUE) + dbinom(x=observed$n_pos,size=observed$n_trials,prob=p,log=TRUE))
  
}
minus_log_REML_p()

model_REML_p = stats4::mle(minus_log_REML_p, start = list(p1=0.5,p2=0.5,p3=0.5,p4=0.5,
                                                        p5=0.5,p6=0.5,p7=0.5,p8=0.5,p9=0.5,p10=0.5),
                         method='L-BFGS-B',
                         upper=c(0.99999,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999,0.99999),
                         lower=c(0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001))

summary(model_REML_p)

#REML estimate works
REML_params=data.frame(true=c(m=m,s=s,observed$p), estimate =c(coef(model_REML),coef(model_REML_p)), se = c(sqrt(diag(vcov(model_REML))),sqrt(diag(vcov(model_REML_p))))) |>
  mutate( z = (true-estimate)^2/se^2,
          lower = estimate - 1.96*se,
          upper = estimate + 1.96*se) |>
  rownames_to_column(var='param') |>
  select(param, everything()) |> 
  mutate(param=factor(param,levels=c('m','s','p1','p2','p3','p4','p5','p6','p7','p8','p9','p10'))) |>
  mutate(level=factor(c('group_mean','group_dispersion',rep('subject',10)))) |>
  mutate(param_index = 1:12)

REML_params


### what about trying to estimate the p's using total likelihood?





## show observed


# this works great!

#+ fig.width=9, fig.height=3
ggplot() +
  geom_point(data=full_lik_params,aes(x=param_index,y=estimate),color='red') +
  geom_segment(data=full_lik_params,aes(x=param_index,y=lower,yend=upper,group=param),color='red') +
  geom_point(data=REML_params,aes(x=param_index+0.1,y=estimate),color='blue') +
  geom_segment(data=REML_params,aes(x=param_index+0.1,y=lower,yend=upper,group=param),color='blue') +
  geom_point(data=full_lik_params,aes(x=param_index-0.1,y=true)) +
  theme_bw()+
  facet_wrap('level',scales='free') +
  scale_x_continuous(breaks=seq(1:12),labels=levels(full_lik_params$param),minor_breaks = NULL) +
  xlab('') +
  geom_text(data=data.frame(param_index=rep(1.945,3),y=full_lik_params$upper[2]*c(0.95,0.85,0.75),label=c('true','full likelihood','REML'),level=rep('group_dispersion',3)),
            aes(x=param_index,y=y,label=label),color=c('black','red','blue'))
# /*
ggsave('scratch/empirical_bayes_random_effects_via_optimization_example_parameters.png',units='in',width=7,height=3)
# */

# spin it to a blog post!
# https://deanattali.com/2015/03/24/knitrs-best-hidden-gem-spin/



# /*
rmarkdown::render(input='./scratch/empirical_bayes_random_effects_via_optimization.R',
                  output_format = c('md_document','html_document'))
# */
               