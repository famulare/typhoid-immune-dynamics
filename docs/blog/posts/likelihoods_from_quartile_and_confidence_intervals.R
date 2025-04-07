#' ---
#' title: "Likelihoods from quartiles and confidence intervals, or how I learned to like order statistics"
#' output:
#'   md_document
#' knit: (function(input, ...) {
#'   out <- rmarkdown::render(input, output_dir = "./docs/blog/posts", ...)
#'   return(out)})
#' ---
#'
#' 
#' # Likelihoods from quartile and confidence intervals, or how I learned to like order statistics
#' 
#' #' This is a quick experiment to learn how to formalize a likelihood function from reported quartiles
#' and confidence intervals, using order statistics as explained in this [Cross-Validated answer](https://stats.stackexchange.com/a/390981/252994) 
#' by the prolific [Kjetil B Halvorsen](https://stats.stackexchange.com/users/11887/kjetil-b-halvorsen).
#' 
#' **Why?** Because quartiles and confidence intervals are the standard outputs from studies, and it's often much more difficult
#' or impossible to reconstruct the individual-level study design and analysis to reuse the data for my purposes. Since model calibration
#' is still not a standardized thing, I can't just pull metastudy regression tools off the shelf and need my own thing. 
#' 
#' <!-- more -->
#' 
#' **What?** Construct likelihood functions from quartiles or confidence intervals, assuming log-relative-risks for vaccine efficacy 
#' or lognormals for antibody concentration measures. Questions are how to deal with sample sizes and censoring. 
#' 
#' **What did I learn?** 
#'
#' **Note about how this post was generated.** This post was generated directly from the commented R script using the `knitr::spin` functionality
#' of `rmarkdown::render`. This gives all the advantages of scripts and r-markdown, without the disadvantages. 
#' Check out this [blog post by Dean Attali](https://deanattali.com/2015/03/24/knitrs-best-hidden-gem-spin/) to learn more.
#' 
#' ## Setting up the study
#' 
#' First, 
#'
#+ echo=TRUE, message=FALSE, results = 'hide'
# set up the environment
library(tidyverse)
library(stats4)
library(knitr)

#' first experiment: an easy example with no censoring

# assumed truth
  log_mu = 6.9
  log_sd = 0.5
# not directly reported
  log_titer_summary = data.frame(gmt = log_mu, lower = log_mu-1.96*log_sd, upper = log_mu+1.96*log_sd) 

# reported (exact, for now)
N = 100
titer_summary = round(exp(log_titer_summary)) # reported transformation
titer_summary

# likelihood with CDFs only
titer_order_stats_likelihood = function(mu=6.9,sigma=0.5){
  mlogL = -N*(
              0.025*log(pnorm((log(titer_summary$lower)-mu)/sigma)) +
              0.475*log(pnorm((log(titer_summary$gmt)-mu)/sigma)-pnorm((log(titer_summary$lower)-mu)/sigma)) +
              0.475*log(pnorm((log(titer_summary$upper)-mu)/sigma)-pnorm((log(titer_summary$gmt)-mu)/sigma)) +
              0.025*log(1-pnorm((log(titer_summary$upper)-mu)/sigma))
            ) -
          # log(dnorm((log(titer_summary$lower)-mu)/sigma)) -
          # log(dnorm((log(titer_summary$gmt)-mu)/sigma)) -
          # log(dnorm((log(titer_summary$upper)-mu)/sigma)) - 
          # 3*log(sigma) -
          0
    
  names(mlogL)=NULL
  return(mlogL)
}

# find MLE and show that it gets back what was assumed at the top
mod = mle(titer_order_stats_likelihood,
          start=list(mu=6,sigma=0.5),
          method= 'L-BFGS-B',
          lower=list(mu=-Inf,sigma=0.01),
          upper=list(mu=Inf,sigma=3))
summary(mod) 

data.frame(truth=c(log_mu,log_sd),mle=round(coef(mod),4),
           lower=round(coef(mod)-1.96*sqrt(diag(vcov(mod))),4),
           upper=round(coef(mod)+1.96*sqrt(diag(vcov(mod))),4)) |> kable()




#' Next, let's look at an example with censoring at the lower level...
#' 

# assumed truth
log_mu = 1
log_sd = 0.5
# not directly reported
log_titer_summary = data.frame(gmt = log_mu, lower = log_mu-1.96*log_sd, upper = log_mu+1.96*log_sd) 

# reported (exact, for now)
N = 1000
titer_summary = (exp(log_titer_summary)) # reported transformation
titer_summary

titer_summary$lower=2

# censored likelihood
titer_order_stats_likelihood = function(mu=1,sigma=0.5){

  l=1
  u=2
  integrand = function(x){pnorm((log(x)-mu)/sigma)/(log(u)-log(l))}  
# 
#   mlogL = -N*(
#       0.025*log(integrate(integrand,lower=l,upper=u)$value) +
#       0.475*log(pnorm((log(titer_summary$gmt)-mu)/sigma)-integrate(integrand,lower=l,upper=u)$value) +
#       0.475*log(pnorm((log(titer_summary$upper)-mu)/sigma)-pnorm((log(titer_summary$gmt)-mu)/sigma)) +
#       0.025*log(1-pnorm((log(titer_summary$upper)-mu)/sigma)))
  
  mlogL = -N*(
      0.5*log(pnorm((log(titer_summary$gmt)-mu)/sigma)) +
      0.475*log(pnorm((log(titer_summary$upper)-mu)/sigma)-pnorm((log(titer_summary$gmt)-mu)/sigma)) +
      0.025*log(1-pnorm((log(titer_summary$upper)-mu)/sigma))
  )
  names(mlogL)=NULL
  return(mlogL)
}
titer_order_stats_likelihood()

# find MLE and show that it gets back what was assumed at the top
mod = mle(titer_order_stats_likelihood,
          start=list(mu=1,sigma=0.5))
summary(mod) 

data.frame(truth=c(log_mu,log_sd),mle=round(coef(mod),4),
           lower=round(coef(mod)-1.96*sqrt(diag(vcov(mod))),4),
           upper=round(coef(mod)+1.96*sqrt(diag(vcov(mod))),4)) |> kable()



#' Next, let's look at an example with censoring thru mean...
#' 

# assumed truth
log_mu = 1
log_sd = 0.5
# not directly reported
log_titer_summary = data.frame(gmt = log_mu, lower = log_mu-1.96*log_sd, upper = log_mu+1.96*log_sd) 

# reported (exact, for now)
N = 1000
titer_summary = round(exp(log_titer_summary)) # reported transformation
titer_summary

titer_summary$lower=5
titer_summary$gmt=5

# uncensored likelihood
titer_order_stats_likelihood = function(mu=1,sigma=0.5){
  mlogL = -N*(
    0.975*log(pnorm((log(titer_summary$upper)-mu)/sigma)) +
      0.025*log(1-pnorm((log(titer_summary$upper)-mu)/sigma)) 
  )
  names(mlogL)=NULL
  return(mlogL)
}

# find MLE and show that it gets back what was assumed at the top
mod = mle(titer_order_stats_likelihood,
          start=list(mu=1,sigma=0.5))
summary(mod) 

data.frame(truth=c(log_mu,log_sd),mle=round(coef(mod),4),
           lower=round(coef(mod)-1.96*sqrt(diag(vcov(mod))),4),
           upper=round(coef(mod)+1.96*sqrt(diag(vcov(mod))),4)) |> kable()



#' so let's make a function that builds likelihoods from a table
data = data.frame(interval = c('95CI','95CI','quartile','quartile'),
                middle = c(100,110,90,80),
                lower=c(70,72,50,45),
                upper = c(135,150,200,290),
                censoring = c('none','none','none','none'),
                N= c(100,50,200,150))
data = data[3,]

summary_statistics_likelihood = function(mu,sigma){

  
  # initialize minus log likelihood
  mlogL = 0
  
  for(row in 1:nrow(data)){
    if (data$interval[row] == '95CI'){
      weights = c(0.025,0.475,0.475,0.025)
    } else if (data$interval[row] == 'quartile'){
      weights = c(0.25,0.25,0.25,0.25)
    } else {
      errorCondition('invalid interval type')
    }
    
    if(data$censoring[row] == 'none'){
      
      mlogL = mlogL - 
        data$N[row]*(
          weights[1]*log(pnorm((log(data$lower[row])-mu)/sigma)) +
          weights[2]*log(pnorm((log(data$middle[row])-mu)/sigma)-pnorm((log(data$lower[row])-mu)/sigma)) +
          weights[3]*log(pnorm((log(data$upper[row])-mu)/sigma)-pnorm((log(data$middle[row])-mu)/sigma)) +
          weights[4]*log(1-pnorm((log(data$upper[row])-mu)/sigma))
        )
    
    } else if (data$censoring[row] == 'lower'){
      
      mlogL = mlogL - 
        data$N[row]*(
          (weights[1]+weights[2])*log(pnorm((log(data$middle[row])-mu)/sigma)) +
            weights[3]*log(pnorm((log(data$upper[row])-mu)/sigma)-pnorm((log(data$middle[row])-mu)/sigma)) +
            weights[4]*log(1-pnorm((log(data$upper[row])-mu)/sigma))
        )
    } else if (data$censoring[row] == 'middle'){
      
      mlogL = mlogL - 
        data$N[row]*(
          (weights[1]+weights[2]+weights[3])*log(pnorm((log(data$upper[row])-mu)/sigma)) +
            weights[4]*log(1-pnorm((log(data$upper[row])-mu)/sigma))
        )
      
    } else if (data$censoring[row] == 'upper'){
      
      mlogL = mlogL - 
        data$N[row]*(
          (weights[1]+weights[2]+weights[3] + weights[4])*log(pnorm((log(data$upper[row])-mu)/sigma)) 
        )
      
    }
  }
  return(mlogL)
}

summary_statistics_likelihood(mu=2,sigma=0.5)

mod = mle(summary_statistics_likelihood,
          start=list(mu=4,sigma=0.5))
summary(mod) 

coef(mod)
c(exp(coef(mod)[1]),exp(coef(mod)[1]-2*coef(mod)[2]),exp(coef(mod)[1]+2*coef(mod)[2]),exp(qnorm(0.25,coef(mod)[1],coef(mod)[2])),exp(qnorm(0.75,coef(mod)[1],coef(mod)[2])))
data


# /* back matter for exporting as a blog post

source('./docs/docs_helper_functions.R')
render_blog_post(input_file = './docs/blog/posts/likelihoods_from_quartile_and_confidence_intervals.R',
                 categories_list = c('Calibration','Self-study'),
                 date_created = '2025-04-07')

# */
