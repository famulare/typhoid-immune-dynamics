# nlmixr2_scratch
# exploring suitability of existing nonlinear mixed models from PKPD literature 
# as basis for typhoid immunity model

# I suspect the challenge will be I want to have stuff that isn't quite
# simple ODEs, but superpositions of simple odes. But let's see.

# following example here: https://nlmixr2.org/reference/nlmixr2-package/#example

library(nlmixr2)
library(xpose)
library(xpose.nlmixr2)
library(tidyverse)

theme_set(theme_bw())

## The basic model consists of an ini block that has initial estimates
one.compartment <- function() {
  ini({
    tka <- log(1.57); label("Ka")
    tcl <- log(2.72); label("Cl")
    tv <- log(31.5); label("V")
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  # and a model block with the error specification and model specification
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

# data 
data(theo_sd)
theo_sd = theo_sd |>
  mutate(ID=factor(ID))
# View(theo_sd)
names(theo_sd)

p1=ggplot(theo_sd) +
  geom_line(aes(x=TIME,y=DV,group=ID,color=ID)) +
  geom_point(aes(x=TIME,y=DV,group=ID,color=ID)) 
p1
p1+
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')


# fit model
fit <- nlmixr2(one.compartment, theo_sd,  est="saem", saemControl(print=0))

print(fit)

## ggplot fits
p2 = p1 + 
  geom_line(data=data.frame(DV=fit$IPRED,TIME=fit$TIME,ID=fit$ID),
            aes(x=TIME,y=DV,group=ID,color=ID),linetype='dashed') +
  geom_line(data=data.frame(DV=fit$PRED,TIME=fit$TIME,ID=fit$ID),
            aes(x=TIME,y=DV,group=ID,color=ID),linetype='dotted') +
  facet_wrap('ID')
p2

p2 +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')


# do some cool plots with xpose
xpdb = xpose_data_nlmixr(fit)


# https://uupharmacometrics.github.io/xpose/articles/plot_list.html
dv_vs_pred(xpdb)
dv_vs_ipred(xpdb)

dv_vs_idv(xpdb)
pred_vs_idv(xpdb)
ipred_vs_idv(xpdb)

dv_preds_vs_idv(xpdb)
ind_plots(xpdb, page = 1)


prm_vs_iteration(xpdb)

prm_distrib(xpdb)
eta_distrib(xpdb)

res_qq(xpdb)

amt_vs_idv(xpdb)

# these don't work. Some bug in database
get_prm(xpdb)

xpdb %>% 
  vpc_data(vpc_type = 'continuous',opt = vpc_opt(n_bins = 'auto')) %>% 
  vpc()
