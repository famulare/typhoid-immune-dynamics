# Fit summary: t2-indiv-vax-prior

**Date:** 2026-08-01 01:42

## Sampler diagnostics

| Metric | Value |
|---|---|
| Divergent transitions | 0 / 1000 |
| Divergence rate |    0% |
| Max-treedepth hits | 0 |
| E-BFMI (min) | 0.949 |

## Parameters

| param | mean | median | sd | 5% | 95% | ess_bulk | ess_tail | rhat |
|---|---|---|---|---|---|---|---|---|
| `log10_N50_inf` | 2.305 | 2.293 | 0.7681 | 1.065 | 3.597 | 1e+03 | 7e+02 | 1.011 |
| `d_fev` | 0.7121 | 0.6118 | 0.5119 | 0.08667 | 1.708 | 8e+02 | 5e+02 | 1.001 |
| `alpha_inf` | 0.302 | 0.2128 | 0.3319 | 0.05608 | 0.8276 | 1e+03 | 7e+02 | 1.001 |
| `alpha_fevginf` | 0.3208 | 0.231 | 0.293 | 0.05971 | 0.8717 | 2e+03 | 7e+02 | 1.001 |
| `gamma_inf` | 0.2614 | 0.2123 | 0.2017 | 0.06162 | 0.6174 | 1e+03 | 7e+02 | 0.999 |
| `gamma_fevginf` | 0.2555 | 0.2005 | 0.1891 | 0.06492 | 0.6797 | 2e+03 | 6e+02 | 1.004 |
| `log10_delta` | 3.493 | 3.492 | 0.6929 |  2.34 | 4.584 | 1e+03 | 7e+02 |     1 |
| `pi_susc` | 0.6357 | 0.646 | 0.1401 | 0.3875 | 0.8412 | 1e+03 | 8e+02 | 1.001 |
| `CoP_imm` | 12.86 | 8.628 | 13.22 | 0.7205 |  38.1 | 1e+03 | 6e+02 | 1.002 |
| `CoP_susc` | 1.047 | 1.014 | 0.3181 | 0.6168 | 1.607 | 1e+03 | 8e+02 |     1 |
| `phi0_a` | 1.033 |     1 | 1.562 | -1.593 | 3.646 | 1e+03 | 8e+02 | 1.002 |
| `phi0_b` | 1.536 | 1.246 | 1.282 | 0.06975 |  4.05 | 7e+02 | 4e+02 | 0.9989 |
| `eta_lo` | 0.4963 | 0.4946 | 0.1521 | 0.2421 | 0.754 | 1e+03 | 8e+02 |     1 |
| `kappa` | 1.678 | 0.9911 | 2.447 | 0.1916 | 5.234 | 2e+03 | 6e+02 | 0.9997 |
| `grand_overdispersion_rho` | 0.02027 | 0.01407 | 0.0205 | 0.00113 | 0.06025 | 1e+03 | 4e+02 | 1.011 |
| `log_V_M01ZH09` | 0.01299 | 0.02082 | 0.9725 | -1.624 | 1.629 | 2e+03 | 5e+02 | 1.002 |
| `log_V_Ty21a` | 0.04102 | 0.04314 | 0.9774 | -1.548 | 1.652 | 9e+02 | 8e+02 | 1.003 |
| `psi_stool` | 0.5978 | 0.6064 | 0.1963 | 0.2571 | 0.8837 | 2e+03 | 5e+02 | 1.004 |
| `frac_late` | 0.7271 | 0.7562 | 0.1785 | 0.378 | 0.9632 | 1e+03 | 4e+02 | 1.011 |
| `log10_N50_fevginf` | 3.017 | 3.001 | 0.7725 | 1.774 | 4.349 | 2e+03 | 9e+02 | 1.011 |

## Strong pairwise correlations (|r| > 0.7)

| pair | r |
|---|---|
| `log10_N50_inf` - `log10_N50_fevginf` | 0.779 |

## priorsense power-scaling sensitivity

| variable | prior | likelihood | diagnosis |
|---|---|---|---|
| `log10_N50_inf` | 0.0738 | 8.91 | potential prior-data conflict |
| `d_fev` | 0.216 |  8.8 | potential prior-data conflict |
| `alpha_inf` |  0.6 | 9.73 | potential prior-data conflict |
| `alpha_fevginf` | 0.59 | 8.57 | potential prior-data conflict |
| `gamma_inf` | 0.483 | 12.7 | potential prior-data conflict |
| `gamma_fevginf` | 0.405 | 5.51 | potential prior-data conflict |
| `log10_delta` | 0.106 | 13.7 | potential prior-data conflict |
| `pi_susc` | 0.0881 | 6.22 | potential prior-data conflict |
| `CoP_imm` | 0.473 | 8.87 | potential prior-data conflict |
| `CoP_susc` |  0.2 | 5.75 | potential prior-data conflict |
| `phi0_a` | 0.111 | 8.27 | potential prior-data conflict |
| `phi0_b` | 0.322 | 10.8 | potential prior-data conflict |
| `eta_lo` | 0.127 | 4.52 | potential prior-data conflict |
| `kappa` | 0.666 | 16.1 | potential prior-data conflict |
| `grand_overdispersion_rho` | 0.529 | 33.3 | potential prior-data conflict |
| `log_V_M01ZH09` | 0.0862 |  5.3 | potential prior-data conflict |
| `log_V_Ty21a` | 0.105 |  7.9 | potential prior-data conflict |
| `psi_stool` | 0.118 |  6.7 | potential prior-data conflict |
| `frac_late` | 0.104 | 5.64 | potential prior-data conflict |
| `log10_N50_fevginf` | 0.0997 | 11.5 | potential prior-data conflict |

## Figures

![calibration_targets](calibration_targets.png)

![cop_response_milk_doses](cop_response_milk_doses.png)

![delta_bridge](delta_bridge.png)

![density](density.png)

![dose_cop_surface](dose_cop_surface.png)

![dose_response_fit](dose_response_fit.png)

![dose_response_grid_all](dose_response_grid_all.png)

![dose_response_grid_maryland](dose_response_grid_maryland.png)

![dose_response_grid_oxford](dose_response_grid_oxford.png)

![energy](energy.png)

![eta_detection](eta_detection.png)

![maryland_mixture](maryland_mixture.png)

![neff](neff.png)

![pairs](pairs.png)

![phi_severity](phi_severity.png)

![ppc](ppc.png)

![prior_posterior](prior_posterior.png)

![rank](rank.png)

![rhat](rhat.png)

![titre_protection](titre_protection.png)

![trace](trace.png)


### grid/

![md_gilman_imm](grid/md_gilman_imm.png)

![md_gilman_mix](grid/md_gilman_mix.png)

![md_gilman_susc](grid/md_gilman_susc.png)

![md_hornick](grid/md_hornick.png)

![md_levine](grid/md_levine.png)

![ox_darton_m01](grid/ox_darton_m01.png)

![ox_darton_plac](grid/ox_darton_plac.png)

![ox_darton_ty21a](grid/ox_darton_ty21a.png)

![ox_jin_ctrl](grid/ox_jin_ctrl.png)

![ox_jin_vips](grid/ox_jin_vips.png)

![ox_jin_vitt](grid/ox_jin_vitt.png)

![ox_naive](grid/ox_naive.png)

