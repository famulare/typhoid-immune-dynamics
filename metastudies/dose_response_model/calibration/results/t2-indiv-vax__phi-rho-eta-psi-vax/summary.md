# Fit summary: t2-indiv-vax

**Date:** 2026-08-01 01:43

## Sampler diagnostics

| Metric | Value |
|---|---|
| Divergent transitions | 1 / 4000 |
| Divergence rate | 0.025% |
| Max-treedepth hits | 0 |
| E-BFMI (min) | 0.947 |

## Parameters

| param | mean | median | sd | 5% | 95% | ess_bulk | ess_tail | rhat |
|---|---|---|---|---|---|---|---|---|
| `log10_N50_inf` | 1.565 | 1.601 | 0.4609 | 0.7653 | 2.264 | 2e+03 | 2e+03 | 1.001 |
| `d_fev` | 0.5549 | 0.4889 | 0.3935 | 0.05245 | 1.301 | 2e+03 | 2e+03 | 1.001 |
| `alpha_inf` | 0.2187 | 0.2122 | 0.05042 | 0.1478 | 0.3076 | 3e+03 | 2e+03 |     1 |
| `alpha_fevginf` | 0.2767 | 0.2668 | 0.06834 | 0.1822 | 0.4026 | 4e+03 | 3e+03 | 1.001 |
| `gamma_inf` | 0.08627 | 0.08118 | 0.03579 | 0.03686 | 0.1508 | 4e+03 | 3e+03 | 0.9999 |
| `gamma_fevginf` | 0.1957 | 0.1921 | 0.06603 | 0.08909 | 0.3075 | 4e+03 | 2e+03 | 1.003 |
| `log10_delta` |  2.02 | 2.005 | 0.3366 | 1.493 | 2.585 | 3e+03 | 2e+03 | 1.001 |
| `pi_susc` | 0.657 | 0.6705 | 0.1341 | 0.4216 | 0.8583 | 5e+03 | 3e+03 |     1 |
| `CoP_imm` | 13.04 | 9.226 | 12.78 | 1.035 | 37.36 | 3e+03 | 2e+03 | 1.001 |
| `CoP_susc` | 0.9531 | 0.9134 | 0.2865 | 0.5619 | 1.478 | 5e+03 | 3e+03 |     1 |
| `phi0_a` | 1.255 | 1.232 | 0.4093 | 0.6119 | 1.947 | 3e+03 | 3e+03 | 1.001 |
| `phi0_b` | 1.834 | 1.814 | 0.5894 | 0.9054 | 2.844 | 3e+03 | 2e+03 | 1.001 |
| `eta_lo` | 0.7108 | 0.7108 | 0.05783 | 0.6141 | 0.8062 | 4e+03 | 2e+03 | 1.002 |
| `kappa` | 1.619 | 0.9779 | 2.114 | 0.1725 |  5.09 | 6e+03 | 3e+03 | 1.002 |
| `grand_overdispersion_rho` | 0.02655 | 0.02435 | 0.01564 | 0.005141 | 0.05468 | 3e+03 | 1e+03 | 1.001 |
| `log_V_M01ZH09` | 0.3635 | 0.3609 | 0.1868 | 0.05962 | 0.6718 | 4e+03 | 3e+03 |     1 |
| `log_V_Ty21a` | 0.6368 | 0.6328 | 0.209 | 0.2999 | 0.9864 | 4e+03 | 3e+03 |     1 |
| `psi_stool` | 0.8141 | 0.8156 | 0.06497 | 0.7061 | 0.9169 | 4e+03 | 2e+03 | 1.001 |
| `frac_late` | 0.8589 | 0.8773 | 0.09512 | 0.6811 | 0.9797 | 5e+03 | 2e+03 | 1.003 |
| `log10_N50_fevginf` |  2.12 | 2.142 | 0.3882 | 1.454 | 2.715 | 3e+03 | 3e+03 | 1.002 |
| `N50_inf` | 60.05 | 39.93 | 65.66 | 5.826 | 183.6 | 2e+03 | 2e+03 | 1.001 |
| `N50_fevginf` | 189.2 | 138.6 | 171.2 | 28.47 |   519 | 3e+03 | 3e+03 | 1.003 |
| `delta` | 144.4 |   101 | 150.2 | 31.09 | 384.5 | 3e+03 | 2e+03 | 1.001 |
| `grand_concentration_k` | 85.86 | 40.07 | 501.3 | 17.29 | 193.5 | 3e+03 | 1e+03 | 1.001 |
| `V_M01ZH09` | 1.464 | 1.435 | 0.2779 | 1.061 | 1.958 | 4e+03 | 3e+03 |     1 |
| `V_Ty21a` | 1.932 | 1.883 | 0.4111 |  1.35 | 2.681 | 4e+03 | 3e+03 |     1 |

## Strong pairwise correlations (|r| > 0.7)

| pair | r |
|---|---|
| `log_V_M01ZH09` - `V_M01ZH09` | 0.991 |
| `log_V_Ty21a` - `V_Ty21a` | 0.989 |
| `phi0_a` - `phi0_b` | 0.766 |

## priorsense power-scaling sensitivity

| variable | prior | likelihood | diagnosis |
|---|---|---|---|
| `log10_N50_inf` | 0.0608 | 0.105 | potential prior-data conflict |
| `d_fev` | 0.134 | 0.0947 | potential prior-data conflict |
| `alpha_inf` | 0.0327 | 0.102 | - |
| `alpha_fevginf` | 0.155 | 0.0568 | potential prior-data conflict |
| `gamma_inf` | 0.0741 | 0.199 | potential prior-data conflict |
| `gamma_fevginf` | 0.0787 | 0.0803 | potential prior-data conflict |
| `log10_delta` | 0.149 | 0.132 | potential prior-data conflict |
| `pi_susc` | 0.0976 | 0.0144 | potential strong prior / weak likelihood |
| `CoP_imm` | 0.425 | 0.082 | potential prior-data conflict |
| `CoP_susc` | 0.159 | 0.0354 | potential strong prior / weak likelihood |
| `phi0_a` | 0.0531 | 0.0664 | potential prior-data conflict |
| `phi0_b` | 0.0751 | 0.0707 | potential prior-data conflict |
| `eta_lo` | 0.0696 | 0.127 | potential prior-data conflict |
| `kappa` | 0.631 | 0.0498 | potential strong prior / weak likelihood |
| `grand_overdispersion_rho` | 0.172 | 0.143 | potential prior-data conflict |
| `log_V_M01ZH09` | 0.0124 | 0.0968 | - |
| `log_V_Ty21a` | 0.014 | 0.0802 | - |
| `psi_stool` | 0.03 | 0.0875 | - |
| `frac_late` | 0.0765 | 0.202 | potential prior-data conflict |
| `log10_N50_fevginf` | 0.0879 | 0.0868 | potential prior-data conflict |
| `N50_inf` | 0.218 | 0.107 | potential prior-data conflict |
| `N50_fevginf` | 0.255 | 0.0803 | potential prior-data conflict |
| `delta` | 0.179 | 0.318 | potential prior-data conflict |
| `grand_concentration_k` | 0.253 | 0.123 | potential prior-data conflict |
| `V_M01ZH09` | 0.0124 | 0.131 | - |
| `V_Ty21a` | 0.0194 | 0.111 | - |

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

