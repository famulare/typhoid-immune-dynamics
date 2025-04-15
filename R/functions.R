# --
library(assertthat)

# Beta-Poisson function (no immune scaling)
beta_poisson <- function(dose, alpha, N50) {
    1 - (1 + dose * (2^(1 / alpha) - 1) / N50)^(-alpha)
}

# Beta-Poisson dose-response function with immune scaling
beta_poisson_CoP <- function(dose, alpha, N50, gamma, CoP) {
    1 - (1 + dose * (2^(1 / alpha) - 1) / N50)^(-alpha / CoP^gamma)
}

p_outcome_given_dose <- function(dose = 1e4, CoP = 1, outcome = "fever | dose", param_list) {
    if (outcome == "fever | inf") {
        # fever | info is calculated from probability of each outcome
        p_fever <- p_outcome_given_dose(dose = dose, CoP = CoP, outcome = "fever | dose", param_list = param_list)
        p_inf <- p_outcome_given_dose(dose = dose, CoP = CoP, outcome = "inf | dose", param_list = param_list)
        return(p_fever / p_inf)
    }

    assert_that(outcome %in% names(param_list), msg = paste("Missing outcome:", outcome))

    pars <- param_list[[outcome]]

    return(beta_poisson_CoP(
        dose  = dose,
        alpha = pars$alpha,
        N50   = pars$N50,
        gamma = pars$gamma,
        CoP   = CoP
    ))
}

make_CoP_labels <- function(x) {
    sprintf("CoP = %d", x)
}
