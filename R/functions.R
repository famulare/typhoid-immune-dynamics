# --


# Beta-Poisson dose-response function with immune scaling
p_outcome_given_dose <- function(dose, alpha, N50, gamma, CoP) {
    scale_factor <- (2^(1 / alpha) - 1)
    base <- 1 + dose * scale_factor / N50
    exponent <- -alpha / (CoP^gamma)
    p <- 1 - base^exponent
    return(p)
}
