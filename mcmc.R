# Example with a Student t ================================================
# Constructing synthetic mcmc output
mu <- c(0.5, 6)
mu_mat <- matrix(rep(mu, 100) + rnorm(200, 0, 0.1),
                 ncol = 2, byrow = TRUE
)
omega <- c(1, 2)
sigma_mat <- matrix(rep(omega, 100) + rnorm(200, 0, 0.1),
                    ncol = 2, byrow = TRUE
)
nu <- c(5, 5)
nu_mat <- matrix(rep(nu, 100) + rnorm(200, 0, 0.1),
                 ncol = 2, byrow = TRUE
)
eta <- c(0.8, 0.2)
eta_mat <- matrix(rep(eta[1], 100) + rnorm(100, 0, 0.05),
                  ncol = 1
)
eta_mat <- cbind(eta_mat, 1 - eta_mat)
xi_mat <- matrix(0, 100, 2)
fit <- cbind(eta_mat, mu_mat, sigma_mat, nu_mat, xi_mat)
colnames(fit) <- c(
  "eta1", "eta2", "mu1", "mu2",
  "omega1", "omega2", "nu1", "nu2", "xi1", "xi2"
)
# sampling observations
data <- c(
  sn::rst(eta[1] * 1000, mu[1], omega[1], nu = nu[1]),
  sn::rst(eta[2] * 1000, mu[2], omega[2], nu = nu[2])
)

pdf_func <- function(x, pars) {
    sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
}

dist_type <- "continuous"

BM <- bayes_mixture(fit, data,
                    burnin = 50,
                    pdf_func = pdf_func, dist_type = dist_type,
                    vars_to_rename = c("sigma" = "omega"), loc = "xi"
)

plot(BM)

