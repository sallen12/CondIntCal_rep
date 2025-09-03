# wrapper to get prediction intervals from the 6 forecast distributions in the simulation study
get_ints <- function(mu, tau, alpha) {
  n <- length(mu)

  # Ideal forecaster: F = N(mu, 1)
  L_id <- qnorm(alpha/2, mu)
  U_id <- qnorm(1 - alpha/2, mu)
  int_id <- data.frame(Lower = L_id, Upper = U_id)


  # Unconditional forecaster: F = N(0, 2)
  L_cl <- qnorm(alpha/2, 0, sqrt(2)) |> rep(n)
  U_cl <- qnorm(1 - alpha/2, 0, sqrt(2)) |> rep(n)
  int_cl <- data.frame(Lower = L_cl, Upper = U_cl)

  # Unfocused forecaster: F = 0.5 N(mu, 1) + 0.5 N(mu + tau, 1), tau = -1/+1 w.e.p
  F_un <- function(q, m, t) 0.5*pnorm(q, m) + 0.5*pnorm(q, m + t)
  Q_un <- function(p, m, t) uniroot(function(x) F_un(x, m, t) - p, interval = c(-10, 10))$root
  L_un <- sapply(1:n, function(i) Q_un(alpha/2, mu[i], tau[i]))
  U_un <- sapply(1:n, function(i) Q_un(1 - alpha/2, mu[i], tau[i]))
  int_un <- data.frame(Lower = L_un, Upper = U_un)

  # mean biased forecaster: F = N(mu + tau, 1), tau = -1/+1 w.e.p
  F_me <- function(q, m, t) pnorm(q, m + t)
  L_me <- qnorm(alpha/2, mu + tau)
  U_me <- qnorm(1 - alpha/2, mu + tau)
  int_me <- data.frame(Lower = L_me, Upper = U_me)

  # sign biased forecaster: F = N(-mu, 1)
  L_si <- qnorm(alpha/2, mean = -mu, sd = 1)
  U_si <- qnorm(1 - alpha/2, mean = -mu, sd = 1)
  int_si <- data.frame(Lower = L_si, Upper = U_si)

  # mixed (unconditional/signed) forecaster: F = N(0, 2)/N(-mu,1) w.e.p
  L_1 <- qnorm(alpha/2, 0, sqrt(2)) |> rep(n)
  U_1 <- qnorm(1 - alpha/2, 0, sqrt(2)) |> rep(n)
  L_2 <- qnorm(alpha/2, mean = -mu, sd = 1)
  U_2 <- qnorm(1 - alpha/2, mean = -mu, sd = 1)
  r <- sample(c(0, 1), n, replace = T)
  L_mi <- r*L_1 + (1 - r)*L_2
  U_mi <- r*U_1 + (1 - r)*U_2
  int_mi <- data.frame(Lower = L_mi, Upper = U_mi)

  return(list("Ideal" = int_id, "Uncon." = int_cl, "Unfoc." = int_un, "MB" = int_me, "SB" = int_si, "Mixed" = int_mi))
}
