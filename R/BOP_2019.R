
# Formulas under (3)
alpha_star_n_BOP19 <- function(y_n_aver, Sigma_n_inv, mu_n, mu_0) {

  numerator <- t(y_n_aver) %*% Sigma_n_inv %*% mu_n %*% 
               t(mu_0) %*% Sigma_n_inv %*% mu_0-
               t(mu_n) %*% Sigma_n_inv %*% mu_0 %*% 
               t(y_n_aver) %*% Sigma_n_inv %*% mu_0

  denomen   <- t(y_n_aver) %*% Sigma_n_inv %*% y_n_aver %*% 
               t(mu_0) %*% Sigma_n_inv %*% mu_0-
               (t(y_n_aver) %*% Sigma_n_inv %*% mu_0)^2

  as.numeric(numerator / denomen)
}


beta_star_n_BOP19 <- function(y_n_aver, Sigma_n_inv, mu_n, mu_0) {

  numerator <- t(y_n_aver) %*% Sigma_n_inv %*% y_n_aver %*% 
               t(mu_n) %*% Sigma_n_inv %*% mu_0-
               t(y_n_aver) %*% Sigma_n_inv %*% mu_0 %*% 
               t(y_n_aver) %*% Sigma_n_inv %*% mu_n

  denomen   <- t(y_n_aver) %*% Sigma_n_inv %*% y_n_aver %*% 
               t(mu_0) %*% Sigma_n_inv %*% mu_0-
               (t(y_n_aver) %*% Sigma_n_inv %*% mu_0)^2

  as.numeric(numerator / denomen)
}


# From Theorem 1
alpha_star_BOP19 <- function(c, mu_n, Sigma_n_inv, mu_0) {

  I1 <- t(mu_n) %*% Sigma_n_inv %*% mu_n
  I2 <- t(mu_0) %*% Sigma_n_inv %*% mu_0
  I3 <- t(mu_n) %*% Sigma_n_inv %*% mu_0

  numerator <- I1*I2 - I3^2

  I4 <- t(mu_n) %*% Sigma_n_inv %*% mu_n
  I5 <- t(mu_0) %*% Sigma_n_inv %*% mu_0
  I6 <- t(mu_n) %*% Sigma_n_inv %*% mu_0

  denomen <- (c+I4)*I5 - I6^2

  as.numeric(numerator / denomen)
}


beta_star_BOP19 <- function(alpha_star, mu_n_t, Sigma_n_inv, mu_0) {

  I1 <- mu_n_t %*% Sigma_n_inv %*% mu_0
  I2 <- t(mu_0) %*% Sigma_n_inv %*% mu_0

  II <-I1 / I2
  as.numeric((1-alpha_star)*II)
}


# From Theorem 3
alpha_star_hat_BOP19 <- function(n, p, y_n_aver, Sigma_n_inv, mu_0) {

  I1 <- t(y_n_aver) %*% Sigma_n_inv %*% y_n_aver - p/(n-p)
  I2 <- t(mu_0) %*% Sigma_n_inv %*% mu_0
  I3 <- t(y_n_aver) %*% Sigma_n_inv %*% mu_0

  numerator <- I1*I2 - I3^2

  I4 <- t(y_n_aver) %*% Sigma_n_inv %*% y_n_aver
  I5 <- t(mu_0) %*% Sigma_n_inv %*% mu_0
  I6 <- t(y_n_aver) %*% Sigma_n_inv %*% mu_0

  denomen <- I4*I5 - I6^2

  as.numeric(numerator / denomen)
}


beta_star_hat_BOP19 <- function(n, p, alpha_star_hat, 
                                y_n_aver_t, Sigma_n_inv, mu_0) {

  I1 <- y_n_aver_t %*% Sigma_n_inv %*% mu_0
  I2 <- t(mu_0) %*% Sigma_n_inv %*% mu_0

  II <-I1 / I2
  as.numeric((1-alpha_star_hat)*II)
}


# From Theorem 4

s_BOP19 <- function(mu_0, Sigma_n_inv, mu_n){

  I1 <- t(mu_n) %*% Sigma_n_inv %*% mu_n
  I2 <- (t(mu_0) %*% Sigma_n_inv %*% mu_n)^2
  I3 <- t(mu_0) %*% Sigma_n_inv %*% mu_0

  as.numeric(I1 - I2/I3)
}

R_BOP19 <- function(mu_0, Sigma_n_inv, mu_n){

  I1 <- t(mu_0) %*% Sigma_n_inv %*% mu_n
  I2 <- t(mu_0) %*% Sigma_n_inv %*% mu_0
  as.numeric(I1 / I2)
}

sigma_s_square_BOP19 <- function(c, s) 2*(c+2*s) + 2/(1-c)*(c+s)^2


Omega_BOP19 <- function(c, sigma_s_square, Sigma_n_inv, mu_0, s, R){

  El11 <- c^2 * sigma_s_square / (c + s)^4
  El12 <- El21 <- c^2 * sigma_s_square * R / (c + s)^4
  El22 <- c^2 * sigma_s_square * R^2 / (c + s)^4 +
          c^2 * (c + s)^(-2) * (1 + (c + s)/(1-c)) / 
          as.numeric(t(mu_0) %*% Sigma_n_inv %*% mu_0)
  
  matrix(data = c(El11, El12, El21, El22), 2, 2)
}


