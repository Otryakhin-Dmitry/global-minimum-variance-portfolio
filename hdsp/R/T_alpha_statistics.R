
d_0 <- function(gamma, p, n){

  c((1+1/(1-p/n))/gamma,
    -1,
    1/(gamma^2*(1-p/n)),
    (-1-1/(1-p/n))/gamma,
    1)
}


#
# t(d_0) * Omega_hat_al_c * d_0

# This function contains both c and c_n. Is that right?
Omega_hat_al_c <- function(x, c, b){

  c_n <- nrow(x)/ncol(x)
  M <- matrix(data=rep(0,25), nrow=5, ncol=5)
  s_hat_c <- s_hat_c(x)
  V_hat_c <- V_hat_GMV(x)/(1-c_n)
  V_hat_b <- V_hat_b(x, b)

  diag(M) <- c(V_hat_c*(s_hat_c+1)/(1-c), (2*V_hat_c^2)/(1-c),
               2*((s_hat_c+1)^2+c-1)/(1-c), V_hat_b, 2*V_hat_b^2)

  R_hat_b <- R_hat_b(x=x, b=b)
  R_hat_GMV<-R_hat_GMV(x=x)

  M[4,1] <- M[1,4] <- V_hat_c
  M[5,1] <- M[1,5] <- -2*V_hat_c*(R_hat_b-R_hat_GMV)
  M[5,2] <- M[2,5] <- 2*V_hat_c^2
  M[4,3] <- M[3,4] <- 2*(R_hat_b-R_hat_GMV)
  M[5,3] <- M[3,5] <- -2*(R_hat_b-R_hat_GMV)^2
  M
}


# T_alpha, formula (44)
#' @export
T_alpha <- function(gamma, x, w_0, beta=0.05) {

  n <- ncol(x)
  p <- nrow(x)

  Omega_hat_al_c <- Omega_hat_al_c(x=x, c=p/n, b=w_0)
  d_0<-d_0(gamma, p, n)
  B_hat <- B_hat(gamma=gamma, x=x, b=w_0)

  T_alpha <- as.numeric(sqrt(n) * alpha_hat_star_c(gamma=gamma, x=x, b=w_0) *
                        B_hat / sqrt(t(d_0) %*% Omega_hat_al_c %*% d_0))

  p_value <- 2*(1-pnorm(abs(T_alpha), mean = 0, sd = 1))
  z <- qnorm(p=1-beta/2 , mean = 0, sd = 1)
  alpha_lower <- as.numeric(z/sqrt(n) * sqrt(t(d_0) %*% Omega_hat_al_c %*% d_0) / B_hat)
  alpha_higher <- -alpha_lower
  list(T_alpha=T_alpha, p_value=p_value,
       alpha_lower=alpha_lower, alpha_higher=alpha_higher)
}


