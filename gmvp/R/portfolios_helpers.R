
# James-Stein mus and alphas

mu_alp_hat_JS <- function(x, mu_0=1){

  Sigma_hat <- Sigma_sample_estimator(x)
  mu_hat <- rowMeans(x, na.rm = TRUE)
  p <- length(mu_hat)
  I_vect <- rep(1, times=p)

  alp_JS_hat <- as.numeric((p+2) / (p+2 + n*t(mu_hat-mu_0*I_vect)%*%Sigma_hat%*%(mu_hat-mu_0*I_vect)))
  mu_hat_JS <- (1-alp_JS_hat) * mu_hat + alp_JS_hat * mu_0 * I_vect
  list(alp_JS_hat=alp_JS_hat, mu_hat_JS=mu_hat_JS)
}













