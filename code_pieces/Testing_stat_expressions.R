
n <-1e3
x <-matrix(data = rnorm(n*p, mean=mu), nrow = p, ncol = n)


abs(V_hat_GMV(x) - V_GMV(Sigma))

r<-replicate(n=1000, {
  x <-matrix(data = rnorm(n*p, mean=mu), nrow = p, ncol = n)
  abs(R_GMV(mu, Sigma) - R_hat_GMV(x))
})

max(r)
