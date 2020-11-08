



nonlin_shrinkLW = function(X){ # 
  n = nrow(X)
  p = ncol(X)
  sampleC = t(X) %*% X / n  # Doesn't match to the def (no mean subtraction)
  eig = eigen(sampleC)
  u = eig$vectors[,p:1]
  lambda = rev(eig$values)
  lambda = lambda[max(1, p-n+1):p]
  L = matrix(rep(lambda, min(p, n)), nrow = length(lambda))
  h = n^(-1/3)
  H = h * t(L)
  x = (L - t(L)) / H
  ftilde = (3/4/sqrt(5)) * rowMeans(pmax(1-x^2/5, 0) / H)
  Hftemp = (-3/10/pi) * x + (3/4/sqrt(5)/pi) * (1 - x^2./5) * log(abs((sqrt(5) - x)/(sqrt(5) + x)))
  Hftemp[abs(x) == sqrt(5)] = (-3/10/pi) * x[abs(x) == sqrt(5)]
  Hftilde = rowMeans(Hftemp / H)
  
  if(p<=n){
    dtilde = lambda / ((pi*(p/n)*lambda*ftilde)^2 + (1-(p/n)-pi*(p/n)*lambda*Hftilde)^2);
  }else{
    Hftilde0 = (1/pi)*(3/10/h^2+3/4/sqrt(5)/h*(1-1/5/h^2) *log((1+sqrt(5)*h)/(1-sqrt(5)*h))) * mean(1/lambda)
    dtilde0=1/(pi*(p-n)/n*Hftilde0) 
    dtilde1 = lambda/(pi^2*lambda^2 * (ftilde^2+Hftilde^2))
    dtilde = c(dtilde0 * rep(1, p-n), dtilde1)
  }
  
  # Hftilde0 = (1/pi)*(3/10/h^2+3/4/sqrt(5)/h*(1-1/5/h^2) *log((1+sqrt(5)*h)/(1-sqrt(5)*h))) * mean(1/lambda)
  # dtilde0=(pi*Hftilde0) 
  #  dtilde0
  u %*% diag(dtilde) %*% t(u)
} # analytical nonlinear shrinkage