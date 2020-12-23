cov.m = function(X, e.mean){n = dim(X)[1]; t(X) %*% X / (n-1) - e.mean %*% t(e.mean) * n / (n-1)}
cov.olse = function(e.cov, sigma0, n){
    alpha_hat = 1-(sum(diag(e.cov))^2 * norm(sigma0, type = "F")^2 / n) / (norm(e.cov, type = "F")^2 * norm(sigma0, type = "F")^2 - sum(diag(e.cov %*% 
    sigma0))^2)
    beta_hat = (1-alpha_hat) * sum(diag(e.cov %*% sigma0)) / norm(sigma0, type = "F")^2
    alpha_hat * e.cov + beta_hat * sigma0
}

js.classic.mean.f = function(X, e.mean, e.cov, a){
    t1 = Sys.time()
    p = dim(X)[2]; 
    n = dim(X)[1]; 
    if(p < n){c(as.vector(as.double(1 - ((p - 2) / (n - p + 3))/(t(e.mean) %*% ginv(e.cov, tol = tol.ginv) %*% e.mean)) * e.mean), difftime(Sys.time(), t1, units = "secs"))}else{NA}
} # classical James-Stein mean, works for p \leq n

js.largep.mean.f = function(X, e.mean, e.cov, a){
    t1 = Sys.time()
    p = dim(X)[2]; 
    n = dim(X)[1]; 
    ge.cov = ginv(e.cov, tol = tol.ginv)
    c(as.vector((diag(1, p) - (a * e.cov %*% ge.cov) / as.double(t(e.mean) %*% ge.cov %*% e.mean)) %*% e.mean), difftime(Sys.time(), t1, units = "secs"))
} # James-Stein Estimator with large p and small n


js.largep.plus.mean.f = function(X, e.mean, e.cov, a){
    t1 = Sys.time()
    p = dim(X)[2]; 
    n = dim(X)[1]; 
    ge.cov = ginv(e.cov, tol = tol.ginv)
    c(as.vector((diag(1, p) - e.cov %*% ge.cov) %*% e.mean + max(0, 1 - a / as.double(t(e.mean) %*% ge.cov %*% e.mean)) * e.cov %*% ge.cov %*% e.mean), difftime(Sys.time(), t1, units = "secs"))
} # James-Stein Estimator with large p and small n, positive parttype

our.oracle = function(e.mean, t.mean, mu0, icov.matr, params = FALSE){
    t1 = Sys.time()
    # i.cov = ginv(cov.matr, tol = tol.ginv)
    i.cov = icov.matr
    eie = t(e.mean) %*% i.cov %*% e.mean
    ei0 = t(e.mean) %*% i.cov %*% mu0
    ni0 = t(mu0)    %*% i.cov %*% mu0
    eit = t(e.mean) %*% i.cov %*% t.mean
    ti0 = t(t.mean) %*% i.cov %*% mu0
    denominator = eie * ni0 - ei0^2
    alpha_n = as.double((eit %*% ni0 - ti0 %*% ei0) / denominator)
    beta_n  = as.double((eie %*% ti0 - ei0 %*% eit) / denominator)
    if(!params){
        c(as.vector(alpha_n * e.mean + beta_n * mu0), difftime(Sys.time(), t1, units = "secs"))
    }else{
        c(alpha = alpha_n, beta = beta_n)
    }
}

our.lim = function(n, t.mean, it.cov, mu0, e.mean, params = FALSE, Lgamma = 0.0){
    t1 = Sys.time()
    p  = dim(it.cov)[1]
    cc = p^(1-Lgamma) / n
    # i.cov = ginv(t.cov, tol = tol.ginv)
    i.cov = it.cov
    mu0icovmu0 = t(mu0) %*% i.cov %*% mu0
    alpha_hat = as.double(1 - cc * mu0icovmu0 / ((cc + t(t.mean) %*% i.cov %*% t.mean)%*% mu0icovmu0 - (t(t.mean) %*% i.cov %*% mu0)^2))
    beta_hat = as.double((1 - alpha_hat) * (t(t.mean) %*% i.cov %*% mu0) / mu0icovmu0)
    
    if(!params){
        c(as.vector(alpha_hat * e.mean + beta_hat * mu0), difftime(Sys.time(), t1, units = "secs")) # I think here should be e.mean !!!
    }else{
        c(alpha = alpha_hat, beta = beta_hat)
    }
}

# our.bonafide.old = function(X, e.mean, e.cov, mu0, params = FALSE, Lgamma = 0){ ## gooood!
    # t1 = Sys.time()
    # p     = dim(X)[2]
    # n     = dim(X)[1]
    # cc    = (p^(1 - Lgamma)) / n
    # i.cov = ginv(e.cov, tol = tol.ginv)
    
    # q00 = t(mu0) %*% i.cov %*% mu0
    # qnn = t(e.mean) %*% i.cov %*% e.mean
    # q0n = t(e.mean) %*% i.cov %*% mu0
    
    # if(cc < 1){
    		# alpha_hat = as.double(1 - p^(Lgamma) * cc / (1 - cc) * q00 / (((cc * (p^Lgamma - 1))/(1 - cc) + qnn) * q00 - (q0n)^2))
    	# }else if(cc > 1){
   		# alpha_hat = as.double(1 - p^(Lgamma) / (cc - 1) * q00 / (((p^Lgamma - 1)/(cc - 1) + qnn) * q00 - q0n^2))
    	# }
# #    c(as.vector(alpha_hat * e.mean + beta_hat * mu0), difftime(Sys.time(), t1, units = "secs"))
	# beta_hat = as.double((1 - alpha_hat) * (t(e.mean) %*% i.cov %*% mu0) / (t(mu0) %*% i.cov %*% mu0))
    # if(!params){
        # c(as.vector(alpha_hat * e.mean + beta_hat * mu0), difftime(Sys.time(), t1, units = "secs"))
    # }else{
        # c(alpha = alpha_hat, beta = beta_hat)
    # }
# }

our.bonafide = function(X, e.mean, ie.cov, mu0, params = FALSE){ ## gooood!
    t1 = Sys.time()
    p     = dim(X)[2]
    n     = dim(X)[1]
    cc    = p / n
    i.cov = ie.cov
    # i.cov = ginv(e.cov, tol = tol.ginv)
    
    q00 = t(mu0) %*% i.cov %*% mu0
    qnn = t(e.mean) %*% i.cov %*% e.mean
    q0n = t(e.mean) %*% i.cov %*% mu0
    
    # gcc  = if(cc < 1){cc / (1 - cc)}else{1 / (cc - 1)}
    gcc  = if(cc < 1){cc / (1 - cc)}else{1 / (cc - 1)}
	alpha_hat = as.double(((qnn - gcc) * q00 - q0n^2) / (qnn * q00 - q0n^2))
	beta_hat = as.double((1 - alpha_hat) * (q0n / q00))
    if(!params){
        c(as.vector(alpha_hat * e.mean + beta_hat * mu0), difftime(Sys.time(), t1, units = "secs"))
    }else{
        c(alpha = alpha_hat, beta = beta_hat)
    }
}


wtcm.mean.f = function(X, e.mean, e.cov, a, param = FALSE){
    t1 = Sys.time()
	c.dim = ncol(X)
	n  = nrow(X)
	e  = rep(1, c.dim)
	Y0 = n^2 / (n-1) * e.mean %*% t(e.mean)-t(X) %*% X / (n-1)
	Y1 = sum(diag(a %*% Y0)) / c.dim
	Y2 = sum(diag(a %*% e.cov)) / c.dim
	Y3 = t(e) %*% a %*% Y0 %*% a %*% e / c.dim / (t(e) %*% a %*% e)
	Y4 = (t(e) %*% a %*% e.mean)/(t(e) %*% a %*% e)
    alpha = as.double((Y1 - Y3) / (Y1 + Y2 - Y3))
    beta  = as.double(Y2 * Y4 / (Y1 + Y2 - Y3))
    if(!param){
	    c(as.vector(alpha * as.vector(e.mean) + beta * e), difftime(Sys.time(), t1, units = "secs"))
	# 	obj = (Y1-Y3)/(Y1+Y2-Y3)*bX+Y2/(Y1+Y2-Y3)*Y4*e
	}else{
		c(alpha = alpha, beta = beta)
	}
}


normality.oracle = function(Lcc, Lmu.not, Ltrue.mean, Ltrue.cov, Lgamma){
# Lcc       = cc 
# Lmu.not   = mu.not
# Ltrue.mean= true.mean
# Ltrue.cov = true.cov
# Lgamma    = gammap
    i.cov = solve(Ltrue.cov)
    q00 = t(Lmu.not)    %*% i.cov %*% Lmu.not
    qnn = t(Ltrue.mean) %*% i.cov %*% Ltrue.mean
    q0n = t(Lmu.not) %*% i.cov %*% Ltrue.mean
    d = q00 * qnn - q0n^2
    
    denom = (Lcc * q00 + d)^4
    sigma_alpha = ((Lcc * q00 - d)^2 * q00 * d + Lcc * d^2 * q00^2) / denom
    sigma_beta = ((d - Lcc * q00)^2 * q0n^2 * qnn + (Lcc * q0n^2 - Lcc * d - d * qnn)^2 * q00 + 2 * (Lcc * q0n^2 - Lcc * d - d * qnn) * (d - Lcc * q00) * q0n^2 + Lcc * d^2 * q0n^2) / denom
    c(sigma_alpha, sigma_beta)
}

normality.bf = function(Lcc, Lmu.not, Ltrue.mean, Ltrue.cov, Lgamma){
#Lcc       = cc 
#Lmu.not   = mu.not
#Ltrue.mean= true.mean
#Ltrue.cov = true.cov
#Lgamma    = gammap

    i.cov = solve(Ltrue.cov)
    q00 = t(Lmu.not)    %*% i.cov %*% Lmu.not
    qnn = t(Ltrue.mean) %*% i.cov %*% Ltrue.mean
    q0n = t(Lmu.not) %*% i.cov %*% Ltrue.mean
    
    s = qnn - q0n^2 / q00
    R = q0n / q00
    
    sigmas = 2 * (Lcc + 2 * s) + 2 / abs(1 - Lcc) * (Lcc + s)^2
    sigma_alpha = (Lcc^2 * sigmas) / ((Lcc + s)^4)
    sigma_beta = sigma_alpha * R^2 + (Lcc^2 / ((Lcc + s)^2)) * (1 + (s + Lcc) / abs(1-Lcc)) / q00
    c(sigma_alpha, sigma_beta)
}
