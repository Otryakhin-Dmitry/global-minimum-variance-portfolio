rm(list=ls(all=TRUE))
# setwd("D:\\Private\\IdeaDescription\\Nestor\\")
setwd("/Users/ookhrin/Dropbox/HighMean_BodOkhPar/2016_code/")
# install.packages(c("snow", "foreach", "doSNOW", "MASS"))
library(snow)
library(foreach)
library(doSNOW)
library(MASS)
library(parallel)
source("est_methods.r")

ei.vals      = c(1, 3, 10)
ei.val.probs = c(1/5, 2/5) # should be by one element smaller, than the previous vector
repss = 100
tol.ginv = sqrt(.Machine$double.eps) # tolerance level for the generalized inverse

ps = c(seq(10, 100, 10), seq(125, 250, 25), seq(300, 500, 50)) # dimensions
cc.m = c(0.1,0.5,0.9,1.5,2.0)
model_setup = "Usual" # LambdaInf" # "
gammap = 0.0
servers = max(detectCores()-1, 1)

# cc = 1.5
for(cc in rev(cc.m))
for(model_setup in c("Usual", "LambdaInf")){
	#cc = 2.0
        dists = list()
        print(paste("c = ", cc, sep = ""))
        # for(p in ps[(which(c.dims == 225)+1):length(c.dims)]){
        for(p in ps){
        #	p = 200
            big.t = Sys.time()
            # cc = 2
            n = round(p^(1-gammap) / cc)
            print(paste("p = ", p, ", n = ", n, sep = ""))
           
            tss = try(qr.Q(qr(array(rnorm(p), dim=c(p, p))))) 
            while(class(tss) == "try-error"){tss = try(qr.Q(qr(array(rnorm(p), dim=c(p,p)))))}
            eigen.vec = tss
            if(model_setup == "Usual"){ # usual setup, with moderate eigenvalues        
                eigen.val = 0;for(i in 1:length(ei.val.probs))eigen.val = c(eigen.val, rep(ei.vals[i], round(p * ei.val.probs[i])));eigen.val = c(eigen.val, rep(ei.vals[length(ei.vals)], p - sum(round(p * ei.val.probs))));eigen.val = eigen.val[-1]
            }else if(model_setup == "LambdaInf"){ # one eigenvalue is \lambda_p = p
                eigen.val = p;for(i in 1:length(ei.val.probs))eigen.val = c(eigen.val, rep(ei.vals[i], round((p - 1) * ei.val.probs[i])));eigen.val = c(eigen.val, rep(ei.vals[length(ei.vals)], (p-1) - sum(round((p-1) * ei.val.probs))));eigen.val = sort(eigen.val)
            }
            true.cov = eigen.vec %*% diag(eigen.val) %*% t(eigen.vec) 
            
            # eigen.vec2 = eigen(rWishart(1, p^2, diag(p))[, , 1])$vectors
            # true.cov2 = eigen.vec2 %*% diag(eigen.val) %*% t(eigen.vec2) 

            #norms = 0    
            #for(i in 1:1000){
            #    norms = c(norms, sum((runif(i,-1,1)/sqrt(i))^2))
            #}
            #norms = norms[-1]
            #plot(norms)
            
            
            ### big norm
            # true.mean = sample(c(-1,1), p, replace = TRUE) # runif(p,-1,1) / sqrt(p) # runif(p,-0.2,0.2)
            ### small norm
            true.mean = runif(p,-1,1) / sqrt(p) # runif(p,-0.2,0.2)
            # mu.not = as.matrix(rep(1, p))
            ### big norm
            # mu.not    = rep(1, p) # runif(p,-1,1) / sqrt(p) # runif(p,-0.2,0.2) # rep(1/p, p) # rep(sqrt(sum(true.mean^2)/p), p) # 
            ### small norm
            # mu.not    = runif(p,-1,1) / sqrt(p) # runif(p,-0.2,0.2) # rep(1/p, p) # rep(sqrt(sum(true.mean^2)/p), p) # 

            mu.not    = true.mean


            sigma.not = diag(1, p, p)
        
            cl = makeCluster(rep("localhost", servers), type = "SOCK")
            clusterEvalQ(cl, library(MASS))
            registerDoSNOW(cl)
            
            l.dist = foreach(i=1:repss, .combine=rbind) %dopar% {
                X = matrix(rnorm(p * n), ncol = p) %*% chol(true.cov) + matrix(rep(true.mean, each = n), ncol = p)
                while(class(try(ginv((1/n)*(t(X)-matrix(rep(colMeans(X), n), ncol = n))%*%t(t(X)-matrix(rep(colMeans(X), n), ncol = n))))) == "try-error"){
                    X = matrix(rnorm(p * n), ncol = p) %*% chol(true.cov) + matrix(rep(true.mean, each = n), ncol = p)
                }
                
                t1 = Sys.time()            
                    sample.mean = colMeans(X)
                t.s.mean = difftime(Sys.time(), t1)
                sample.mean = as.matrix(sample.mean)
                sample.cov = (1/n)*(t(X)-matrix(rep(sample.mean, n), ncol = n))%*%t(t(X)-matrix(rep(sample.mean, n), ncol = n)) # Nestor
                # olse.cov = cov.olse(sample.cov, sigma.not, n)
                est.means = cbind(
                    "True" = c(true.mean, 0), 
                    "Sample" = c(as.vector(sample.mean), t.s.mean), 
                    "JS" = if(cc < 1)                 {as.vector(js.classic.mean.f(    X, sample.mean, sample.cov))}else{}, 
                    "JS p>n" = if(cc >= 1)            {as.vector(js.largep.mean.f(     X, sample.mean, sample.cov, a = (n - 2)/(p - n + 3)))}else{}, 
                    "JS +" = if(cc >= 1)              {as.vector(js.largep.plus.mean.f(X, sample.mean, sample.cov, a = (n - 2)/(p - n + 3)))}else{}, 
                    "Wang et al. (Q=1)" = if(cc >= 1) {as.vector(wtcm.mean.f(          X, sample.mean, sample.cov, a = ginv(sample.cov)))}else{},
                    "Wang et al. (Q=SigmaInv)" = if(cc >= 1) {as.vector(wtcm.mean.f(          X, sample.mean, sample.cov, a = diag(1, p)))}else{},
                    "BOP (oracle)" =                   as.vector(our.oracle(           sample.mean, true.mean, mu.not, true.cov)),
                    "BOP (as)" =                       as.vector(our.lim(              n, true.mean, true.cov, mu.not, sample.mean, Lgamma = gammap)),
                    "BOP (bf)" =                       as.vector(our.bonafide(         X, sample.mean, sample.cov, mu.not, Lgamma = gammap))
                )
                
                ##### alpha and beta from our methods
                # aaa.as = as.vector(our.lim(              n, true.mean, true.cov, mu.not, sample.mean, Lgamma = gammap))[-(p+1)]
                # our.lim(              n, true.mean, true.cov, mu.not, sample.mean, Lgamma = gammap, params = TRUE)
                # aaa.bf = as.vector(our.bonafide(         X, sample.mean, sample.cov, mu.not, Lgamma = gammap))[-(p+1)]
                # our.bonafide(         X, sample.mean, sample.cov, mu.not, Lgamma = gammap, params = TRUE)
                # t(aaa.as-true.mean) %*% solve(true.cov) %*% t(t(aaa.as-true.mean))
                # t(aaa.bf-true.mean) %*% solve(true.cov) %*% t(t(aaa.bf-true.mean))
                
                rbind(diag(t(est.means[1:p,-1] - matrix(rep(true.mean, dim(est.means)[2]-1), nrow = p)) %*% solve(true.cov) %*% (est.means[1:p,-1] - matrix(rep(true.mean, dim(est.means)[2]-1), nrow = p)))/p, times = tail(est.means, 1)[-1])
        		}
        		where.times = seq(2,dim(l.dist)[1], by = 2)
        		stopCluster(cl)
        		dists[[which(p == ps)]] = rbind(l.dist[-where.times,], "AET" = colMeans(l.dist[where.times,]))
        		print(paste(p, " (", difftime(Sys.time(), big.t, units = "secs"), " sec.)", sep = ""))
    		}    
		# save(dists, file = paste("MeanEst_Dists_new_", model_setup, "_c", round(cc * 10), ".dat", sep = ""))
		save(dists, file = paste("MeanEst_Dists_new_mu0EQmun_", model_setup, "_c", round(cc * 10), ".dat", sep = ""))
}



        # ps = c(seq(10, 100, 10), seq(125, 250, 25), seq(300, 500, 50))
        
        means = 0
        for(i in 1:length(dists)){
            means = rbind(means, colMeans(dists[[i]][-dim(dists[[i]])[1], ]))
        }
        means = means[-1,]
        rownames(means) = ps
        
        means = means * matrix(rep(ps, dim(means)[2]), nrow = length(ps))

        yrange = range(means, na.rm = TRUE)
        yrange[2] = min(5, yrange[2])
        plot(ps, means[,"Sample"], ylim = yrange, lwd = 3, col = "black", type = "l", xlab = "p", ylab = expression((hat(mu)-mu[n])^T*Sigma[n]^(-1)*(hat(mu)-mu[n])))
        lines(ps, means[,"JS p>n"], lwd = 3, col = "green3", lty = "dotted")
        lines(ps, means[,"JS +"], lwd = 3, col = "green3", lty = "dashed")
        lines(ps, means[,"Wang et al. (Q=I)"], lwd = 3, col = "yellow3", lty = "dotted")
		lines(ps, means[,"Wang et al. (Q=SigmaInv)"], lwd = 3, col = "yellow3", lty = "dashed")
		lines(ps, means[,"BOP (oracle)"], lwd = 3, col = "blue3", lty = "dashed")
		lines(ps, means[,"BOP (as)"], lwd = 3, col = "red3", lty = "dotdash")
		lines(ps, means[,"BOP (bf)"], lwd = 3, col = "plum", lty = "dotdash")
		
		legend("topright", bg ="white", legend = c("Sample", "JS p>n", "JS +", "Wang et al. (Q=I)", "Wang et al. (Q=SigmaInv)", "BOP (optimal)", "BOP (asymptotic)", "BOP (bona fide)"), 
                                    col = c("black", "green3", "green3", "yellow3", "yellow3", "blue3",  "red3", "plum"), 
                                    lty = c("solid", "dotted", "dashed", "dotted",  "dashed",  "dashed", "dotdash", "dotdash"), lwd = 3)

