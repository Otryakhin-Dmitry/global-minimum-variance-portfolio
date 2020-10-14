

n<-8e2 # number of realizations
p<-0.3*n # number of assets
w_0 <- rep(1/p,p)
gamma<-5e4

#################################################


#### Check limiting distribution of alpha_hat_star_c - alpha_star ####
# Remark 1
Sigma<-matrix(data=0, nrow=p, ncol=p)
diag(Sigma) <- 1

vect_al <-
  replicate(n=5e2, {

    x <-matrix(data = rnorm(n=n*p), nrow = p, ncol = n)
    alpha_hat_star_c(gamma=gamma, x=x, b=w_0) -
      alpha_star(gamma=gamma, mu=rep(0,p), Sigma=Sigma, b=w_0, c=p/n)
  })

Lb <- V_b(Sigma=Sigma, b=w_0)/V_GMV(Sigma=Sigma) - 1
v_ex_up <- 2*(1-p/n)*(p/n)^2*(Lb+1)
v_ex_dwn<- ((1-p/n)*R_b(mu=rep(0,p), b=w_0) + p/n)^4
v_ex_l  <- ((2-p/n)*Lb+p/n)

var_al <- v_ex_up/v_ex_dwn*v_ex_l

# Plot densities
plot(density(vect_al), xlim=c(-1,1))
points(x=seq(-1,1, by=0.1), y=dnorm(x=seq(-1,1, by=0.1), sd=sqrt(var_al)))

# simplified alphas
vect_al_simple <- sqrt(n)*
    replicate(n=1e3, {

        x <-matrix(data = rnorm(n=n*p), nrow = p, ncol = n)
        V_GMV <- V_hat_GMV(x=x)
        V_hat_c <- V_GMV/(1-p/n)
        V_hat_b <- V_hat_b(x=x, b=w_0)

        V_b <- V_b(Sigma=Sigma, b=w_0)
        V_GMV <- V_GMV(Sigma=Sigma)

        res<-
            (1-p/n)*(V_hat_b-V_hat_c)/(p/n+(1-p/n)*(V_hat_b-V_hat_c))-
            (1-p/n)*(V_b-V_GMV)/(p/n+(1-p/n)*(V_b-V_GMV))
        res
    })

plot(density(vect_al_simple), xlim=c(-0.01,0.01))

#### Check convergence for T_alpha and (44) ####

vect_T_al <-
replicate(n=1e3,{
    x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)
    T_alpha(gamma=gamma, x=x, w_0=w_0, c=p/n)
})

vect_plot <- vect_T_al[is.na(vect_T_al)==FALSE & vect_T_al>-10 & vect_T_al<10]

plot(density(vect_plot), xlim=c(-5,5))
points(x=seq(-5,5, by=0.1), y=dnorm(x=seq(-5,5, by=0.1)))












