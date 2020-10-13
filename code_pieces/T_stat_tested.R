


mult <- 10
p<-5*mult # number of assets
n<-1e1*mult # number of realizations
w_0 <- rep(1/5,p)
gamma<-0.8
#################################################

vect_T_al <-
replicate(n=1e3,{
    x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)
    T_alpha(gamma, x, w_0, c=p/n)
})

vect_plot <- vect_T_al[is.na(vect_T_al)==FALSE & vect_T_al>-10 & vect_T_al<10]

hist(vect_plot, freq=FALSE)
plot(density(vect_plot))
points(x=seq(-2,2, by=0.1), y=dnorm(x=seq(-2,2, by=0.1)))




