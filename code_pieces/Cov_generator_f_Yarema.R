



eig = 0.1*exp(q*cc*seq(0,1,length=p)) # eig.v[i.eig] was replaced

Z  <- array(rnorm(p*p,0,1),dim=c(p,p))
QR <- qr(Z)
Q  <- qr.Q(QR)
R  <- qr.R(QR)
D  <- diag(as.vector(diag(R)))
D  <-D%*%solve(chol(D%*%D))
H  <-Q%*%D
E  <-diag(eig)
Sigma<-H%*%E%*%t(H)







