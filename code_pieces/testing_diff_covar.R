 p<-8e2 # number of assets
 n<-8e2 # number of realizations
 reps<-200
 x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)


 tt<-Sys.time()
 replicate(reps,{ a<-Sigma_sample_estimator(x)})

 t0<-Sys.time()-tt


 Sigma_sample_estimator_new <- function(x) {

   p <- nrow(x)
   n <- ncol(x)
   a <- .rowMeans(x, m=p, n=n, na.rm = TRUE)
   a_x_size <- matrix(rep(a,n),nrow=p, ncol=n)
   tcrossprod(x-a_x_size)/(ncol(x)-1)
 }

 t1<-Sys.time()
 replicate(reps, {a<-Sigma_sample_estimator_new(x)})

 t2<-Sys.time()-t1

 t3<-Sys.time()
 replicate(reps, {a<-fcprd(x)})

 t4<-Sys.time()-t3

 t0
 t2
 t4
