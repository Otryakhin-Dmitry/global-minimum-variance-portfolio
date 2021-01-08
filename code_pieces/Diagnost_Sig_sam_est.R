
n<-8e2 # number of realizations
p<-.5*n # number of assets

x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)

t1 <- Sys.time()
  replicate(n=9e2, {test1 <- cov(t(x))})
t2 <- Sys.time()

x1 <- t(x)
t1 <- Sys.time()
  replicate(n=9e2, {test1 <- cov(x1)})
t2 <- Sys.time()

t3 <- Sys.time()
  replicate(n=9e2, {test2 <- Sigma_sample_estimator(x)})
t4 <- Sys.time()

t2-t1
t4-t3

test1[10:16,10:16]

# Result: Sigma_sample_estimator is almost 2 times faster than cov.


#####

Sigma_sample_estimator_df <- function(x) {

  p <- nrow(x)
  n <- ncol(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  a <- .rowMeans(x, m=p, n=n, na.rm = TRUE)
  a_x_size <- matrix(rep(a,n),nrow=p, ncol=n)
  tcrossprod(x-a_x_size)/(ncol(x)-1)
}


Sigma_sample_estimator <- function(x) {

  p <- nrow(x)
  n <- ncol(x)

  a <- .rowMeans(x, m=p, n=n, na.rm = TRUE)
  a_x_size <- matrix(rep(a,n),nrow=p, ncol=n)
  tcrossprod(x-a_x_size)/(ncol(x)-1)
}


x1 <- t(x)
t1 <- Sys.time()
replicate(n=9e2, {test1 <- Sigma_sample_estimator(x)})
t2 <- Sys.time()

t3 <- Sys.time()
replicate(n=9e2, {test2 <- Sigma_sample_estimator_df(x)})
t4 <- Sys.time()

t2-t1
t4-t3



#####



#### Exper with dataframes

p<-5 # number of assets
n<-1e1 # number of realizations
x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)
x_df <- as.data.frame(x)
rownames(x_df) <- c('a', 'b', 'c', 'd', 'e')
Sigma_sample_estimator(x)
Sigma_sample_estimator(x_df)


p <- nrow(x_df)
n <- ncol(x_df)
if (is.data.frame(x_df)) x_df <- as.matrix(x_df)
a <- .rowMeans(x_df, m=p, n=n, na.rm = TRUE)
a_x_size <- matrix(rep(a,n),nrow=p, ncol=n)
tcrossprod(x_df-a_x_size)/(ncol(x_df)-1)
