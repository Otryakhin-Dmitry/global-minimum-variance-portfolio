
#### Sigma simple estimator (5) in EU tests paper
#' Simple sample covariance estimator
#'
#' It computes covariance of matrix x as follows:
#' \deqn{\hat\Sigma_n = \frac{1}{n-1} \sum_{j=1}^n (x_j - \bar x_n)(x_j - \bar x_n)^T
#' ,\quad \bar x_n = \frac{1}{n} \sum_{j=1}^n x_j}
#'
#' @param x a numeric matrix or a data frame in which columns are independent realizations of asset returns
#' @references \insertRef{BDOPS2020}{hdsp}
#' @examples
#' p<-5 # number of assets
#' n<-1e1 # number of realizations
#'
#' x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)
#' Sigma_sample_estimator(x)
#' @importFrom stats rnorm
#' @export
Sigma_sample_estimator <- function(x) {

  p <- nrow(x)
  n <- ncol(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  a <- .rowMeans(x, m=p, n=n, na.rm = TRUE)
  a_x_size <- matrix(rep(a,n),nrow=p, ncol=n)
  tcrossprod(x-a_x_size)/(ncol(x)-1)
}


#### Q estimator (page 3, IEEE)
Q_hat_n <- function(x){

  SS <- Sigma_sample_estimator(x)
  invSS <- solve(SS)
  Ip <- rep.int(1, nrow(x))

  invSS - (invSS %*% Ip %*% t(Ip) %*% invSS)/as.numeric(t(Ip) %*% invSS %*% Ip)
}


Q <- function(Sigma){

  invSS <- solve(Sigma)
  Ip <- rep.int(1, nrow(Sigma))
  invSS - (invSS %*% Ip %*% t(Ip) %*% invSS)/as.numeric(t(Ip) %*% invSS %*% Ip)
}



#### S's

s_hat <- function(x) {

  a <- rowMeans(x, na.rm = TRUE)
  as.numeric(t(a) %*% Q_hat_n(x) %*% a)
}

s <- function(mu, Sigma) {

  as.numeric(t(mu) %*% Q(Sigma) %*% mu)
}

s_hat_c <- function(x) as.numeric((1-nrow(x)/ncol(x))*s_hat(x) - nrow(x)/ncol(x))

#### R_GMV (page 5, IEEE)
# R_GMV. The deterministic value. The expected return of the GMV portfolio

R_GMV <- function(mu, Sigma){

  p <- length(mu)
  invSS <- solve(Sigma)
  Ip <- rep.int(1, p)

  as.numeric((t(Ip) %*% invSS %*% mu)/(t(Ip) %*% invSS %*% Ip))
}

# Its estimator
R_hat_GMV <- function(x){

  a <- rowMeans(x, na.rm = TRUE)
  SS <- Sigma_sample_estimator(x)
  invSS <- solve(SS)
  Ip <- rep.int(1, nrow(x))

  as.numeric((t(Ip) %*% invSS %*% a)/(t(Ip) %*% invSS %*% Ip))
}


R_b <- function(mu, b) as.numeric(b %*% mu)

R_hat_b <- function(x, b) as.numeric(b %*% rowMeans(x, na.rm = TRUE))

#### V's

V_b <- function(Sigma, b) as.numeric(t(b) %*% Sigma %*% b)

V_hat_b <- function(x, b) {

  Sigma <- Sigma_sample_estimator(x)
  as.numeric(t(b) %*% Sigma %*% b)
}

V_GMV <- function(Sigma){

  as.numeric(1/(rep.int(1, nrow(Sigma)) %*% solve(Sigma) %*% rep.int(1, nrow(Sigma))))
}

V_hat_GMV <- function(x){

  Sigma <- Sigma_sample_estimator(x)
  as.numeric(1/(rep.int(1, nrow(Sigma)) %*% solve(Sigma) %*% rep.int(1, nrow(Sigma))))
}

V_hat_c <- function(x) {V_hat_GMV(x)/(1-nrow(x)/ncol(x))}

#### alphas, B and A expressions

# In case of GMV portfolio one needs to set gamma=infty
alpha_star <- function(gamma, mu, Sigma, b, c){

  R_GMV <- R_GMV(mu, Sigma)
  R_b <- R_b(mu, b)
  V_GMV <- V_GMV(Sigma)
  V_b <- V_b(Sigma, b)
  s <- s(mu, Sigma)

  Exp1 <- (R_GMV-R_b)*(1+1/(1-c))/gamma
  Exp2 <- (V_b-V_GMV)
  Exp3 <- s/(gamma^2)/(1-c)
  numerator <- Exp1 + Exp2 + Exp3

  Exp4 <- V_GMV/(1-c)
  Exp5 <- -2*(V_GMV + (R_b - R_GMV)/(gamma*(1-c)))
  Exp6 <- ((s+c)/(1-c)^3)/(gamma^2)
  denomenator <- Exp4 + Exp5 + Exp6 + V_b

  as.numeric(numerator/denomenator)
}


# In case of GMV portfolio one needs to set gamma=infty
alpha_hat_star_c <- function(gamma, x, b){

  R_GMV <- R_hat_GMV(x)
  R_b <- R_hat_b(x, b)
  V_GMV <- V_hat_GMV(x)
  V_b <- V_hat_b(x, b)

  c <- nrow(x)/ncol(x)
  s <- s_hat_c(x)

  V_c <- V_GMV/(1-c)

  Exp1 <- (R_GMV-R_b)*(1+1/(1-c))/gamma
  Exp2 <- (V_b-V_c)
  Exp3 <- s/(gamma^2)/(1-c)
  numerator <- Exp1 + Exp2 + Exp3

  Exp4 <- V_c/(1-c)
  Exp5 <- -2*(V_c + (R_b - R_GMV)/(gamma*(1-c)))
  Exp6 <- ((s+c)/(1-c)^3)/gamma^2
  denomenator <- Exp4 + Exp5 + Exp6 + V_b

  as.numeric(numerator/denomenator)
}

alpha_hat_star_c_fast <- function(gamma, c, s, b, R_GMV, R_b, V_GMV, V_b){

  V_c <- V_GMV/(1-c)

  Exp1 <- (R_GMV-R_b)*(1+1/(1-c))/gamma
  Exp2 <- (V_b-V_c)
  Exp3 <- s/(gamma^2)/(1-c)
  numerator <- Exp1 + Exp2 + Exp3

  Exp4 <- V_c/(1-c)
  Exp5 <- -2*(V_c + (R_b - R_GMV)/(gamma*(1-c)))
  Exp6 <- ((s+c)/(1-c)^3)/gamma^2
  denomenator <- Exp4 + Exp5 + Exp6 + V_b

  as.numeric(numerator/denomenator)
}

alpha_star_GMV <- function(Sigma, b, c){

  V_GMV <- V_GMV(Sigma)
  V_b <- V_b(Sigma, b)

  numer <- (1-c)*(V_b-V_GMV)
  as.numeric(numer/(numer + c*V_GMV))
}

# simplified alphas when gamma=Inf
alpha_hat_star_c_GMV <- function(x, b, c = nrow(x)/ncol(x)){

  V_GMV <- V_hat_GMV(x)
  V_b <- V_hat_b(x, b)
  c <- nrow(x)/ncol(x)
  V_c <- V_GMV/(1-c)

  numer <- (1-c)*(V_b-V_c)
  as.numeric(numer/(numer + c*V_c))
}

# A and B expressions
B_hat <- function(gamma, x, b){

  R_GMV <- R_hat_GMV(x)
  R_b <- R_hat_b(x, b)
  V_GMV <- V_hat_GMV(x)
  V_b <- V_hat_b(x, b)

  c <- nrow(x)/ncol(x)
  s <- s_hat_c(x)

  V_c <- V_GMV/(1-c)

  Exp4 <- V_c/(1-c)
  Exp5 <- -2*(V_c + (R_b - R_GMV)/(gamma*(1-c)))
  Exp6 <- ((s+c)/(1-c)^3)/gamma^2
  denomenator <- Exp4 + Exp5 + Exp6 + V_b

  as.numeric(denomenator)
}

Var_alpha_simple <- function(Sigma, b, mu, n){

  c <- nrow(Sigma)/n
  V_b <- V_b(Sigma, b)
  V_GMV <- V_GMV(Sigma)
  Lb <- V_b/V_GMV - 1
  # R_b <- R_b(mu, b)

  numer <- 2*(1-c)*c^2*(Lb+1)
  denom <- ((1-c)*Lb+c)^4
  multip<- (2-c)*Lb +c

  numer / denom * multip
}

# BDOPS2020, under formula 16
Omega.Lest <- function(s.est, cc, gamma, V.est, L, Q.est, eta.est){

  (((1-cc)/(s.est+cc)+ (s.est+cc)/gamma)/gamma + V.est)*(1-cc)*L%*%Q.est%*%t(L)+
  gamma^{-2}*(2*(1-cc)*cc^3/(s.est+cc)^2+ 4*(1-cc)*cc*s.est*(s.est+2*cc)/(s.est+cc)^2 +
              2*(1-cc)*cc^2*(s.est+cc)^2/(s.est^2)-s.est^2)*eta.est%*%t(eta.est)
}



