# Shrinkage portfolio
EUShrinkPortfolio <- function(x, gamma, type, subtype, ...) {

  if(!is.numeric(gamma) || is.na(gamma)) stop("gamma is not numeric")
  if(gamma==Inf) stop("GMVP methods are absent")

  if(type=='mean') {
    if(subtype=='James-Stein') {
      output <- new_ExUtil_portfolio_mean_JamesStein(x=x, gamma=gamma, ...)
    }
    if(subtype=='Bayes-Stein') {
      output <- new_ExUtil_portfolio_mean_BayesStein(x=x, gamma=gamma, ...)
    }
  }

  if(type=='weights') {
    output <- new_ExUtil_portfolio(x=x, gamma=gamma, ...)
  }
  output
}

# Examples
n<-3e2 # number of realizations
p<-.5*n # number of assets
b<-rep(1/p,p)
gamma<-1

x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)

test <- EUShrinkPortfolio(x=x, gamma=gamma, type='weights', b=b)
str(test)

test <- EUShrinkPortfolio(x=x, gamma=gamma, type='mean', subtype='Bayes-Stein')
str(test)

test <- EUShrinkPortfolio(x=x, gamma=gamma, type='mean', subtype='James-Stein')
str(test)



##########################################################################


fdisp <- cxxfunction(signature(AA = "matrix"),
                     'using namespace std;
                      double gamma=100;
                      std::string type;

                      if( gamma == 0 ) {

                      } else {
                      if( type == 'mean' ) {

                      } else {
                      }
                      }

                      return List::create(Named("crossprod(A)") = 4);'
                     , "RcppEigen")





