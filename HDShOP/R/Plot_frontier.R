#' Plot the Bayesian efficient frontier \insertCite{bauder21}{HDShOP} and the
#' provided portfolios.
#'
#' The plotted  Bayesian efficient frontier is provided by \insertCite{@Eq. (8) in @bauder21;textual}{HDShOP}.
#' It is the set of optimal portfolios obtained by employing the posterior predictive
#' distribution on the asset returns. This efficient frontier can be used to assess
#' the mean-variance efficiency of various estimators of the portfolio weights.
#' The standard deviation of the portfolio return is plotted in the \eqn{x}-axis and the
#' mean portfolio return in the \eqn{y}-axis. The portfolios with the weights \eqn{\rm{w}} are added
#' to the plot by computing  \eqn{\sqrt{\rm{w}^\prime  S  w} }  and  \eqn{\rm w^\prime \bar x}.
#'
#' @inheritParams MVShrinkPortfolio
#' @param weights.eff matrix of portfolio weights. Each column contains p values of
#' the weights for a given portfolio. Default: equally weighted portfolio.
#'
#' @returns a ggplot object
#'
#' @references \insertAllCited{}
#' @examples
#' p = 150
#' n = 300
#' gamma <- 10
#' mu = seq(0.2,-0.2, length.out=p)
#' Sigma = RandCovMtrx(p=p)
#'
#' x <- t(MASS::mvrnorm(n=n , mu=mu, Sigma=Sigma))
#'
#' EW_port <- rep(1/p, length=p)
#' MV_shr_port <- new_MV_portfolio_weights_BDOPS21(x=x, gamma=gamma, b=EW_port, beta=0.05)$weights
#' GMV_shr_port <- new_MV_portfolio_weights_BDOPS21(x=x, gamma=Inf, b=EW_port, beta=0.05)$weights
#' MV_trad_port <- new_MV_portfolio_traditional(x=x, gamma=gamma)$weights
#' GMV_trad_port <- new_MV_portfolio_traditional(x=x, gamma=Inf)$weights
#'
#' weights.eff = cbind(EW_port, MV_shr_port, GMV_shr_port, MV_trad_port, GMV_trad_port)
#' colnames(weights.eff) <- c("EW", "MV_shr", "GMV_shr", "MV_trad", "GMV_trad")
#'
#'
#' Fplot <- plot_frontier(x, weights.eff)
#' Fplot
#' @export
plot_frontier <- function(x, weights.eff = rep(1/nrow(x), length=nrow(x)))
{
  portfolio.sd <- NULL
  portfolio.return <- NULL

  p.eff <- nrow(x)
  n.eff <- ncol(x)

  ones <- rep.int(1, p.eff); tones <- t(ones)

  Sigma.eff <- Sigma_sample_estimator(x)
  mu.eff  <- .rowMeans(x, m=p.eff, n=n.eff, na.rm = TRUE)

  points.weights = apply(weights.eff, 2, function(x) c(as.double(sqrt(t(x)%*%Sigma.eff%*%x)), as.double(t(x)%*%mu.eff)))
  points.weights = t(points.weights)

  points.weights = data.frame(portfolio.return = points.weights[,2],
                              portfolio.sd = points.weights[,1],
                              names = colnames(weights.eff))

  print(points.weights)

  Sigma.bayes <- Sigma.eff*(n.eff-1);
  c.bayes = 1/(n.eff-p.eff-1) + (2*n.eff-p.eff-1) /(n.eff*(n.eff-p.eff-1)*(n.eff-p.eff-2));

  iSigma.bayes<-solve(Sigma.bayes)
  V.est.bayes<-c.bayes/sum(tones%*%iSigma.bayes%*%ones) #estimated variance (cons.)
  Q.est.bayes<- iSigma.bayes-(iSigma.bayes%*%ones%*%tones%*%iSigma.bayes)/sum(tones%*%iSigma.bayes%*%ones)
  R.est.bayes<- (tones%*% iSigma.bayes %*% mu.eff)/sum(tones%*%iSigma.bayes%*%ones) #returns of EU portfolio

  V.front.bayes <- seq(sqrt(V.est.bayes), max(1.5*points.weights$portfolio.sd), length=100)
  R.front.bayes <- c(R.est.bayes) + sqrt( as.double(t(mu.eff)%*% Q.est.bayes %*% mu.eff)*(V.front.bayes^2-V.est.bayes)/c.bayes)

  points.eff =  data.frame(portfolio.return = R.front.bayes, portfolio.sd = V.front.bayes)

  pic <- ggplot2::ggplot(points.eff, ggplot2::aes(x=portfolio.sd, y=portfolio.return)) +
         ggplot2::geom_line(col = "royalblue4", linetype = 1, size=1) +
         ggplot2::theme_bw() +
         ggplot2::theme(legend.position = "none") +
         ggplot2::geom_point(data=points.weights, ggplot2::aes(x=portfolio.sd, y=portfolio.return), shape=2, size=3) +
         ggplot2::geom_text(data=points.weights, ggplot2::aes(label=names), col = "royalblue4", hjust=0, vjust=0)

  return(pic)
}
