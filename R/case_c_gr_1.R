
V_hat_c_pgn_fast <- function(ones, invSS, tones, c) {

  V_hat_GMV <- as.numeric(1/(tones %*% invSS %*% ones))
  V_hat_GMV/c/(c-1)
}


# s_hat is the same for the 2 cases
s_hat_c_pgn <- function(x) {

  cc <- nrow(x)/ncol(x)
  as.numeric(cc*(cc-1)*s_hat(x) - cc)
}

# alpha in MV case, p>n
alpha_hat_star_c_pgn_fast <- function(gamma, c, s, R_GMV, R_b, V_c, V_b){

  Exp1 <- (R_GMV-R_b)*( 1 + 1/(c*(c-1)) )/gamma
  Exp2 <- (V_b-V_c)
  Exp3 <- s/(gamma^2)/(c*(c-1))
  numerator <- Exp1 + Exp2 + Exp3

  Exp4 <- c^2*V_c/(c-1)
  Exp5 <- -2*(V_c + (R_b - R_GMV)/(gamma*c*(c-1)))
  Exp6 <- ((s+c^2)/(c-1)^3)/gamma^2
  denomenator <- Exp4 + Exp5 + Exp6 + V_b

  as.numeric(numerator/denomenator)
}

# alpha in GMV case, p>n
alpha_hat_star_c_GMV_pgn_fast <- function(c, V_c, V_b){

  as.numeric( (V_b-V_c)/(c^2*V_c/(c-1) - 2*V_c + V_b) )
}








