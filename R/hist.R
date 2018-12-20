#' hctrial: A package for designing phase 2 clinical trials adjusting for heterogeneous populations.
#'
#' The hctrial package provides functions for designing phase 2 clinical trials that adjust for the
#' heterogeneity in the population.
#'
#' Two different ways are considered for designing a trial: based on known subgroups or based on historical data.
#'
#' For initializing a stratified trial, use \code{strat_start}.
#'
#' At interim, \code{strat_interim} should be used to adjust the trial.
#'
#' At the end of the study, \code{strat_end} is used to adjust the trial again.
#'
#' \code{hist_start}, \code{hist_interim} and \code{hist_end} work analogously, but are based on historical controls.
#'
#' @docType package
#' @name hctrial
NULL



#' Initializes a design based on historical controls before the start of the study.
#' @export
#' @param hist_data A data frame containing covariates and binary responses for historical controls.
#' @param formula A formula which is used for fitting a logistic regression model on the historical data.
#' @param phi The relation between the response rate under the null and the response rate under the interesting alternative.
#' "odds_ratio" assumes that the odds ratio (OR) between these response rates is constant with OR = \code{c1+1}.
#' "difference" assumes that the response rate under the alternative is \code{c1} higher than under the null.
#'  Can also be specified by the user by providing a function with arguments \code{c} and \code{x}.
#' @param c1 parameter for obtaining the response rate under the alternative, see description of phi.
#' @param modelfit Can be used instead of \code{formula} and \code{hist_data} to provide an arbitrary fitted model that is
#' compatible with \code{predict(modelfit, type="response")}. \code{formula} and \code{hist_data} are ignored if \code{modelfit} is
#' specified.
#' @param mean0 Optional: Can be used to overwrite the estimated average response rate under the null of the fitted model.
#' @param mean1 Optional: Can be used to overwrite the estimated average response rate under the alternative of the fitted model.
#' @param alpha Specified type I error of the trial.
#' @param beta Specified type II error of the trial.
#' @return A list returning the arguments of the function and the preliminary design for starting the stratified trial.
#' @examples
#' X <- abs(rnorm(1000, 0, 1))
#' Y <- rbinom(1000, 1, 1-exp(-X))
#' mydata <- data.frame("X" = X, "Y" = Y)
#' hist_start(mydata, Y~X, c1 = 2)


hist_start <- function(hist_data, formula, phi = "odds_ratio", c1, modelfit = NULL, mean0 = NULL, mean1 = NULL, alpha = 0.05, beta = 0.2) {

  if (is.null(modelfit)) {
    model <- as.formula(formula)
    modelfit <- glm(formula = formula, family = binomial, data = hist_data)
  }

  if (c1 <= 0) {
    stop("c1 must be greater than 0")
  }

  if (phi == "odds_ratio") {
    phic1 <- function(x) {(c1 * x + x) / (c1 * x + 1)}
  } else if (phi == "difference") {
    phic1 <- function(x) {min(x + c1, 1)}
  } else {
    phic1 <- function(x) {match.fun(phi)(c = c1, x = x)}
  }

  predhist <- predict(modelfit, newdata=hist_data, type="response")

  if (is.null(mean0)) {
    mean0 <- mean(predhist)
  }

  if (is.null(mean1)) {
    mean1 <- mean(sapply(predhist, phic1))
  }



  design_mat <- clinfun::ph2simon(pu = mean0, pa = mean1, ep1 =  alpha, ep2 = beta)
  ind_opt <- which.min(design_mat$out[,5])

  out <- design_mat$out[ind_opt,]

  out[4] <- out[4] - out[2]
  names(out)[4] <- "n2"


  r1 <- out[1]
  n1 <- out[2]
  r <- out[3]
  n2 <- out[4]

  sumup <- min(r,n1)
  sumindx <- ((r1+1):sumup)
  summands <- pbinom(r-sumindx, n2, mean0) * dbinom(sumindx, n1, mean0)
  type1 <- 1-(pbinom(r1, n1, mean0) + sum(summands))

  summands <- pbinom(r-sumindx, n2, mean1)*dbinom(sumindx, n1 , mean1)
  type2 <- pbinom(r1, n1, mean1)+sum(summands)

  out <- c(out, type1, type2)

  names(out)[c(7,8)] <- c("type1_est", "type2_est")

  return(list("modelfit" = modelfit, "phic1" = phic1, "phi" = phi, "c1" = c1, "mean_null" = mean0, "mean_alt" = mean1, "alpha" = alpha, "beta" = beta, "des_start" = out))

}

#' Adjust a design based on historical controls at interim using the covariate data of the patients accrued in stage 1.
#' @export
#' @param start An initialized design based on historical controls as returned by \code{hist_start()}.
#' @param stageone_data A dataframe containing the relevant covariate data of the patients accrued in stage 1.
#' @return A list returning the arguments of the function and the preliminary design of a trial based on historical controls adjusted at interim.
#' @examples
#' X <- abs(rnorm(1000, 0, 1))
#' Y <- rbinom(1000, 1, 1-exp(-X))
#' mydata <- data.frame("X" = X, "Y" = Y)
#' start <- hist_start(mydata, Y~X, c1 = 2)
#' n1 <- start$des_start[2]
#' X1 <- abs(rnorm(n1, 0, 1))
#' dataone <- data.frame("X" = X1)
#' hist_interim(start, dataone)


hist_interim <- function (start, stageone_data) {

  n1 <- start$des_start[2]



  if (nrow(stageone_data) != n1)
    stop("length of observations must equal n1")

  alpha <- start$alpha
  beta <- start$beta

  p0_observ <- predict(start$modelfit, newdata=stageone_data, type="response")
  p1_observ <- sapply(p0_observ, start$phic1)
  mean0 <- start$mean_null
  mean1 <- start$mean_alt

  vecbin <- GenBinomApps::pgbinom(0:n1, rep(1,n1), p1_observ)

  if (min(vecbin) >= beta)
    r1_up <- 0
    else
    r1_up <- max(which(vecbin < beta))-1

  r_cand <- n2_cand <- rep(NA,r1_up+1)

  for (r1 in 0:r1_up) {
    n2 <- 1
    stopp <- 0

    if ((1-GenBinomApps::pgbinom(r1, rep(1,n1), p0_observ) <= alpha)) {
      stopp <- 2
      r_cand[r1+1] <- r1
      n2_cand[r1+1] <- 0
    }


    while (stopp<1) {
      for (r in (r1+1):(r1+n2)) {
        sumup <- min(r,n1)
        sumindx <- ((r1+1):sumup)
        summands <- pbinom(r-sumindx, n2, mean0) * GenBinomApps::dgbinom(sumindx, size = rep(1,n1), prob = p0_observ)
        type1 <- 1-(GenBinomApps::pgbinom(r1, rep(1,n1), p0_observ) + sum(summands))

        summands <- pbinom(r-sumindx, n2, mean1)*GenBinomApps::dgbinom(sumindx, rep(1,n1) , p1_observ)
        type2 <- GenBinomApps::pgbinom(r1, rep(1,n1), p1_observ)+sum(summands)

        if (alpha >= type1 & beta >= type2) {
          stopp<-2
          r_cand[r1+1] <- r
          n2_cand[r1+1] <- n2
        }
      }
        n2 <- n2+1
    }
  }
  EN0_cand <- n1+n2_cand*(1-GenBinomApps::pgbinom(0:r1_up, rep(1,n1), p0_observ))

  EN0 <- min(EN0_cand)
  r1 <- which.min(EN0_cand)-1
  PET <- GenBinomApps::pgbinom(r1, rep(1,n1), p0_observ)

  n2 <- n2_cand[which.min(EN0_cand)]
  r <- r_cand[which.min(EN0_cand)]

  sumup <- min(r,n1)
  sumindx <- ((r1+1):sumup)
  summands <- pbinom(r-sumindx, n2, mean0) * GenBinomApps::dgbinom(sumindx, size = rep(1,n1), prob = p0_observ)
  type1 <- 1-(GenBinomApps::pgbinom(r1, rep(1,n1), p0_observ) + sum(summands))

  summands <- pbinom(r-sumindx, n2, mean1)*GenBinomApps::dgbinom(sumindx, rep(1,n1) , p1_observ)
  type2 <- GenBinomApps::pgbinom(r1, rep(1,n1), p1_observ)+sum(summands)


  interim <- c(r1, n1, r, n2, EN0, PET, type1, type2)
  names(interim) <- c("r1*" , "n1",  "r*", "n2*","EN*(p0)","PET*(p0)","type1_est","type2_est")

  out <- start

  out$p0_stageone <- p0_observ
  out$p1_stageone <- p1_observ
  out$des_interim <- interim

  return(out)

}

#' Adjust a design based on historical controls at the end of the study using the covariate data of the patients accrued in stage 2.
#' @export
#' @param interim An design based on historical controls and adjusted at interim as returned by \code{hist_interim()}.
#' @param stagetwo_data A dataframe containing the relevant covariate data of the patients accrued in stage 2.
#' @return A list returning the arguments of the function and the final design of the trial.
#' @examples
#' X <- abs(rnorm(1000, 0, 1))
#' Y <- rbinom(1000, 1, 1-exp(-X))
#' mydata <- data.frame("X" = X, "Y" = Y)
#' start <- hist_start(mydata, Y~X, c1 = 2)
#' n1 <- start$des_start[2]
#' X1 <- abs(rnorm(n1, 0, 1))
#' dataone <- data.frame("X" = X1)
#' interim <- hist_interim(start, dataone)
#' n2 <- interim$des_interim[4]
#' X2 <- abs(rnorm(n2, 0, 1))
#' datatwo <- data.frame("X" = X2)
#' hist_end(interim, datatwo)


hist_end <- function(interim, stagetwo_data) {

  r1 <- interim$des_interim[1]
  n1 <- interim$des_interim[2]
  n2 <- interim$des_interim[4]

  p0_stageone <- interim$p0_stageone
  p1_stageone <- interim$p1_stageone

  p0_stagetwo <- predict(interim$modelfit, newdata=stagetwo_data, type="response")
  p1_stagetwo <- sapply(p0_stagetwo, interim$phic1)

  alpha <- interim$alpha
  beta <- interim$beta

  type1 <- 1
  r <- r1
  while (type1 > alpha)
  {
    r <- r+1
    sumup <- min(r,n1)
    sumindx <- ((r1 + 1):sumup)
    summands <- GenBinomApps::pgbinom(r-sumindx, rep(1,n2), p0_stagetwo) * GenBinomApps::dgbinom(sumindx, rep(1,n1), p0_stageone)
    type1 <- 1-(GenBinomApps::pgbinom(r1, rep(1,n1), p0_stageone)+sum(summands))
  }

  sumup <- min(r,n1)
  sumindx <- ((r1 + 1):sumup)
  summands <- GenBinomApps::pgbinom(r-sumindx, rep(1,n2), p1_stagetwo) * GenBinomApps::dgbinom(sumindx, rep(1,n1), p1_stageone)
  type2 <- GenBinomApps::pgbinom(r1, rep(1,n1), p1_stageone)+sum(summands)


  des_end <- interim$des_interim

  des_end[3] <- r
  names(des_end)[3] <- "r**"
  des_end[c(7,8)] <- c(type1, type2)

  out <- interim

  out$p0_stagetwo <- p0_stagetwo
  out$p1_stagetwo <- p1_stagetwo

  out$des_end <- des_end

  return(out)


}






