#' Initializes a subspace stratified design before the start of the study.
#' @import stats
#' @export
#' @param p0_sub A vector, where the $i$-th entry corresponds to the response rate under the null for the $i$-th subtype.
#' @param p1_sub A vector, where the $i$-th entry corresponds to the response rate under the alternative for the $i$-th subtype.
#' @param distr_sub A vector, where the $i$-th entry corresponds to the prevalence of the $i$-th subtype in the population.
#' @param alpha Specified type I error of the trial.
#' @param beta Specified type II error of the trial.
#' @return A list returning the arguments of the function and the preliminary design for starting the stratified trial.
#' @examples
#' p0_sub <- c(0.1, 0.3, 0.5)
#' p1_sub <- c(0.3, 0.5, 0.7)
#' distr_sub <- c(1/3, 1/3, 1/3)
#' strat_start(p0_sub, p1_sub, distr_sub)

strat_start <- function(p0_sub, p1_sub, distr_sub, alpha = 0.05, beta = 0.2) {

  if (sum(distr_sub) != 1)
    stop("distr_sub is not a valid probability mass function")

  if (any(p0_sub<0) | any(p0_sub>1))
    stop("probabilities must be in the interval [0,1]")

  if (any(p1_sub<0) | any(p1_sub>1))
    stop("probabilities must be in the interval [0,1]")

  if (any(p1_sub-p0_sub<=0))
    stop("probabilites under alternative must be higher than under null")


  m <- length(distr_sub)

  if (length(p0_sub) != m | length(p1_sub) != m)
    stop("vector p0_sub, p1_sub and distr_sub must have the same length")

  mean0 <- sum(p0_sub * distr_sub)
  mean1 <- sum(p1_sub * distr_sub)


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
  type2 <- GenBinomApps::pgbinom(r1, n1, mean1)+sum(summands)

  out <- c(out, type1, type2)

  names(out)[c(7,8)] <- c("type1", "type2")

  return(list("p0_sub" = p0_sub, "p1_sub" = p1_sub, "distr_sub" = distr_sub, "alpha" = alpha, "beta" = beta, "des_start" = out))


}


#' Adjust a subspace stratified design at interim.
#' @export
#' @param start An initialized stratified design as returned by \code{strat_start()}.
#' @param sub_stageone The subtypes observed for the patients accrued in stage 1.
#' @return A list returning the arguments of the function and the preliminary design of a stratified trial adjusted at interim.
#' @examples
#' p0_sub <- c(0.1, 0.3, 0.5)
#' p1_sub <- c(0.3, 0.5, 0.7)
#' distr_sub <- c(1/3, 1/3, 1/3)
#' start <- strat_start(p0_sub, p1_sub, distr_sub)
#' n1 <- start$des_start[2]
#' subone <- sample(c(1,2,3), n1, TRUE)
#' strat_interim(start, subone)




strat_interim <- function (start, sub_stageone) {

  n1 <- start$des_start[2]



  if (length(sub_stageone) != n1)
    stop("length of observations must equal n1")

  alpha <- start$alpha
  beta <- start$beta
  p0_observ <- start$p0_sub[sub_stageone]
  p1_observ <- start$p1_sub[sub_stageone]
  mean0 <- sum(start$p0_sub * start$distr_sub)
  mean1 <- sum(start$p1_sub * start$distr_sub)

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
  names(interim) <- c("r1*" , "n1",  "r*", "n2*","EN*(p0)","PET*(p0)","type1","type2")

  out <- start

  out$sub_stageone <- sub_stageone
  out$des_interim <- interim

  return(out)

}


#' Adjust a subspace stratified design at the end of the study.
#' @export
#' @param interim A preliminary stratified design adjusted at interim as returned by \code{strat_interim()}.
#' @param sub_stagetwo The subtypes observed for the patients accrued in stage 2.
#' @return A list returning the arguments of the function and the final design of the stratified trial.
#' @examples
#' p0_sub <- c(0.1, 0.3, 0.5)
#' p1_sub <- c(0.3, 0.5, 0.7)
#' distr_sub <- c(1/3, 1/3, 1/3)
#' start <- strat_start(p0_sub, p1_sub, distr_sub)
#' n1 <- start$des_start[2]
#' subone <- sample(c(1,2,3), n1, TRUE)
#' interim <- strat_interim(start, subone)
#' n2 <- interim$des_interim[4]
#' subtwo <- sample(c(1,2,3), n2, TRUE)
#' strat_end(interim, subtwo)



strat_end <- function(interim, sub_stagetwo) {

  r1 <- interim$des_interim[1]
  n1 <- interim$des_interim[2]
  n2 <- interim$des_interim[4]

  p0_stageone <- interim$p0_sub[interim$sub_stageone]
  p1_stageone <- interim$p1_sub[interim$sub_stageone]

  p0_stagetwo <- interim$p0_sub[sub_stagetwo]
  p1_stagetwo <- interim$p1_sub[sub_stagetwo]

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

  out$sub_stagetwo <- sub_stagetwo
  out$des_end <- des_end

  return(out)


}






