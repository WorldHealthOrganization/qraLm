#' @title Effect of a potential intervention step post-packaging
#'
#' @description
#' The function [fvReductionPostpack()] describes the effect of any potential intervention step that would
#' be applied on packed frozen vegetables in order to reduce the occurrence of \emph{L. monocytogenes}.
#' It is intended that the user proposes as inputs the minimum, mode and maximum of the log10 reduction for this intervention (MRA46).
#'
#' @param data a list of:
#' \describe{
#'           \item{N}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units representing the numbers of \emph{L. monocytogenes} per pack unit, before the intervention;}
#'           \item{ProbUnitPos}{Probability of individual lots being contaminated before the intervention (vector);}
#'           \item{P}{Mean prevalence of contaminated lots before the intervention (scalar).}
#'           }
#' @param nLots Number of lots sampled or size of the Monte Carlo simulation (scalar).
#' @param sizeLot Number of units or portions produced in a lot (scalar).
#' @param logIntervMin (log10) minimum value of the triangular distribution representing log10 reduction (scalar).
#' @param logIntervMode (log10) mode value of the triangular distribution representing log10 reduction (scalar).
#' @param logIntervMax (log10) maximum value of the triangular distribution representing log10 reduction (scalar).
#'
#' @return A list of four elements:
#'     \describe{
#'              \item{N}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units representing the numbers of \emph{L. monocytogenes}
#'               per pack unit after post-packaging intervention;}
#'              \item{ProbUnitPos}{Probability that the lot is contaminated after post-packaging intervention (a lot is considered contaminated if at least
#'               one pack remained contaminated after the intervention) (vector);}
#'              \item{P}{Mean prevalence of contaminated lots after post-packaging intervention (scalar);}
#'              \item{pSurvReduction}{Probability of a microbial cell to survive the post-packaging intervention (vector).}
#'              }
#'
#' @author Laurent Guillier
#'
#' @keywords inactivation reduction
#'
#' @references
#' \insertRef{stats}{qraLm}
#'
#' \insertRef{mc2d}{qraLm}
#'
#' @importFrom mc2d rtriang
#' @importFrom stats rbinom
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
#' library(mc2d)
#' dat <- list(
#'   N = matrix(10^rnorm(1000 * 500, 1, 0), ncol = 1000, nrow = 500), P = 0.10,
#'   ProbUnitPos = rep(0.5, 500)
#' )
#' out <- fvReductionPostpack(dat,
#'   logIntervMin = 0,
#'   logIntervMode = 1,
#'   logIntervMax = 3
#' )
#' hist(out$N)
fvReductionPostpack <- function(data = list(),
                                nLots = NULL,
                                sizeLot = NULL,
                                logIntervMin,
                                logIntervMode,
                                logIntervMax) {
  Nbefore <- data$N
  ifelse(exists("nLots", data) == TRUE, nLots <- data$nLots, nLots <- nrow(Nbefore))
  ifelse(exists("sizeLot", data) == TRUE, sizeLot <- data$sizeLot, sizeLot <- ncol(Nbefore))

  #  nLots <- nrow(Nbefore)
  #  sizeLot <- ncol(Nbefore)


  log10_Intervention <- mc2d::rtriang(
    nLots,
    logIntervMin, logIntervMode, logIntervMax
  )
  # Intervention effect is variable from batch to batch
  pSurvive <- 10^(-log10_Intervention)

  # Trick to get only positive lots
  # apply the log-reduction to the sum and redistribute, proportionally to the original distribution
  # Watch out: overall equivalent, but not serving to serving (may lead to an increase for a given portion)
  # but not a problem here
  sumN <- rowSums(Nbefore)
  sumN1 <- extraDistr::rtbinom(n = nLots, size = sumN, prob = pSurvive, a = 0)
  # Protect
  sumN1 <- pmax(1, sumN)
  N1 <- mc2d::rmultinomial(n = nLots, size = sumN1, prob = Nbefore)

  # Evaluate the probabilities of at least one survive
  atLeastOneSurvive <- 1 - (1 - pSurvive)^sumN

  ProbUnitPos <- data$ProbUnitPos * atLeastOneSurvive

  P1 <- data$P * mean(atLeastOneSurvive)

  data$N <- N1
  data$P <- P1
  data$ProbUnitPos <- P1
  data$pSurvReduction <- pSurvive

  return(data)
}
