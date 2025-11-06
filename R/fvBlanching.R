#' @title  Effect of blanching on \emph{L. monocytogenes} in processed vegetables
#'
#' @description
#' The function [fvBlanching()] describes the effect of blanching on the numbers of \emph{L. monocytogenes}
#' present in vegetables such as broccoli, mushroom, onions, peas and pepper. The function is based on the Bigelow
#' model, which describes the decimal reduction time (`D`) as a function of temperature, with parameters `z` and reference `D` (`Dref`)
#' at 70 \eqn{^\circ C}.
#' This function assumes that the extent of reduction is global for the aforementioned vegetables, and depends on the duration and
#' temperature of blanching. Whereas the algorithm considers the `z` value as fixed, `Dref` is assumed to be strain-specific,
#' and thus different between lots.
#'
#' @param data a list of:
#' \describe{
#'   \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units containing the numbers  of
#'    \emph{L. monocytogenes} in pre-blanching units or portions from contaminated lots;}
#'   \item{`ProbUnitPos`}{Probability that the lot is contaminated (a lot is considered contaminated
#'   if at least one unit or portion is contaminated) (vector);}
#'   \item{`P`}{Mean prevalence of contaminated lots (scalar).}
#' }
#' @param nLots Number of lots sampled or size of the Monte Carlo simulation (scalar).
#' @param tempBlanch (\eqn{^\circ C}) Temperature of blanching (scalar or vector).
#' @param timeBlanch (`min`) Duration of blanching (scalar or vector).
#' @param logDrefMean  Mean parameter of the normal distribution representing the variability in the log 10 of the
#'  reference `D` (`min`) at 70 \eqn{^\circ C} (scalar or vector).
#' @param logDrefSd Standard deviation of the normal distribution representing the variability in
#'  the log 10 of the reference `D` (min) at 70 \eqn{^\circ C} (scalar or vector).
#' @param zT fixed `z` value characterising the effect of temperature on `D` for \emph{L. monocytogenes} in vegetables (scalar).
#' @return A list of four elements:
#'     \describe{
#'              \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units representing the numbers of \emph{L. monocytogenes}
#'               per unit or portion after blanching;}
#'              \item{`ProbUnitPos`}{Probability that the lot is contaminated after blanching (a lot is considered contaminated if at least
#'               one unit or portion remained contaminated after blanching) (vector);}
#'              \item{`P`}{Mean prevalence of contaminated lots after blanching (scalar);}
#'              \item{`pSurviveBlanching`}{Probability of a microbial cell to survive the blanching process (vector).}
#'              }
#'
#' @author Laurent Guillier
#'
#' @keywords blanching cooling heat treatment
#'
#' @references
#' \insertRef{stats}{qraLm}
#'
#' \insertRef{Mazzotta2001}{qraLm}
#'
#' @importFrom stats rnorm rbinom
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @note If no more bacteria are in a lot, one cell will be set in one unit of that lot. The heat resistance parameters suggested for use
#' (\eqn{zT=6.06\ ^{\circ} C}, and \eqn{log10\ Dref\_mean=-1.78} and \eqn{log10Dref\_sd=0.251} at the reference temperature of 70 \eqn{^{\circ} C})
#'  were obtained from fitting a Bigelow model to survival data of \emph{L. monocytogenes} in broccoli, mushrooms, onions, peas and pepper
#' extracted from \insertCite{Mazzotta2001;textual}{qraLm}
#'
#' @examples
#' dat <- Lot2LotGen(
#'   nLots = 1000, sizeLot = 1000, unitSize = 500,
#'   betaAlpha = 0.5112, betaBeta = 9.959, C0MeanLog = 1.023, C0SdLog = 0.3267,
#'   propVarInter = .7
#' )
#' res <- fvBlanching(dat,
#'   tempBlanch = 85,
#'   timeBlanch = 1,
#'   logDrefMean = -1.78,
#'   logDrefSd = 0.251,
#'   zT = 6.06
#' )
#' res$P
#' hist(c(res$N))
#'
fvBlanching <- function(data = list(),
                        nLots = NULL,
                        tempBlanch = 83,
                        timeBlanch = 0.75,
                        logDrefMean = -1.78,
                        logDrefSd = 0.251,
                        zT = 6.06){
  ifelse(exists("nLots", data) == TRUE, nLots <- data$nLots, nLots <- nrow(data$N))
  #  if (exists("nLots", data)==TRUE) {
  #    nLots <- data$nLots
  #    } else {
  #     print("Add 'nLots=#' to function arguments")
  #   }
  #  nLots <- nrow(data$N)
  #  sizeLot <- ncol(data$N)
  N0 <- data$N

  # Draw heat inactivation values at batch level
  log10Dref <- stats::rnorm(
    nLots,
    logDrefMean,
    logDrefSd
  )
  D <- timeBlanch * (log10Dref - ((tempBlanch - 70) / zT))
  pSurvive <- 10^D

  # Trick to get only positive lots
  # apply the log-reduction to the sum and redistribute, proportionally to the original distribution
  # Watch out: overall equivalent, but not serving to serving (may lead to an increase for a given portion)
  # but not a problem here
  sumN <- rowSums(data$N)
  sumN1 <- extraDistr::rtbinom(n = nLots, size = sumN, prob = pSurvive, a = 0)
  # To avoid special cases with extremely low pSurvive
  # (check extraDistr::rtbinom(n= 10, size=10, prob = 1E-18, a=0))
  sumN1 <- pmax(1, sumN)

  N1 <- mc2d::rmultinomial(n = nLots, size = sumN1, prob = data$N)

  # Evaluate the probabilities of at least one survive
  atLeastOneSurvive <- 1 - (1 - pSurvive)^sumN

  ProbUnitPos <- data$ProbUnitPos * atLeastOneSurvive

  P1 <- data$P * mean(atLeastOneSurvive)

  lotMeans <- rowMeans(N1 / data$unitSize, na.rm = TRUE)
  unitsCounts <- c((ProbUnitPos/mean(ProbUnitPos)) * (N1/data$unitSize)) 
  
  data$lotMeans <- lotMeans
  data$unitsCounts <- unitsCounts
  data$nLots <- nLots
  data$N <- N1
  data$P <- P1
  data$ProbUnitPos <- ProbUnitPos
  data$pSurviveBlanching <- pSurvive
  return(data)
}
