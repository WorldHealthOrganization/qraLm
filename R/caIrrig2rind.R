#' @title Route of contamination of cantaloupe from irrigation water to rind
#'
#' @description
#' The function [caIrrig2rind()] evaluates the contamination of cantaloupes through irrigation water only.
#' It considers water contamination characteristics: prevalence `pIrrig` and concentration `cIrrig`.
#' `pIrrig` and `cIrrig` have to be chosen by the user according to existing data of prevalence. `cIrrig` is
#' conditional to water sources contaminated with \emph{L. monocytogenes}.
#'
#' @param nLots Size of the Monte Carlo simulation (scalar).
#' @param sizeLot Number of cantaloupes per cultivation lot or field (scalar).
#' @param pIrrig Prevalence of contamination in irrigation water (provided by the user, \eqn{default=0.131} according to \insertCite{Raschle2021;textual}{qraLm}).
#' @param cIrrigLogMin (log CFU/L) Minimum value of the uniform distribution (provided by the user, \eqn{default=-1.52} according to \insertCite{Sharma2020;textual}{qraLm}).
#' @param cIrrigLogMax (log CFU/L)  Maximum value of the uniform distribution (provided by the user, \eqn{default=1.04} according to \insertCite{Sharma2020;textual}{qraLm}).
#' @param cantaWeight (g) Weight of a cantaloupe.
#' @param pWaterGainMin Minimum value of the fraction of water gain (ml) relative to the cantaloupe weight in g (\eqn{default=0.0}).
#' @param pWaterGainMax Maximum value of the fraction of water gain (ml) relative to the cantaloupe weight in g (\eqn{default=0.004} according to \insertCite{Richards2004;textual}{qraLm}).
#' @return
#'  A list of two elements:
#'     \describe{
#'              \item{N}{(`CFU`) A matrix of size `nLots` by `sizeLot` containing the numbers of \emph{L. monocytogenes} cells on cantaloupe from the irrigation water route;}
#'              \item{P}{Prevalence of field lots of cantaloupes contaminated with \emph{L. monocytogenes} from the irrigation water route.}
#'              }
#'
#' @author
#' Laurent Guillier
#'
#' @keywords
#' irrigation cross-contamination transfer
#'
#' @references
#' \insertRef{stats}{qraLm}
#'
#' \insertRef{mc2d}{qraLm}
#'
#' \insertRef{Background2022}{qraLm}
#'
#' \insertRef{Raschle2021}{qraLm}
#'
#' \insertRef{Richards2004}{qraLm}
#'
#' \insertRef{Sharma2020}{qraLm}
#'
#' @importFrom stats rbinom runif rpois
#'
#' @export
#'
#' @note
#' \insertCite{Background2022;textual}{qraLm} lists many estimates of \emph{L. monocytogenes} in water environments.
#' The estimate of \eqn{0.131} provided by \insertCite{Raschle2021;textual}{qraLm} has been chosen as default.
#' The distribution about the concentration of \emph{L. monocytogenes} in irrigation water is represented as a uniform distribution,
#' using the minimum and maximum values from \insertCite{Sharma2020;textual}{qraLm}, who reported \eqn{<0.03\ to\ 11} MPN \emph{L. monocytogenes}/ L water.
#' The algorithm assumes that the amount of water deposited on the cantaloupe rind after the last irrigation is as a percentage of the cantaloupe weight,
#' and is sampled from a uniform distribution of parameters between a minimum value of zero and a maximum value of \eqn{0.004},
#' taken from \insertCite{Richards2004;textual}{qraLm}.
#'
#' @examples
#' dat <- caPrimaryProduction(
#'   nLots = 100,
#'   sizeLot = 100,
#'   pSoil = 0.089,
#'   pManure = 0.1,
#'   pIrrigRaining = 0.05,
#'   pFoil = 0.5,
#'   rFoil = 0.75,
#'   pIrrig = 0.131,
#'   cantaWeight = 800
#' )
#' 
#' caIrrig2rind(
#'   nLots = 5,
#'   sizeLot = 10,
#'   pIrrig = 0.5,
#'   cIrrigLogMin = 1,
#'   cIrrigLogMax = 2
#' )
#'
caIrrig2rind <- function(nLots,
                         sizeLot,
                         pIrrig = 0.131,
                         cIrrigLogMin = -1.52,
                         cIrrigLogMax = 1.04,
                         cantaWeight = 1000,
                         pWaterGainMin = 0,
                         pWaterGainMax = 0.004) {
  
  # Quantity of water by cantaloupe (function of their weights) in L
  water_gain <- cantaWeight * stats::runif(n = nLots * sizeLot, min = pWaterGainMin, pWaterGainMax) / 1000

  # Concentration of L. monocytogenes in contaminated water variability (batch to batch)
  c_irrig <- 10^stats::runif(n = nLots, min = cIrrigLogMin, max = cIrrigLogMax)

  # Quantity of L. monocytogenes on the whole cantaloupe rind
  N0 <- water_gain * c_irrig
  N0 <- stats::rpois(nLots * sizeLot, N0)
  N0 <- matrix(N0, nrow = nLots, ncol = sizeLot)

  # Resample zeroes
  zeroes <- rowSums(N0) == 0
  if (all(zeroes, TRUE)) {
    N0[, 1] <- 1
    zeroes <- rowSums(N0) == 0
  } else {
    # Nice trick to avoid sample(10,2, replace=TRUE) vs sample(c(10,11),2, replace=TRUE)
    N0[zeroes, ] <- N0[sample(rep(which(!zeroes), 2), sum(zeroes), replace = TRUE), ]
  }


  # Prevalence of field batches is conditional to the prevalence of contamination in irrigated water
  P0 <- pIrrig * mean(!zeroes)
  return(list(N = N0, P = P0))
}
