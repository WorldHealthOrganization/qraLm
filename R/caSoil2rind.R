#' @title Route of contamination of cantaloupe from soil to rind
#'
#' @description
#' The function [caSoil2rind()] evaluates the contamination of cantaloupes through soil only.
#' It considers soil contamination characteristics such as prevalence `pSoil` and concentration `c_soil` as well as the
#' quantity of soil deposited on the cantaloupe rind. Furthermore, the algorithm allows for the use of foils to protect
#' cantaloupe from soil contamination, and the user must provide the probability of growing cantaloupes using such protective foils.
#' `pSoil` and `c_soil` have to be chosen by the user according to existing data; whereas `pSoil` is conditional to several risk factors.
#' The identified risk factors \insertCite{Background2022}{qraLm} are:
#'   \enumerate{
#'      \item irrigation in the previous days before harvesting; and
#'      \item the use of organic fertilizer.
#'      }
#' which affect `pSoil` as associated odds-ratios (`F_irrig_rain` and `fManure`). In this context, the user has to define the proportion of fields
#' practising irrigation prior to harvest or undergoing rain (`p_irrig_rain`), and the proportion of fields using organic amendments (`pManure`).
#'
#' @param nLots Size of the Monte Carlo simulation (scalar).
#' @param sizeLot Number of cantaloupes per cultivation lot or field (scalar).
#' @param pSoil Probability of contamination of soil (provided by the user, \eqn{default=0.089} according to \insertCite{Strawn2013a;textual}{qraLm})
#' @param fManure Odds-ratio estimate associated to use or organic amendment in soil (\eqn{default=7.0} according to \insertCite{Strawn2013b;textual}{qraLm})
#' @param pManure Proportion of fields using organic amendments (provided by user, \eqn{default=0.5})
#' @param fIrrigRaining Odds-ratio estimate associated to use of irrigation and raining events up to 2 days before harvest
#' (\eqn{default=25.0} according to \insertCite{Weller2015;textual}{qraLm}).
#' @param pIrrigRaining Proportion of fields that undergo irrigation or raining just previous harvest (provided by user, \eqn{default=0.1}).
#' @param cSoilLogMin (log10 CFU/g) Minimum value of the triangular distribution describing variability of concentration
#' (according to \insertCite{Dowe1997;textual}{qraLm} \eqn{default=-1\ log10\ CFU/g}).
#' @param cSoilLogMode (log10 CFU/g) Mode value of the triangular distribution describing variability of concentration
#' (according to \insertCite{Dowe1997;textual}{qraLm} \eqn{default=0.6\ log10\ CFU/g}).
#' @param cSoilLogMax (log10 CFU/g) Maximum value of the triangular distribution describing variability of concentration
#' (according to  \insertCite{Dowe1997;textual}{qraLm}: \eqn{default=1.48\ log10\ CFU/g}).
#' @param qSoilMin (`g`) Minimum value of the triangular distribution describing variability of quantity of soil deposited on cantaloupe (\eqn{default=0.05\ g}).
#' @param qSoilMode (`g`) Mode value of the triangular distribution describing variability of quantity of soil deposited on cantaloupe (\eqn{default=0.5\ g}).
#' @param qSoilMax (`g`) Maximum value of the triangular distribution describing variability of quantity of soil deposited on cantaloupe (\eqn{default=5\ g}).
#' @param pFoil Proportion of fields grown in foil (e.g. plastic mulch) (EKE: \eqn{default=0.5})
#' @param rFoil Reduction fraction of the quantity of soil transferred to rind (EKE: \eqn{default=0.9})
#'
#' @return A list of two elements:
#'     \describe{
#'              \item{`N`}{(`CFU`) A matrix of size `nLots` by `sizeLot` containing the numbers of \emph{L. monocytogenes} cells on cantaloupe
#'                 from the soil route;}
#'              \item{`P`}{Prevalence of field lots of cantaloupes contaminated with \emph{L. monocytogenes} from the soil route (scalar).}
#'              }
#' @author Laurent Guillier
#'
#' @keywords soil cross-contamination transfer
#'
#' @references
#' \insertRef{stats}{qraLm}
#'
#' \insertRef{mc2d}{qraLm}
#'
#' \insertRef{Background2022}{qraLm}
#'
#' \insertRef{Dowe1997}{qraLm}
#'
#' \insertRef{Strawn2013a}{qraLm}
#'
#' \insertRef{Strawn2013b}{qraLm}
#'
#' \insertRef{Weller2015}{qraLm}
#'
#' @importFrom stats rbinom runif rpois
#' @importFrom mc2d rtriang
#' @export
#'
#' @note The prevalence in soil, `pSoil`, must be provided by the user according to existing data. \insertCite{Background2022;textual}{qraLm}
#' lists various prevalence estimates in different regions of the world. A `default=0.089` \insertCite{Strawn2013a}{qraLm} is taken.
#' `pSoil` is conditional to several risk factors; yet this functions is based upon two risk factors:
#'   \enumerate{
#'     \item irrigation or rain occurring prior to harvesting, whose odds-ratio has been estimated at 25.0 when taking place 24 hours before
#'      harvesting \insertCite{Weller2015}{qraLm}; and
#'     \item application of organic fertilizer, whose odds-ratio has been estimated at 7.0 when manure is applied within 1 year \insertCite{Strawn2013b}{qraLm}.
#'     }
#' The distribution about the concentration of \emph{L. monocytogenes} in soil is represented as a triangular distribution,
#' using data from \insertCite{Dowe1997;textual}{qraLm}, who reported a mean of 4.0 MPN \emph{L. monocytogenes}/ g soil, with a 95%
#' confidence interval of \eqn{<1.0-28\ MPN/g}.
#'
#' @examples
#' # Considering all fields with foil: less cultivation lots contaminated than pSoil
#' caSoil2rind(nLots = 20, sizeLot = 10, pSoil = 0.089, pManure = 0.0, pIrrigRaining = 0, pFoil = 1.0)
#'
#' # Effect of manure (risk factor): more batches contaminated than pSoil
#' caSoil2rind(nLots = 20, sizeLot = 10, pSoil = 0.089, pManure = 0.5, pIrrigRaining = 0, pFoil = 0.0)
#'
caSoil2rind <- function(nLots,
                        sizeLot,
                        pSoil = 0.089,
                        fManure = 7,
                        pManure = 0.5,
                        fIrrigRaining = 25,
                        pIrrigRaining = 0.1,
                        cSoilLogMin = -1.0,
                        cSoilLogMode = 0.6,
                        cSoilLogMax = 1.48,
                        qSoilMin = 0.05,
                        qSoilMode = 0.5,
                        qSoilMax = 5,
                        pFoil = 0.5,
                        rFoil = 0.9) {
  # Quantity of soil by cantaloupe (cantaloupe to cantaloupe variability) (g)
  soil_gain <- mc2d::rtriang(nLots * sizeLot, qSoilMin, qSoilMode, qSoilMax)

  # For fields with foil protection the quantity of soil is lower
  # Was by melons:

  # soil_foil <- sample(nLots*sizeLot, nLots*sizeLot*pFoil, replace = FALSE)
  # soil_gain[soil_foil] <- soil_gain[soil_foil] * (1-rFoil)

  # I suggest (RP) by lot
  soil_foil <- stats::rbinom(n = nLots, size = 1, prob = pFoil)
  soil_gain <- soil_gain * (1 - soil_foil * rFoil)

  # Concentration of L. monocytogenes in contaminated soil variability (field to field)
  c_soil <- 10^mc2d::rtriang(n = nLots, min = cSoilLogMin, mode = cSoilLogMode, max = cSoilLogMax)

  # Quantity of L. monocytogenes on the whole cantaloupe rind
  N0 <- soil_gain * c_soil
  N0 <- stats::rpois(nLots * sizeLot, N0)
  N0 <- matrix(N0, nrow = nLots, ncol = sizeLot)

  # Resample zeroes
  # We keep it there
  zeroes <- rowSums(N0) == 0
  if (all(zeroes, TRUE)) {
    N0[, 1] <- 1
    zeroes <- rowSums(N0) == 0
  } else {
    # Trick to avoid bad behavior of sample
    N0[zeroes, ] <- N0[sample(rep(which(!zeroes), 2), sum(zeroes), replace = TRUE), ]
  }

  ## Prevalence of field batches is conditional to the prevalence of soil

  # pSoil
  p_soil_MC <- rep(pSoil, nLots)

  # pSoil considering manure for a certain proportion of fields
  p_soil_manure_MC <- p_soil_MC
  bmanure <- as.logical(stats::rbinom(nLots, 1, pManure))
  p_soil_manure_MC[bmanure] <- (fManure * pSoil) / (1 - pSoil + fManure * pSoil)

  # pSoil considering manure for a certain proportion of fields and raining/irrigation before harvest
  p_soil_irrig_raining_MC <- p_soil_manure_MC
  birrig_raining <- stats::rbinom(nLots, 1, pIrrigRaining)
  whichOnes <- as.logical(birrig_raining)
  p_soil_irrig_raining_MC[whichOnes] <- (fIrrigRaining * p_soil_manure_MC[whichOnes]) /
    (1 - p_soil_manure_MC[whichOnes] + fIrrigRaining * p_soil_manure_MC[whichOnes])

  p_soil_sim <- mean(p_soil_irrig_raining_MC)

  # Account that is some situation contamination of soil does not lead to cantaloupe contamination
  P0 <- p_soil_sim * mean(!zeroes)
  return(list(N = N0, P = P0))
}
