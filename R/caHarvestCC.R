#' @title Cross-contamination of cantaloupe rind by \emph{L. monocytogenes} during harvest
#'
#' @description
#' The function [caHarvestCC()] simulates cross-contamination that might occur at the moment of harvesting from elements such as conveyors,
#' crates or plastic surfaces. This function can be used for cantaloupes harvested in large farms (intended for RTE and
#' sold as whole cantaloupe in formal retail) and in small farms (intended for traditional markets). Parameters to be changed accordingly.
#'
#' @param data a list of
#' \describe{
#'       \item{`N`}{(`CFU`) A matrix of size: `nLots`  by `sizeLot` containing the numbers of \emph{L. monocytogenes} on the rind
#'        of cantaloupes from contaminated cultivation lots}
#'       \item{`P`}{Prevalence of contaminated cultivation lots (scalar)}
#'       }
#' @param probCCH Probability of cross-contamination during harvest
#' @param trMean Mean of the normal distribution describing the variability in the transfer coefficient from the
#'  surface of elements to rind (scalar)
#' @param trSd Standard deviation of the normal distribution describing the variability in the transfer coefficient
#'  from the surface of elements to rind (scalar)
#' @param nPlas (`CFU`) Numbers of \emph{L. monocytogenes} on food contact surfaces (conveyors or crates at harvest)
#'  touching the cantaloupes (scalar)
#' @param nLots Number of cultivation lots
#' @param sizeLot Number of cantaloupes cultivated and harvested in a lot

#' @return A list of two elements:
#'     \describe{
#'        \item{`N`}{(CFU) A matrix of size `nLots`  by `sizeLot` of \emph{L. monocytogenes} numbers on the rind of cantaloupes
#'         from contaminated harvested lots}
#'        \item{`P`}{Prevalence of harvested contaminated lots}
#'             }
#'
#' @author Ursula Gonzales-Barron \email{ubarron@ipb.pt} and Regis Pouillot \email{rpouillot.work@gmail.com}
#'
#' @keywords surface contamination conveyors crates
#'
#' @references
#' \insertRef{stats}{qraLm}
#'
#' @importFrom stats rbinom
#' @importFrom extraDistr rtnorm
#' @import Rdpack
#' @export
#'
#' @note The algorithm evaluates the cross-contamination event at cantaloupe level (and not at cultivation lot level), assuming that
#' every cantaloupe has the same probability `probCCH` of being contaminated from conveyors, crates or other surfaces during harvesting.
#'
#' @examples
#' # Example of use in LARGE production farms
#' dat <- caPrimaryProduction(
#' nLots = 100,
#'   sizeLot = 100)
#' probCCH <- 0.005
#' LargeFarms <- caHarvestCC(dat,
#'                           probCCH = probCCH, 
#'                           trMean = -1.42,
#'                           trSd = 0.52, 
#'                           nPlas = 200
#'                           )
#' hist(LargeFarms$N)
# Example of use in SMALL production farms
#' probCCH <- 0.02
#' nPlas <- 500
#' SmallFarms <- caHarvestCC(dat,
#'                           probCCH = probCCH, trMean = -1.42,
#'                           trSd = 0.52, nPlas = 500,
#'                           nLots = 100, sizeLot = 100
#'                           )
#'                           hist(SmallFarms$N)
caHarvestCC <- function(data = list(),
                        probCCH,
                        trMean,
                        trSd,
                        nPlas,
                        nLots = NULL, # Unneeded
                        sizeLot = NULL) {
  ##############################################################################
  # Batches of contaminated melons pre-harvest (coming from the original matrix)
  ##############################################################################
  if (missing(nLots)) nLots <- data$nLots # test if nLots was defined
  if (is.null(nLots)) warning("Add 'nLots=#' to function arguments") # test again if nLots is defined
  
  if (missing(sizeLot)) sizeLot <- data$sizeLot # test if sizeLot was defined
  if (is.null(sizeLot)) warning("Add 'sizeLot=#' to function arguments") # test again if sizeLot is defined
  # if(missing(P)) P <- data$P #test if P was defined
  # if(is.null(P)) warning("Add 'P=#' to function arguments") #test again if P is defined
  
  Nt_batch <- data$N
  pcc_flag <- matrix(stats::rbinom(nLots * sizeLot, size = 1, prob = probCCH), nrow = nLots, ncol = sizeLot)
  p_TR1 <- matrix((10^extraDistr::rtnorm(nLots * sizeLot, mean = trMean, sd = trSd, b = 0)),
                  nrow = nLots, ncol = sizeLot
  )
  N_TR1 <- pcc_flag * stats::rbinom(nLots * sizeLot, size = nPlas, prob = p_TR1)
  N_total1 <- Nt_batch + N_TR1
  
  ################################################################################
  # Batches of non-contaminated melons pre-harvest (batches added from "1-P")
  ################################################################################
  # Do more than needed for a better resampling
  ad_batches <- max(round(nLots * (1 / data$P - 1)), nLots)
  if (length(probCCH) != 1) stop("Need to update caHarvestCC function")
  pcc_flag2 <- matrix(stats::rbinom(ad_batches * sizeLot, size = 1, prob = probCCH),
                      nrow = ad_batches, ncol = sizeLot
  )
  p_TR2 <- matrix((10^extraDistr::rtnorm(ad_batches * sizeLot, mean = trMean, sd = trSd, b = 0)),
                  nrow = ad_batches, ncol = sizeLot
  )
  N_total2 <- pcc_flag2 * stats::rbinom(ad_batches * sizeLot, size = nPlas, prob = p_TR2)
  
  # Probability batches remain non-contaminated (rows with zero counts)
  zero_lines <- rowSums(N_total2) == 0
  prob_clean_batches <- mean(zero_lines)
  # Keeping only batches that ended up contaminated (from the clean pre-harvest fraction)
  N_total2_c <- N_total2[!zero_lines, ,drop=FALSE]
  # Check the number of remaining positive lots
  N_newly_stock <- nrow(N_total2_c)
  
  ###########################################################################
  # calculating probabilities of each contamination fraction and resizing
  ###########################################################################
  p_PosFraction_CCevent <- data$P
  p_NegFraction_CCevent <- (1 - data$P) * (1 - prob_clean_batches)
  overall_p_batches <- p_PosFraction_CCevent + p_NegFraction_CCevent
  
  # re-sizing matrix of contaminated batches
  origin <- sample(c("Previously", "Newly"),
                   size = nLots,
                   prob = c(p_PosFraction_CCevent, p_NegFraction_CCevent), # Note: function will standardize to 1
                   replace = TRUE
  )
  
  N_out <- matrix(NA, nrow = nLots, ncol = sizeLot)
  
  # Lines of previously Positive lots
  N_out[origin == "Previously", ] <- N_total1[origin == "Previously", ]
  
  N_newly <- sum(origin == "Newly")
  if(N_newly != 0){
    if (N_newly_stock < N_newly) {
      warning("have to reshuffle in cc_harvest")
      sample_index2 <- sample(1:N_newly_stock, N_newly, replace = TRUE)
    } else {
      sample_index2 <- 1:N_newly
    }
    N_out[origin == "Newly", ] <- N_total2_c[sample_index2, ]
  }
  
  # hist(N_out) #displays histogram of N (cells)
  N <- N_out
  P <- overall_p_batches
  
  lotMeans <- rowMeans(N / data$cantaWeight, na.rm = TRUE)
  unitsCounts <- c(N / data$cantaWeight)
  
  data$lotMeans <- lotMeans
  data$unitsCounts <- unitsCounts
  data$N <- N
  data$P <- P
  
  return(data)
}
