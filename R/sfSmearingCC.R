#' @title External cross-contamination of fish fillets during smearing with salt, sugar or spices
#'
#' @description
#' The function [sfSmearingCC()] simulates the potential external contamination of fish fillets during the process of
#' dry-salting when producing smoked salmon, or during smearing the fish fillets with salt, sugar and spices when producing
#' gravad or any macerated fish. The algorithm evaluates the cross-contamination event at fish fillet level (and not at lot level),
#' assuming that every fish fillet has the same probability `pccSmearing` of being contaminated when in contact
#' with contaminated environmental elements during the process of smearing with salt, sugar and/or spices.
#'
#' @param data a list of
#'    \describe{
#'      \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units containing the numbers of \emph{L. monocytogenes} on fish fillets;}
#'      \item{`P`}{Mean prevalence of contaminated lots (scalar);}
#'      \item{`ProbUnitPos`}{Probability of individual lots being contaminated (vector).}
#'    }
#' @param pccSmearing Probability that cross-contamination with \emph{L. monocytogenes} occurs during dry salting,
#'  or smearing the fillets with sugar/spices, through tables or other surfaces (scalar).
#' @param trSmearingMean Mean parameter of the normal distribution representing the variability in the log 10
#'  of the transfer coefficient of \emph{L. monocytogenes} cells from surfaces to fish fillets (scalar or vector).
#' @param trSmearingSd Standard deviation parameter of the normal distribution representing the variability in the
#'  log 10 of the transfer coefficient of \emph{L. monocytogenes} from surfaces to fish fillets (scalar or vector).
#' @param nSurface (`CFU`) Numbers of \emph{L. monocytogenes} cells on environmental elements in contact with fish
#'  fillets while dry salting (or while smearing the fillets with sugar/spices) (scalar or vector).
#' @return A list of three elements:
#'     \describe{
#'        \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units containing the numbers of \emph{L. monocytogenes}
#'         in (dry-salted) fish;}
#'        \item{`ProbUnitPos`}{Probability of individual lots being contaminated after smearing with ingredients (vector);}
#'        \item{`P`}{Mean prevalence of contaminated lots after smearing with ingredients (scalar).}
#'            }
#'
#' @author Ursula Gonzales-Barron \email{ubarron@ipb.pt} and Regis Pouillot \email{rpouillot.work@gmail.com}
#'
#' @keywords salt sugar spices smearing external surface contamination transfer

#' @references
#'
#' \insertRef{Hoelzer2012}{qraLm}
#' \insertRef{extraDistr}{qraLm}
#' \insertRef{stats}{qraLm}
#' \insertRef{iRisk2021}{qraLm}
#'
#' @importFrom extraDistr rtnorm
#' @importFrom stats rbinom
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @note The suggested parameters `trSmearingMean=-0.29` and `trSmearingSd=0.31` for the normal distribution about the
#' variability in the log 10 of the transfer coefficient of \emph{L. monocytogenes} were taken from \insertCite{Hoelzer2012;textual}{qraLm}, to
#' represent cross-contamination from board to meat. The values of `pccSmearing` and `nSurface` must be defined by the
#' user and/or assessed in scenarios.
#'
#' @examples
#' 
#' dat <- Lot2LotGen(
#'                   nLots = 50,
#'                   sizeLot = 100,
#'                   unitSize = 500,
#'                   betaAlpha = 0.5112,
#'                   betaBeta = 9.959,
#'                   C0MeanLog = 1.023,
#'                   C0SdLog = 0.3267,
#'                   propVarInter = 0.7
#'                   )
#' Nf <- sfSmearingCC(
#'                    dat,
#'                    pccSmearing = 0.05,
#'                    trSmearingMean = -0.29,
#'                    trSmearingSd = 0.31,
#'                    nSurface = 200
#'                    )
#' str(Nf)                    
#' Nf$P
#' mean(Nf$ProbUnitPos)
#' hist(Nf$N)
#'
sfSmearingCC <- function(data = list(),
                         pccSmearing,
                         trSmearingMean,
                         trSmearingSd,
                         nSurface) {
  Nt_batch <- data$N
  ProbUnitPos <- data$ProbUnitPos
  
  # Get the simulation dimension
  # if(missing(nLots)) nLots <- data$nLots #test if nLots was defined
  # if(is.null(nLots)) warning("Add 'nLots=#' to function arguments") #test again if nLots is defined
  # 
  # if(missing(sizeLot)) sizeLot <- data$sizeLot #test if sizeLot was defined
  # if(is.null(sizeLot)) warning("Add 'sizeLot=#' to function arguments") #test again if sizeLot is defined
  nLots <- nrow(Nt_batch)
  sizeLot <- ncol(Nt_batch)
  
  ##############################################################################
  # Lots of contaminated fish (coming from the original matrix N)
  ##############################################################################
  
  pcc_flag <- matrix(stats::rbinom(nLots * sizeLot, size = 1, prob = pccSmearing), nrow = nLots, ncol = sizeLot)
  p_TR1 <- matrix((10^extraDistr::rtnorm(nLots * sizeLot, mean = trSmearingMean, sd = trSmearingSd, b = 0)),
                  nrow = nLots, ncol = sizeLot
  )
  N_TR1 <- pcc_flag * stats::rbinom(nLots * sizeLot, size = nSurface, prob = p_TR1)
  N_total1 <- Nt_batch + N_TR1
  
  ################################################################################
  # Lots of non-contaminated fish (lots added from "1-P")
  ################################################################################
  # Do more than needed for a better resampling
  ad_batches <- max(round(nLots * (1 / data$P - 1)), nLots)
  if (length(pccSmearing) != 1) stop("Need to update sfSmearingCC function")
  pcc_flag2 <- matrix(stats::rbinom(ad_batches * sizeLot, size = 1, prob = pccSmearing),
                      nrow = ad_batches, ncol = sizeLot
  )
  p_TR2 <- matrix((10^extraDistr::rtnorm(ad_batches * sizeLot, mean = trSmearingMean, sd = trSmearingSd, b = 0)),
                  nrow = ad_batches, ncol = sizeLot
  )
  N_total2 <- pcc_flag2 * stats::rbinom(ad_batches * sizeLot, size = nSurface, prob = p_TR2)
  
  # Probability batches remain non-contaminated (rows with zero counts)
  zero_lines <- rowSums(N_total2) == 0
  prob_clean_batches <- mean(zero_lines)
  # Keeping only batches that ended up contaminated (from the clean fraction)
  N_total2_c <- N_total2[!zero_lines, , drop=FALSE]
  # Check the number of remaining positive lots
  N_newly_stock <- nrow(N_total2_c)
  
  ###########################################################################
  # calculating probabilities of each contamination fraction and resizing
  ###########################################################################
  p_PosFraction_CCevent <- data$P
  # Modified RP (was an error identified by UGB on 18 Jan 2024)
  p_NegFraction_CCevent <- (1 - data$P) * (1 - prob_clean_batches)
  overall_p_batches <- p_PosFraction_CCevent + p_NegFraction_CCevent
  
  # re-arranging matrix of contaminated batches, and updating ProbUnitPos
  origin <- base::sample(c("Previously", "Newly"),
                         size = nLots,
                         prob = c(p_PosFraction_CCevent, p_NegFraction_CCevent), # Note: function will standardize to 1
                         replace = TRUE
  )
  
  N_out <- matrix(NA, nrow = nLots, ncol = sizeLot)
  
  # Lines of previously Positive lots
  N_out[origin == "Previously", ] <- N_total1[origin == "Previously", ]
  
  # Lines of newly Positive lots
  N_newly <- sum(origin == "Newly")
  if (N_newly != 0) {
    if (N_newly_stock < N_newly) {
      warning("have to reshuffle in cc_smearing")
      sample_index2 <- base::sample(1:N_newly_stock, N_newly, replace = TRUE)
    } else {
      sample_index2 <- 1:N_newly
    }
    N_out[origin == "Newly", ] <- N_total2_c[sample_index2, ]
  }
  
  
  # output
  N <- N_out
  # Modified RP (was an error identified by UGB on 18 Jan 2024)
  ProbUnitPos <- ProbUnitPos + (1 - ProbUnitPos) * (1 - prob_clean_batches)
  P <- overall_p_batches
  
  
  #lotMeans <- rowMeans(N / data$unitSize, na.rm = TRUE)
  unitsCounts <- c((ProbUnitPos/mean(ProbUnitPos)) * (N / data$unitSize))
  
  #data$lotMeans <- lotMeans
  data$unitsCounts <- unitsCounts
  data$N <- N
  # Modified RP (was an error identified by UGB on 18 Jan 2024)
  data$ProbUnitPos <- ProbUnitPos
  data$P <- P
  return(data)
}
