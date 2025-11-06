#' @title Cross-contamination during freezing/packaging and partitioning into packed units
#'
#' @description
#' The [fvPartitioningCC()] function simulates the potential re-contamination of vegetables post-blanching, when in direct
#' contact with contaminated surfaces of conveyors, freezer or packaging machine, followed by the partitioning of the bulk
#'  of frozen vegetables into packed units.
#' The cross-contamination algorithm accounts for four possible scenarios:
#' \enumerate{
#'    \item cross-contamination occurring in lots already contaminated;
#'    \item re-contamination occurring in lots that were not contaminated;
#'    \item no cross-contamination occurring in lots already contaminated; and
#'    \item no cross-contamination occurring in lots that were not contaminated.
#'    }
#' Probabilities of occurrence of every event are computed. The partitioning algorithm randomly distributes the total numbers of cells from a
#' contaminated lot into packed units. The dispersion factor `b`, which is a parameter of the beta distribution, indicates the
#' extent of cell clustering in the bulk of frozen vegetables in a lot, and ultimately the heterogeneity in the number of cells
#' distributed among pack units.
#'
#' @param data a list of:
#'   \describe{
#'      \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units containing the numbers of \emph{L. monocytogenes} in blanched
#'            units or portions from contaminated lots;}
#'      \item{`P`}{Mean prevalence of contaminated lots (scalar);}
#'      \item{`ProbUnitPos`}{Probability of individual lots being contaminated (vector).}
#'         }
#' @param nLots Number of lots sampled or size of the Monte Carlo simulation (scalar).
#' @param sizeLot Number of units or portions produced in a lot (scalar).
#' @param probCC Probability of cross-contamination of vegetables from conveyors, freezer or packaging machine  (scalar).
#' @param trMean Mean parameter of the normal distribution representing the variability in the log 10 of the transfer coefficient of \emph{L. monocytogenes}
#' cells from surfaces to (frozen) vegetables  (scalar or vector).
#' @param trSd Standard deviation parameter of the normal distribution representing the variability in the log 10 of the transfer coefficient of
#'  \emph{L. monocytogenes} from surfaces to (frozen) vegetables (scalar or vector).
#' @param nEquip (`CFU`) Numbers of \emph{L. monocytogenes} cells on the surface of conveyors, freezer or packaging machine in contact
#' with the (frozen) vegetables  (scalar or vector).
#' @param bCCFV Dispersion factor of the beta distribution representing the degree of heterogeneity in the number of cells between pack units (scalar).
#'
#' @return A list of three elements:
#'     \describe{
#'        \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units representing the numbers of \emph{L. monocytogenes} per pack unit,
#'          from contaminated lots;}
#'        \item{`ProbUnitPos`}{Probability of individual lots being contaminated after packaging (a lot is considered contaminated if at
#'          least one pack unit is contaminated) (vector);}
#'        \item{`P`}{Mean prevalence of contaminated lots after packaging (scalar).}
#'        }
#'
#' @author Ursula Gonzales-Barron \email{ubarron@ipb.pt} and Regis Pouillot \email{rpouillot.work@gmail.com}
#'
#' @keywords cross-contamination re-contamination transfer
#'
#' @references
#' \insertRef{Truchado2021}{qraLm}
#' \insertRef{Hoelzer2012}{qraLm}
#' \insertRef{Nauta2005}{qraLm}
#' \insertRef{extraDistr}{qraLm}
#' \insertRef{mc2d}{qraLm}
#' \insertRef{Rdpack}{qraLm}
#' \insertRef{stats}{qraLm}
#' \insertRef{iRisk2021}{qraLm}
#'
#' @importFrom extraDistr rtnorm
#' @importFrom mc2d rmultinomial rdirichlet
#' @importFrom stats rbinom
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @note The value of \eqn{probCC=0.125} is taken from \insertCite{Truchado2021;textual}{qraLm}, who analysed the presence of \emph{L. monocytogenes} on
#' conveyor belts, freezing tunnels and packaging machines in frozen food processing plants. The value of beta=1 represents moderate
#' clustering of cells in the bulk of frozen vegetables from a lot \insertCite{Nauta2005;textual}{qraLm}. \insertCite{Hoelzer2012;textual}{qraLm}
#' established the log 10 of the transfer coefficient of \emph{L. monocytogenes} from stainless steel to vegetables as a normal
#' distribution with \eqn{TR\_mean=-0.44} and \eqn{TR\_sd=0.40}.
#'
#' @examples
#' dat <- Lot2LotGen(
#'   nLots = 50,
#'   sizeLot = 100,
#'   unitSize = 500,
#'   betaAlpha = 0.5112,
#'   betaBeta = 9.959,
#'   C0MeanLog = 1.023,
#'   C0SdLog = 0.3267,
#'   propVarInter = 0.7
#' )
#' Nf <- fvPartitioningCC(dat,
#'                        probCC = 0.125, trMean = -0.44,
#'                        trSd = 0.40, nEquip = 500, bCCFV = 1
#' )
#' hist(Nf$N) # displays histogram of cells in a pack from contaminated lots
fvPartitioningCC <- function(data = list(),
                             nLots = NULL,
                             sizeLot = NULL,
                             probCC,
                             trMean=-0.44,
                             trSd=0.40,
                             nEquip,
                             bCCFV) {
  # Shortcut
  if(probCC == 0 | sum(nEquip) == 0){
    return(data)
  }
  
  ifelse(exists("nLots", data) == TRUE, nLots <- data$nLots, nLots <- nrow(data$N))
  ifelse(exists("sizeLot", data) == TRUE, sizeLot <- data$sizeLot, sizeLot <- ncol(data$N))

  # nLots <- nrow(data$N)
  # sizeLot <- ncol(data$N)

  Nt_batch <- rowSums(data$N)
  # Error in case a non-contaminated lot is present in the matrix
  if (any(Nt_batch == 0)) stop("There should not be any non-contaminated lot in fvPartitioningCC")

  # calculated probabilities
  p_BulkPos_CCNeg <- data$ProbUnitPos * (1 - probCC)
  p_BulkPos_CCPos <- data$ProbUnitPos * probCC
  p_BulkNeg_CCPos <- (1 - data$ProbUnitPos) * probCC
  # Add random
  # 1 is p_BulkPos_CCNeg, 2 is p_BulkPos_CCPos, 3 is p_BulkNeg_CCPos
  matProb <- cbind(p_BulkPos_CCNeg, p_BulkPos_CCPos, p_BulkNeg_CCPos)
  n_Bulk <- sapply(1:nLots, function(x) sample(1:3, size = 1, prob = matProb[x, ]))

  # Level of cross-contamination
  p_TR1 <- 10^extraDistr::rtnorm(nLots, mean = trMean, sd = trSd, b = 0)
  N_TR1 <- extraDistr::rtbinom(nLots, nEquip, p_TR1, a = 0)

  Nt_batchPostCC <- rep(0, nLots)
  # Previously contaminated: take the previous values
  previously <- n_Bulk %in% c(1, 2)
  Nt_batchPostCC[previously] <- Nt_batch[previously]
  # New contamination: add to 0 or the previous values the cross contamination
  newly <- n_Bulk %in% c(2, 3)
  Nt_batchPostCC[newly] <- Nt_batchPostCC[newly] + N_TR1[newly]

  # Defining the partition function to generate cells in a pack from the bulk

  if (length(bCCFV == 1)) {
    p <- mc2d::rdirichlet(nLots, rep(bCCFV, sizeLot))
    # mc2d multinomial because p is a matrix
    N_out <- mc2d::rmultinomial(n = nLots, size = Nt_batchPostCC, prob = p)
  } else {
    (stop("Vector of bCCFV non implemented yet."))
  }
  
  ProbUnitPos <- 1 - (1 - data$ProbUnitPos) * (1 - probCC)
  N <- N_out
  P <- 1 - (1 - data$P) * (1 - probCC)
  lotMeans <- rowMeans(N / data$unitSize, na.rm = TRUE)
  unitsCounts <- c((ProbUnitPos/mean(ProbUnitPos)) * (N / data$unitSize))
  
  data$lotMeans <- lotMeans
  data$unitsCounts <- unitsCounts
  data$N <- N
  data$P <- P
  data$ProbUnitPos <- ProbUnitPos

  return(data)
}
