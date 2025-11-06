#' @title Cross-contamination during slicing of raw fish or processed fillets
#'
#' @description
#' The function [sfSlicer()] mimics the transfer of \emph{L. monocytogenes} during the slicing process of raw whole fish or
#' of processed fish fillets, using the compartmental model published in \insertCite{Hoelzer2012;textual}{qraLm}, which is
#' defined by two variability distributions: `a`, the transfer rate between slicer blade and product, and
#' `e`, the transfer rate from the original contamination to the slicing system. Transfer coefficients are sampled for every lot.
#'
#' @param data a list of
#'    \describe{
#'      \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units containing the numbers of \emph{L. monocytogenes}
#'       in/on the pre-slicing units;}
#'      \item{`workDone`}{Work done is optional (scalar, vector or matrix).}
#'    }
#' @param nLots Number of lots sampled or size of the Monte Carlo simulation (scalar).
#' @param sizeLot Number of units or portions produced in a lot (scalar).
#' @param wSlices Weight of a slice that can be obtained from one unit (scalar).
#' @param initSlicer (`CFU`) Initial contamination of the slicer (scalar or vector). Default is 0.
#' @param aParamLoc  Location parameter of the logistic distribution of `a`, the major parameter of the system (scalar).
#' @param aParamScale  Scale parameter of the logistic distribution of `a`, the major parameter of the system (scalar).
#' @param aParamMax  Maximum value used to truncate the logistic distribution of `a`, the major parameter of the  system (scalar).
#' @param eMean Mean of the log10 distribution of the `e` parameter (scalar)
#' @param eSd Standard deviation of the log10 distribution of the `e` parameter (scalar).
#'
#' @return A list of two elements:
#'    \describe{
#'      \item{`N`}{(`CFU`) A matrix of size nLots lots by sizeLot units containing the numbers of \emph{L. monocytogenes}
#'           in salted fish;}
#'      \item{`workDone`}{Updated work done (scalar, vector or matrix).}
#'            }

#' @author Regis Pouillot \email{rpouillot.work@gmail.com}
#'
#' @keywords contamination transfer coefficient
#'
#' @references
#'
#' \insertRef{Aarnisalo2007}{qraLm}
#'
#' \insertRef{extraDistr}{qraLm}
#'
#' \insertRef{Hoelzer2012}{qraLm}
#'
#' \insertRef{stats}{qraLm}
#'
#' @importFrom extraDistr rtnorm
#' @importFrom stats rbinom
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @note
#' The parameters of the distributions about the slicer transfer coefficients `a` and `e` are set to the values
#' determined in \insertCite{Hoelzer2012;textual}{qraLm}; although the user can alter or update any of them.
#' The number of columns of the input matrix `N` is increased by a factor of `nbSlices`.
#' The code is written so that if a workDone object is element of `data`, it is updated as needed.
#' When this function is used for slicing processed fillets, `nbSlices` can be assumed to be 40 (slices per fillet),
#' according to \insertCite{Aarnisalo2007;textual}{qraLm}.
#'
#' @examples
#' nLots <- 1000
#' sizeLot <- 1000
#' N <- matrix(1E5, nc = nLots, nr = sizeLot)
#' data <- list(N = N, unitSize = 1300)
#'
#' wSlices <- 20 # slice/fillet Aarnisalo 2007
#' Res <- sfSlicer(data, wSlices = wSlices, initSlicer = 0)
#' plot(log10(Res$N[1, ]), ylab = "log10 bacteria on each slices")
sfSlicer <- function(data = list(),
                     nLots = NULL,
                     sizeLot = NULL,
                     wSlices, # slice/fillet Aarnisalo 2007
                     initSlicer = 10,
                     aParamLoc = 0.07,
                     aParamScale = 0.03,
                     aParamMax = 0.5,
                     eMean = -2.12,
                     eSd = 0.85) {
  # Get the initial number of bacteria on the product
  initProd <- data$N
  # Need to do as.numeric (from integer) to avoid
  # errors such as "In M1 + N : NAs produced by integer overflow"
  initProd[] <- as.numeric(initProd)
  # Get the initial size of the simulation (lot and mc)

  ifelse(exists("nLots", data) == TRUE, nLots <- data$nLots, nLots <- nrow(initProd))
  ifelse(exists("sizeLot", data) == TRUE, sizeLot <- data$sizeLot, sizeLot <- ncol(initProd))

  #  sizeLot <- ncol(initProd)
  #  nLots <- nrow(initProd)
  # Check if there is a "workDone" object
  anyWorkDone <- !is.null(data$workDone)

  # Get the variability distribution of the slicer parameters
  # see Hoelzer2012
  min.a <- stats::plogis(0, location = aParamLoc, scale = aParamScale)
  max.a <- stats::plogis(aParamMax, location = aParamLoc, scale = aParamScale)
  a <- stats::qlogis(runif(nLots, min = min.a, max = max.a),
    location = aParamLoc,
    scale = aParamScale
  )
  # e is log10normal truncated on 0-1
  e <- 10^extraDistr::rtnorm(n = nLots, mean = eMean, sd = eSd, a = -Inf, b = 0)

  # Init
  cfu.on.slicer <- initSlicer

  # Initialize the matrices
  nbSlices <- floor(data$unitSize / wSlices)
  cfu.on.slices <- matrix(NA, nrow = nLots, ncol = sizeLot * nbSlices)
  if (anyWorkDone) {
    # Modify the workDone object to a comfortable array
    old.work.done <- matrix(data$workDone, nrow = nLots, ncol = sizeLot)
    # prepare the new work.done object (larger)
    new.work.done <- matrix(NA, nrow = nLots, ncol = sizeLot * nbSlices)
  }
  sliceIndex <- 1

  # Loop on the lot
  for (fillet in 1:sizeLot) {
    # Assume half the contamination is on one side of the fillet
    cfu.on.fillet <- cfu.on.fillet.Start <- rbinom(nLots, data$N[, fillet], prob = 0.5)

    for (slice in 1:(nbSlices - 1)) {
      # Lost
      cfu.on.fillet <- rbinom(nLots, cfu.on.fillet, 1 - e)
      # slicer -> slicer, fillet, slice via a multinomial
      # Stay on the slicer
      M1 <- rbinom(nLots, cfu.on.slicer, 1 - 2 * a)
      # Goes on fillet
      M2 <- rbinom(nLots, cfu.on.slicer - M1, 0.5)
      # Goes on slice
      M3 <- cfu.on.slicer - M1 - M2

      # fillet ->  slicer
      N <- rbinom(nLots, cfu.on.fillet, a)

      # Do the transfer
      # slices: remaining from the fillet and from the slicer
      cfu.on.slices[, sliceIndex] <- M3 + (cfu.on.fillet - N)
      # Slicer: Stay on the slicer + from the fillet
      cfu.on.slicer <- as.numeric(M1) + as.numeric(N)
      # Commented out: if(any(is.na(cfu.on.slicer))) browser()
      # fillet: from the slicer
      cfu.on.fillet <- M2

      # Make the workDone follow
      if (anyWorkDone) new.work.done[, sliceIndex] <- old.work.done[, fillet]

      sliceIndex <- sliceIndex + 1
    }

    # Last slice takes it all (remaining of the fillet + half the original contamination)
    cfu.on.slices[, sliceIndex] <- cfu.on.fillet + (data$N[, fillet] - cfu.on.fillet.Start)
    # Make the workDone follow
    if (anyWorkDone) new.work.done[, sliceIndex] <- old.work.done[, fillet]
    sliceIndex <- sliceIndex + 1
    # Goto next fillet
  }

  # Just in case of zeroes (but shouldn't happen)
  zeroes <- rowSums(cfu.on.slices) == 0
  if (any(zeroes)) {
    # warning(paste("Bacteria added in", sum(zeroes), "lines because no bacteria were left."))
    cfu.on.slices[zeroes, 1] <- 1
    Pi_0 <- sum(zeroes) / nLots
    data$P <- data$P * (1 - Pi_0)
  }

  N <- cfu.on.slices
  
  lotMeans <- rowMeans(N / wSlices, na.rm = TRUE)
  unitsCounts <- c((data$ProbUnitPos/mean(data$ProbUnitPos)) * (N / wSlices))
  
  data$lotMeans <- lotMeans
  data$unitsCounts <- unitsCounts
  data$N <- N
  data$sizeLot <- ncol(data$N)
  data$unitSize <- wSlices # update unitSize and add it to the dataset
  # If there is a workDone, repeat the workDone by the nb of slices.
  # Note:
  if (anyWorkDone) data$workDone <- new.work.done
  return(data)
}