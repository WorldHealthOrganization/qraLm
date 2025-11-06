#' @title Generation of contaminated lots with variability in prevalence and concentration
#'
#' @description
#' The function [Lot2LotGen()] generates a matrix of contaminated lots from information on
#' the parameters of a beta distribution (`alpha`, `beta`) representing the between-lot variability in prevalence;
#' the parameters of a normal distribution (`mean`, `sd`) representing the total variability in the microbial concentration
#' in the contaminated units (`log10 CFU/g`); and the proportion of such total variance assigned to between-lot variance.
#' Every row of the matrix corresponds to a production lot, broken down in units or portions, which are the columns of the matrix.
#'
#' @param nLots Number of lots sampled or size of the Monte Carlo simulation (scalar).
#' @param sizeLot Number of units or portions produced in a lot (scalar).
#' @param unitSize Weight (`g`) of single unit or portion from a lot (scalar).
#' @param betaAlpha Parameter `alpha` of the beta distribution representing the
#'   between-batch variability in prevalence (scalar).
#' @param betaBeta Parameter `beta` of the beta distribution representing the
#'   between-batch variability in prevalence (scalar).
#' @param C0MeanLog Mean parameter (`log10 CFU/g`) of the normal distribution representing
#'   the total variability in the microbial concentration of contaminated units (scalar or vector).
#' @param C0SdLog Standard deviation parameter (`log10 CFU/g`) of the normal distribution
#'   representing the total variability in the microbial concentration of contaminated units (scalar or vector).
#' @param propVarInter proportion of the total variance (\eqn{C0\_sdLog10^2}) attributed
#'   to between-lot variance. The remaining (`1-propVarInter`) is the proportion of within-lot variance.
#' @param if Poisson is `TRUE`, a Poisson distribution is used with parameter lambda = -log(1 - prob) / 25, with prob ~ Beta(alpha, beta).
#' Note that, in this case `C0MeanLog` and `C0SdLog` are not used. 
#' @param ... Other options used to control [LotGen()]
#'
#' @return A list of four elements:
#'    \describe{
#'      \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units representing the microbial numbers
#'      in the units or portions from contaminated lots;}
#'      \item{`ProbUnitPos`}{Probability that the lot is contaminated (a lot is considered contaminated if at
#'      least one unit or portion is contaminated) (vector);}
#'      \item{`P`}{Mean prevalence of contaminated lots (scalar);}
#'      \item{`betaGen`}{True prevalence of contaminated units in the lots (sampled from the beta distribution) (vector).}
#'      }
#' @author Vasco Cadavez \email{vcadavez@ipb.pt} and Regis Pouillot \email{rpouillot.work@gmail.com}
#'
#' @keywords between-batch batches contamination numbers beta-binomial
#'
#' @family data generation
#'
#' @references
#' \insertRef{stats}{qraLm}
#' \insertRef{Background2022}{qraLm}
#' \insertRef{Jeya2011}{qraLm}
#' \insertRef{Magdovitz2021}{qraLm}
#' \insertRef{Kuan2017}{qraLm}
#'
#' @importFrom stats rnorm rbeta rbinom dbinom
#' @importFrom extraDistr dbbinom
#' @importFrom matrixStats rowSums2
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @note This function requires that the variability in the prevalence of \emph{L. monocytogenes} between lots is
#' represented by a beta distribution.
#' To determine the parameters of such a beta distribution, a beta-binomial distribution was fitted to survey data
#' (compiled in Table 3.21 of \insertCite{Background2022;textual}{qraLm}, consisting of positive samples and total sample size.
#' The beta-binomial was fitted through a Markov-Chain Monte Carlo approach, with `gamma(0.01, 0.01)` as priors for both parameters,
#' and assuming that every data source (a published paper) represents a different lot. The parameters of the beta distribution
#' were determined at: `betaAlpha=0.5112` and `betaBeta=9.959`.
#' The mean (`C0MeanLog=1.023` log10 CFU/g) and standard deviation (`C0SdLog=0.3267` log10\ CFU/g) of the normal distribution
#' representing the overall variability in \emph{L. monocytogenes} were obtained using data from \insertCite{Jeya2011;textual}{qraLm}
#' (3.0 MPN/g, 3.6 MPN/g in freshly harvested vegetables), \insertCite{Magdovitz2021;textual}{qraLm} (21, 36 and 75 MPN/g in corn and peas
#' arriving at the frozen food facility), and \insertCite{Kuan2017;textual}{qraLm} (3.0 MPN/g in bulk carrots at retail).
#' The algorithm decomposes at the log10 level the total variance in concentration into between-lot variance and within-lot
#' variance according to the proportion `propVarInter`.
#'
#' @examples
#'
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
#' 
#' print(dat$P)
#' summary(c(dat$N))
#'
Lot2LotGen <- function(nLots,
                       sizeLot,
                       unitSize,
                       betaAlpha,
                       betaBeta,
                       C0MeanLog,
                       C0SdLog,
                       propVarInter,
                       Poisson = FALSE,
                       ...){
  
  Lot2LotGenParameters <- c(as.list(environment()), list(...))
  
  # Prevalence (beta distribution)
  prob <- rbeta(nLots, betaAlpha, betaBeta)
  # Prob that the lot is not contaminated at all
  prob0 <- stats::dbinom(0, sizeLot, prob)

  if (!Poisson) { # Nb of positive in the lot (binomial truncated at zero)
    # Double protect the binomial because some prob can be extremely low
    nbPos <- extraDistr::rtbinom(n = nLots, size = sizeLot, prob = prob, a = 0)
    if (any(nbPos == 0)) {
      warning("Some values of intralot prevalence are extremely low, please check if OK")
      nbPos <- pmax(nbPos, 1)
    }
    # Build the matrix of 0/1
    posOrNeg <- t(sapply(nbPos, function(x) sample(c(rep(1, x), rep(0, sizeLot - x)))))


    # Concentration (log10 CFU/g)
    # Inter Lot (part of the variance)
    CMeanLot <- stats::rnorm(nLots,
      mean = C0MeanLog,
      sd = sqrt(propVarInter * C0SdLog^2)
    )
    # Intra lot. Truncated to have at least one bacteria per contaminated unit
    Cest <- matrix(
      extraDistr::rtnorm(nLots * sizeLot,
        mean = CMeanLot,
        sd = sqrt((1 - propVarInter) * C0SdLog^2),
        a = log10(0.501 / unitSize)
      ),
      nrow = nLots, ncol = sizeLot
    )

    # Then multiply by prev to fill only the contaminated units
    Ncounts <- round(posOrNeg * 10^Cest * unitSize)
  } else {
    # protect against log(0) and log(1)
    lambda <- -log(1 - pmax(pmin(prob, 1-10^-8), 10^-8)) / 25
    # Double protect the binomial because some prob can be extremely low
    NcountsLot <- extraDistr::rtpois(n = nLots, lambda = lambda * unitSize * sizeLot, a = 0)
    if (any(NcountsLot == 0)) {
      warning("Some values of intralot prevalence are extremely low, please check if OK")
      NcountsLot <- pmax(NcountsLot, 1)
    }
    Ncounts <- mc2d::rmultinomial(n = nLots, size = NcountsLot, prob = rep(1, sizeLot))
  }

  # So we just consider the beta binomial (because the lognormals are only for
  # contaminated ones)
  P0 <- 1 - extraDistr::dbbinom(x = 0, size = sizeLot, alpha = betaAlpha, beta = betaBeta)
  
  ProbUnitPos = 1 - prob0
  N <- Ncounts
  
  lotMeans <- rowMeans(N / unitSize, na.rm = TRUE)
  unitsCounts <- c((ProbUnitPos/mean(ProbUnitPos)) * (N/unitSize))
  
  data <- list(
               Lot2LotGenParameters=Lot2LotGenParameters,
               lotMeans = lotMeans,
               unitsCounts = unitsCounts,
               N = N,
               ProbUnitPos = ProbUnitPos,
               P = P0,
               betaGen = prob,
               nLots = nLots,
               sizeLot = sizeLot,
               unitSize = unitSize
               )
  class(data) <- "qraLm"

  return(data)
}
