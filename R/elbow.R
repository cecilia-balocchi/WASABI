#' Elbow method for the WASABI algorithm
#'
#' @description Run the WASABI algorithm for different number of particles (from 1 to \code{L_max})
#' to find the optimal value of $L$ using the elbow method.
#'
#' @usage elbow(cls.draw, L_max = 10, psm = NULL,
#'       multi.start = 1, ncores = 1,
#'       method.init = "topvi", add_topvi = FALSE, lb = TRUE,
#'       thin.init = NULL, part.init = NULL,
#'       method = "average",
#'       max.k = NULL, max.iter = 20, eps = 0.01,
#'       mini.batch = NULL, extra.iter = NULL,
#'       swap_countone = FALSE,
#'       suppress.comment = TRUE, seed = NULL)
#'
#'
#' @param cls.draw A matrix of the MCMC samples of partitions of $n$ data points.
#' @param L_max An integer specifying the largest value of $L$ to run the WASABI algorithm for.
#' @param psm The posterior similarity matrix obtained from MCMC samples of partitions stored in \code{cls.draw}.
#' @param multi.start Integer, the number of random initializations for the WASABI algorithm. If \code{multi.start > 1}, the algorithm is run multiple times with different random initializations and the best result is returned.
#' @param ncores Integer, the number of cores to use for parallel processing. If \code{ncores > 1}, the WASABI algorithm is run in parallel.
#' @param method.init Initialization method. Options are "topvi" (default), "average", "complete", "greedy", "salso", and "fixed".
#' @param add_topvi Boolean, default FALSE. If TRUE, the top VI partition is added to the initial particles.
#' @param lb Logical, if TRUE, the lower bound for the VI is used in methods "average" or "complete".
#' @param thin.init Integer, thinning factor for the MCMC samples used to initialize the particles. If NULL, defaults to 10.
#' @param part.init A matrix of size \code{L} x \code{n}, containing the initial particles. Needs to be provided when \code{method.init = "fixed"}.
#' @param method The method used to find the partition with minimum EVI (minVI partition). Options are "average", "complete", "greedy", and "salso".
#' @param max.k Integer, the maximum number of clusters considered in the WASABI approximation (for "average", "complete", "greedy"). If NULL, it is set to the minimum of \code{max(Ks.draw)+10} and \code{ceiling(sqrt(n))}.
#' @param max.iter Integer, the maximum number of iterations for the WASABI algorithm.
#' @param eps Numeric, the convergence threshold for the WASABI algorithm. The algorithm stops when the difference in Wasserstein distance between two consecutive iterations is less than \code{eps}.
#' @param mini.batch Integer, the size of the mini-batch used in the WASABI algorithm. If 0, the full batch is used.
#' @param extra.iter Integer, the number of additional iterations to run after the mini-batch optimization. If NULL, defaults to 1 if \code{mini.batch > 0}. Has to be greater than 0.
#' @param swap_countone Logical, if TRUE, the WASABI algorithm allows swapping of particles with only one sample assigned to them (outlier-check step).
#' @param suppress.comment Logical, if TRUE, suppresses the output comments during the WASABI algorithm execution.
#' @param seed An optional integer to set the random seed for reproducibility. If NULL, no seed is set.
#'
#' @return A list with elements:
#' \describe{
#'   \item{wass_vec}{A numeric vector of length \code{L_max}, containing the Wasserstein distances achieved by the WASABI approximation for each value of \code{L}.}
#'   \item{output_list}{A list of length \code{L_max}, where each element contains the output of the WASABI algorithm for each value of \code{L}. Each element contains:
#'   \itemize{
#'      \item \code{particles} A matrix with \code{L} rows, each containing one of the WASABI particles, ordered by decreasing weight.
#'      \item \code{EVI} A vector of length \code{L}, containing expected VI associated to each particle.
#'      \item \code{wass.dist} A numeric value giving the Wasserstein distance achieved by the WASABI approximation.
#'      \item \code{part.psm} A list of length \code{L}, each element containing the posterior similarity matrix corresponding to the region of attraction of each particle. Only returned if \code{return_psm} is TRUE.
#'      \item \code{part.weights} A vector of length \code{L}, containing weight associated to each particle.
#'      \item \code{draws.assign} A vector containing the assignment of each MCMC sample to its closest particle.
#'   }}
#' }
#'
#'
#' @details Several initialization methods are available:
#' \itemize{
#'    \item \code{"average"} and \code{"complete"} initialize the particles using hierarchical clustering (choosing the $L$ ones with smallest EVI);
#'    \item \code{"fixed"} initializes the algorithm with a set of $L$ partition provided in \code{part.init};
#'    \item \code{"++"} uses an algorithm similar to k-means++, i.e. promotes diversity among initial centers by iteratively choosing the next center with probability proportional to its VI distance from the closest already chosen center;
#'    \item \code{"+++"} implements the k-means++ initialization using by choosing from the MCMC samples and the partitions obtained from hierarchical clustering (average and complete);
#'    \item \code{"topvi"} initializes with the $L$ partitions with smallest EVI, chosen from the ones generated by \code{"average"}, \code{"complete"} and \code{"fixed"} (if \code{part.init} is provided).
#' }
#' The WASABI algorithm iteratively updates the particles by computing the region of attractions (i.e. the assignment of each MCMC sample to the closest particle), in the N-update step, and finding the minEVI particle for each group, in the VI-search step.
#' The VI-search step can rely on different algorithms for finding the minEVI particle, depending on the \code{method} argument:
#' \itemize{
#'    \item \code{"average"} and \code{"complete"} use hierarchical clustering to find the minEVI particle;
#'    \item \code{"greedy"} uses the greedy algorithm implemented in \code{MinimiseEPL} of the \code{GreedyEPL} package;
#'    \item \code{"salso"} uses the \code{salso} package to find the minEVI particle.
#' }
#'
#' @seealso WASABI, WASABI_multistart
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' mu <- c(-1.1, 1.1)
#' prop <- c(0.5, 0.5)
#' n <- 300
#' components <- sample(1:2, size = n, replace = TRUE, prob = prop)
#' y <- rnorm(n, mean = mu[components], sd = 1)
#' est_model <- BNPmix::PYdensity(y = y,
#'                                mcmc = list(niter = 6000,
#'                                            nburn = 5000,
#'                                            model = "LS"),
#'                                output = list(out_type = "FULL", out_param = TRUE))
#' cls.draw = est_model$clust
#' psm=mcclust::comp.psm(cls.draw+1)
#' WASABI_elbow <- elbow(cls.draw, psm = psm)
#' plot(WASABI_elbow$wass_vec, type = "b")
#' }
elbow <- function(cls.draw, L_max = 10, psm = NULL,
                  multi.start = 1, ncores = 1,
                  method.init = "topvi", add_topvi = FALSE, lb = TRUE,
                  thin.init = NULL, part.init = NULL,
                  method = "average",
                  max.k = NULL, max.iter = 20, eps = 0.01,
                  mini.batch = NULL, extra.iter = NULL,
                  swap_countone = FALSE,
                  suppress.comment = TRUE, seed = NULL) {
  if (is.null(psm)) {
    psm <- mcclust::comp.psm(cls.draw)
  }
  # note: if thin.init is null, then default is thin = 10
  if (is.null(mini.batch) == TRUE) {
    mini.batch <- round(nrow(cls.draw) / 5)
  }

  output_list <- list()
  wass_vec <- numeric(L_max)

  for (ell in 1:L_max) {
    output_wvi <- WASABI_multistart(cls.draw, L = ell, psm = psm,
      multi.start = multi.start, ncores = ncores,
      method.init = method.init, add_topvi = add_topvi,
      lb = lb, thin.init = thin.init, part.init = part.init,
      method = method, max.k = max.k,
      max.iter = max.iter, eps = eps,
      mini.batch = mini.batch, extra.iter = extra.iter,
      swap_countone = swap_countone,
      suppress.comment = suppress.comment,
      seed = seed
    )
    output_list[[ell]] <- output_wvi
    wass_vec[ell] <- output_wvi$wass.dist
    cat(paste("Completed ", ell, "/", L_max, "\n"))
  }
  return(list(wass_vec = wass_vec, output_list = output_list))
}
