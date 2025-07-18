#' Run the WASABI algorithm multiple times with different starting points
#'
#' @description The WASABI algorithm can be used to summarize posterior uncertainty in Bayesian clustering by
#' identifying multiple representative partitions (''particles'') rather than a single
#' point estimate. It approximates the full posterior over clusterings with a small
#' set of $L$ weighted partitions, selected to minimize the Wasserstein distance (with
#' Variation of Information loss) between the posterior and its summary.
#' \code{WASABI_multistart} runs the algorithm from multiple starting points and
#' returns the run achieving the best WASABI approximation, i.e. the one minimizing the Wasserstein distance.
#'
#' @usage WASABI_multistart(cls.draw = NULL, psm = NULL, multi.start = 10, ncores = 1,
#'                   method.init = c("average", "complete", "fixed", "++", 
#'                                   "random_partition", "+++", "topvi"),
#'                   add_topvi = TRUE,
#'                   lb = FALSE, thin.init = NULL, part.init = NULL,
#'                   method = c("average", "complete", "greedy", "salso"),
#'                   max.k = NULL, L = 10, max.iter = 100, eps = 0.0001, mini.batch = 0,
#'                   extra.iter = NULL,
#'                   swap_countone = FALSE,
#'                   suppress.comment = TRUE,
#'                   return_psm = TRUE, seed = NULL)
#'
#'
#' @param cls.draw A matrix of the MCMC samples of partitions of $n$ data points.
#' @param psm The posterior similarity matrix obtained from MCMC samples of partitions stored in \code{cls.draw}.
#' @param multi.start Integer, the number of times to run the WASABI algorithm with different starting points.
#' @param ncores Integer, the number of cores to use for parallel processing, relies on parallel::mclapply, so can only be used on MacOs and Linux systems. If 1, no parallel processing is used.
#' @param method.init Initialization method. Options are "average", "complete", "fixed", "++", "random_partition", "+++", and "topvi".
#' @param add_topvi Logical, if TRUE, the "topvi" method is added once to the initialization methods. This method uses the top VI partitions from hierarchical clustering.
#' @param lb Logical, if TRUE, the lower bound for the VI is used in methods "average" or "complete".
#' @param thin.init Integer, thinning factor for the MCMC samples used to initialize the particles. If NULL, defaults to 10.
#' @param part.init A matrix of size \code{L} x \code{n}, containing the initial particles. Needs to be provided when \code{method.init = "fixed"}.
#' @param method The method used to find the partition with minimum EVI (minVI partition). Options are "average", "complete", "greedy", and "salso".
#' @param max.k Integer, the maximum number of clusters considered in the WASABI approximation (for "average", "complete", "greedy"). If NULL, it is set to the minimum of \code{max(Ks.draw)+10} and \code{ceiling(sqrt(n))}.
#' @param L Integer, the number of particles to be used in the WASABI approximation.
#' @param max.iter Integer, the maximum number of iterations for the WASABI algorithm.
#' @param eps Numeric, the convergence threshold for the WASABI algorithm. The algorithm stops when the difference in Wasserstein distance between two consecutive iterations is less than \code{eps}.
#' @param mini.batch Integer, the size of the mini-batch used in the WASABI algorithm. If 0, the full batch is used.
#' @param extra.iter Integer, the number of additional iterations to run after the mini-batch optimization. If NULL, defaults to 1 if \code{mini.batch > 0}.
#' @param swap_countone Logical, if TRUE, the WASABI algorithm allows swapping of particles with only one sample assigned to them (outlier-check step).
#' @param suppress.comment Logical, if TRUE, suppresses the output comments during the WASABI algorithm execution.
#' @param return_psm Logical, if TRUE, returns the posterior similarity matrix for each particle.
#' @param seed An optional integer, or a vector of length \code{multi.start} of integers to set the random seed for reproducibility. If NULL, no seed is set.
#' @return particles A matrix with \code{L} rows, each containing one of the WASABI particles, ordered by decreasing weight.
#' @return EVI A vector of size \code{L}, containing expected VI associated to each particle.
#' @return wass.dist A numeric value giving the Wasserstein distance achieved by the WASABI approximation.
#' @return part.psm A list of size \code{L}, each element containing the posterior similarity matrix corresponding to the region of attraction of each particle. Only returned if \code{return_psm} is set to TRUE.
#' @return part.weights A vector of size \code{L}, containing weight associated to each particle.
#' @return draws.assign A vector containing the assignment of each MCMC sample to its closest particle. Its length is equal to the number of rows in \code{cls.draw}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{particles}{A matrix with \code{L} rows, each containing one of the WASABI particles, ordered by decreasing weight.}
#'   \item{EVI}{A vector of length \code{L}, containing expected VI associated to each particle.}
#'   \item{wass.dist}{A numeric value giving the Wasserstein distance achieved by the WASABI approximation.}
#'   \item{part.psm}{A list of length \code{L}, each element containing the posterior similarity matrix corresponding to the region of attraction of each particle. Only returned if \code{return_psm} is TRUE.}
#'   \item{part.weights}{A vector of length \code{L}, containing weight associated to each particle.}
#'   \item{draws.assign}{A vector containing the assignment of each MCMC sample to its closest particle.}
#' }
#'
#' @details Several initialization methods are available:
#' \itemize{
#'    \item \code{"average"} and \code{"complete"} initialize the particles using hierarchical clustering (choosing the $L$ ones with smallest EVI);
#'    \item \code{"fixed"} initializes the algorithm with a set of $L$ partition provided in \code{part.init};
#'    \item \code{"++"} uses an algorithm similar to k-means++, i.e. promotes diversity among initial centers by iteratively choosing the next center with probability proportional to its VI distance from the closest already chosen center;
#'    \item \code{"random_partition"} initializes by randomly assigning each data point to one of $L$ groups, and the center for each group is chosen based on these assignments;
#'    \item \code{"+++"} implements the k-means++ initialization using by choosing from the MCMC samples and the partitions obtained from hierarchical clustering (average and complete);
#'    \item \code{"topvi"} initializes with the $L$ partitions with smallest EVI, chosen from the ones generated by \code{"average"}, \code{"complete"} and \code{"fixed"} (if \code{part.init} is provided).
#' }
#' It is recommended to use random initializations with multiple starts, such as \code{"++"}, \code{"random_partition"}, or \code{"+++"}, to avoid local minima.
#' The function uses \code{parallel::mclapply} to run the algorithm in parallel, however this is only supported on MacOS and Linux systems. 
#' If you are using Windows, set \code{ncores = 1} to run the algorithm sequentially.
#' 
#' The WASABI algorithm iteratively updates the particles by computing the region of attractions (i.e. the assignment of each MCMC sample to the closest particle), in the N-update step, and finding the minEVI particle for each group, in the VI-search step.
#' The VI-search step can rely on different algorithms for finding the minEVI particle, depending on the \code{method} argument:
#' \itemize{
#'    \item \code{"average"} and \code{"complete"} use hierarchical clustering to find the minEVI particle;
#'    \item \code{"greedy"} uses the greedy algorithm implemented in \code{MinimiseEPL} of the \code{GreedyEPL} package;
#'    \item \code{"salso"} uses the \code{salso} package to find the minEVI particle.
#' }
#' The recommended method is \code{"salso"}, as it is the most accurate while remaining efficient for larger datasets.
#' In case of larger datasets, it is recommended to use method \code{"average"}, or the \code{mini.batch} argument to speed up the algorithm.
#'
#' @seealso WASABI, elbow
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
#' # if running WASABI once, a non-random initialization is recommended, such as "topvi" or "average"
#' out_WASABI <- WASABI_multistart(cls.draw, psm = psm, L = 2, multi.start = 10, 
#'                                 method.init = "++", method = "salso")
#' 
#' # for larger n, it can be faster to use mini.batch
#' out_WASABI <- WASABI_multistart(cls.draw, psm = psm, L = 2, multi.start = 10, 
#'                                 method.init = "topvi", method = "salso", 
#'                                 mini.batch = 200, max.iter = 10, extra.iter = 10)
#'                                 
#' # Parallel processing is currently supported only on MacOS and Linux systems. 
#' # To activate it, set ncores > 1. If you are using Windows, set ncores = 1.
#' out_WASABI <- WASABI_multistart(cls.draw, psm = psm, L = 2, 
#'                                 ncores = 3, multi.start = 20,
#'                                 method.init = "topvi", method = "salso")
#' }
WASABI_multistart <- function(cls.draw = NULL, psm = NULL, multi.start = 10, ncores = 1,
                              method.init = c("average", "complete", "fixed", "++", "random_partition", "+++", "topvi"),
                              add_topvi = TRUE,
                              lb = FALSE, thin.init = NULL, part.init = NULL,
                              method = c("average", "complete", "greedy", "salso"),
                              max.k = NULL, L = 10, max.iter = 100, eps = 0.0001, mini.batch = 0,
                              extra.iter = NULL,
                              swap_countone = FALSE,
                              suppress.comment = TRUE,
                              return_psm = TRUE, seed = NULL) {
  if (!is.null(seed)) {
    if (length(seed) == 1) {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
      seeds <- NULL
    } else if (length(seed) == multi.start) {
      seeds <- seed
    } else {
      seeds <- NULL
      warning("wassersteinVI: seed needs to be either a vector of length equal to multi.start or length 1.
              Using NULL seed.")
    }
  } else {
    seeds <- NULL
  }

  if (ncores > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("The 'parallel' package is required for multi-core processing. Please install it or run with ncores = 1.")
    }
    if (multi.start < ncores) {
      ncores <- multi.start
    }
  }

  if (L == 1) {
    out <- WASABI(cls.draw, psm,
      method.init, lb,
      thin.init, part.init,
      method, max.k,
      L, max.iter, eps,
      mini.batch, extra.iter,
      swap_countone,
      suppress.comment,
      seed = NULL
    )
    return(out)
  }

  if (add_topvi & method.init != "topvi") {
    multi.start <- 1 + multi.start
    method.init_vec <- c("topvi", rep(method.init, length = multi.start))
  } else {
    method.init_vec <- rep(method.init, length = multi.start)
  }
  out_par <- parallel::mclapply(1:multi.start,
    function(g) {
      WASABI(cls.draw, psm,
        method.init_vec[g], lb,
        thin.init, part.init,
        method, max.k,
        L, max.iter, eps,
        mini.batch, extra.iter,
        swap_countone,
        suppress.comment,
        return_psm,
        seed = seeds[g]
      )
    },
    mc.cores = ncores
  )
  i_opt <- which.min(lapply(out_par, function(x) x$wass.dist))
  return(out_par[[i_opt]])
}
