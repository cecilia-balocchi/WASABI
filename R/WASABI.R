#' Run the WASABI algorithm once
#'
#' @description The WASABI algorithm can be used to summarize posterior uncertainty in Bayesian clustering by
#' identifying multiple representative partitions (''particles'') rather than a single
#' point estimate. It approximates the full posterior over clusterings with a small
#' set of $L$ weighted partitions, selected to minimize the Wasserstein distance (with
#' Variation of Information loss) between the posterior and its summary.
#' As the algorithm is only guaranteed convergence to a local optimum, it is
#' recommended to run it multiple times with different initializations.
#'
#' @usage WASABI(cls.draw = NULL, L = 4, psm = NULL,
#'        method.init = c("++", "average", "complete", "fixed", 
#'                        "+++", "topvi"),
#'        lb = FALSE, thin.init = NULL, part.init = NULL,
#'        method = c("salso", "average", "complete", "greedy"),
#'        max.k = NULL, max.iter = 30, eps = 0.0001, mini.batch = 0,
#'        extra.iter = NULL,swap_countone = FALSE,suppress.comment = TRUE,
#'        return_psm = FALSE,seed = NULL, ...)
#' 
#' @param cls.draw A matrix of the MCMC samples of partitions of $n$ data points.
#' @param L Integer, the number of particles to be used in the WASABI approximation.
#' @param psm The posterior similarity matrix obtained from MCMC samples of partitions stored in \code{cls.draw}.
#' @param method.init Initialization method. Options are "++", "average", "complete", "fixed", "+++", and "topvi". Default is "++".
#' @param lb Logical, if TRUE, the lower bound for the VI is used in methods "average" or "complete".
#' @param thin.init Integer, thinning factor for the MCMC samples used to initialize the particles. If NULL, defaults to 10.
#' @param part.init A matrix of size \code{L} x \code{n}, containing the initial particles. Needs to be provided when \code{method.init = "fixed"}.
#' @param method The method used to find the partition with minimum EVI (minVI partition). Options are "salso", "average", "complete", and "greedy". Default is "salso".
#' @param max.k Integer, the maximum number of clusters considered in the WASABI approximation (for "average", "complete", "greedy"). If NULL, it is set to the minimum of \code{max(Ks.draw)+10} and \code{ceiling(sqrt(n))}.
#' @param max.iter Integer, the maximum number of iterations for the WASABI algorithm.
#' @param eps Numeric, the relative convergence threshold for the WASABI algorithm. The algorithm stops when the difference in Wasserstein distance between two consecutive iterations is less than \code{eps * log2(n)}.
#' @param mini.batch Integer, the size of the mini-batch used in the WASABI algorithm. If 0, the full batch is used.
#' @param extra.iter Integer, the number of additional iterations to run after the mini-batch optimization. If NULL, defaults to 3 if \code{mini.batch > 0}. Has to be greater than 0.
#' @param swap_countone Logical, if TRUE, the WASABI algorithm allows swapping of particles with only one sample assigned to them (outlier-check step).
#' @param suppress.comment Logical, if TRUE, suppresses the output comments during the WASABI algorithm execution.
#' @param return_psm Logical, if TRUE, returns the posterior similarity matrix for each particle.
#' @param seed An optional integer to set the random seed for reproducibility. If NULL, no seed is set.
#' @param ... Additional arguments passed to the salso algorithm.
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
#'    \item \code{"++"} uses an algorithm similar to k-means++, i.e. promotes diversity among initial centers by iteratively choosing the next center with probability proportional to its VI distance from the closest already chosen center;
#'    \item \code{"average"} and \code{"complete"} initialize the particles using hierarchical clustering (choosing the $L$ ones with smallest EVI);
#'    \item \code{"fixed"} initializes the algorithm with a set of $L$ partition provided in \code{part.init};
#'    \item \code{"+++"} implements the k-means++ initialization using by choosing from the MCMC samples and the partitions obtained from hierarchical clustering (average and complete). Warning: Initialization method \code{"+++"} is deprecated and will be removed in future versions. Use \code{"++"} instead.
#'    \item \code{"topvi"} initializes with the $L$ partitions with smallest EVI, chosen from the ones generated by \code{"average"}, \code{"complete"} and \code{"fixed"} (if \code{part.init} is provided).
#' }
#' The WASABI algorithm iteratively updates the particles by computing the region of attractions (i.e. the assignment of each MCMC sample to the closest particle), in the N-update step, and finding the minEVI particle for each group, in the VI-search step.
#' The VI-search step can rely on different algorithms for finding the minEVI particle, depending on the \code{method} argument:
#' \itemize{
#'    \item \code{"salso"} uses the \code{salso} package to find the minEVI particle.
#'    \item \code{"average"} and \code{"complete"} use hierarchical clustering to find the minEVI particle;
#'    \item \code{"greedy"} uses the greedy algorithm implemented in \code{MinimiseEPL} of the \code{GreedyEPL} package;
#' }
#'
#' @seealso WASABI_multistart, elbow
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
#'                                mcmc = list(niter = 15000,
#'                                            nburn = 5000,
#'                                            model = "LS"),
#'                                output = list(out_type = "FULL", out_param = TRUE))
#' cls.draw = est_model$clust
#' psm=mcclust::comp.psm(cls.draw+1)
#' # if running WASABI once, a non-random initialization is recommended, such as "topvi" or "average"
#' out_WASABI <- WASABI(cls.draw, psm = psm, L = 2,method.init = "topvi", method = "salso")
#' 
#' # for larger n, it can be faster to use mini.batch
#' out_WASABI <- WASABI(cls.draw, psm = psm, L = 2,method.init = "topvi", method = "salso", 
#'                      mini.batch = 200, max.iter = 20, extra.iter = 10)
#' }
WASABI <- function(cls.draw = NULL, L = 4, psm = NULL,
                   method.init = c("++", "average", "complete", "fixed", "+++", "topvi"),
                   lb = FALSE, thin.init = NULL, part.init = NULL,
                   method = c("salso", "average", "complete", "greedy"),
                   max.k = NULL, max.iter = 10, eps = 0.0001, mini.batch = 200,
                   extra.iter = NULL,
                   swap_countone = FALSE,
                   suppress.comment = TRUE,
                   return_psm = FALSE,
                   seed = NULL, ...) {
  method.init <- match.arg(method.init, choices = method.init)
  method <- match.arg(method, choices = method)
  if (is.null(cls.draw)) {
    stop("cls.draw must be provided")
  }

  if (method == "greedy") {
    if (!requireNamespace("GreedyEPL", quietly = TRUE)) {
      stop("The 'GreedyEPL' package is required for the 'greedy' method. Please install it.")
    }
  }
  if (lb == TRUE) {
    if (!requireNamespace("mcclust.ext", quietly = TRUE)) {
      stop("The 'mcclust.ext' package is required for using lower bound VI. Please install it.")
    }
  }

  psm_needed <- (method.init == "average" | method.init == "complete" | method.init == "topvi" | method.init == "+++")
  psm_needed <- psm_needed | (method == "average" | method == "complete") 
  psm_needed <- psm_needed | (lb == TRUE)
  if (psm_needed) {
    if (is.null(psm)) {
      stop("psm must be provided if method.init = avg or comp and lb = TRUE")
    }
    if (any(psm != t(psm)) | any(psm > 1) | any(psm < 0) | sum(diag(psm)) != nrow(psm)) {
      stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")
    }
  }

  if (method.init == "fixed" | (method.init == "topvi" & !is.null(part.init))) {
    if (dim(part.init)[1] != L) {
      stop("L initial particles must be supplied")
    }
    if (dim(part.init)[2] != dim(cls.draw)[2]) {
      stop("number of data points must be the same in initial particles and draws")
    }
  }
  
  if (method.init == "+++") {
    warning("Initialization method '+++' is deprecated and will be removed in future versions. Use '++' instead.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  cls.draw <- t(apply(cls.draw, 1, relabel_partition))
  cls.draw_relab <- cls.draw - 1
  Ks.draw <- apply(cls.draw_relab, 1, function(x) max(x)) + 1

  n <- dim(cls.draw)[2]
  S <- dim(cls.draw)[1]
  part <- matrix(0, L, n)
  part.evi <- rep(0, L)

  eps <- eps * log2(n)
  
  if (is.null(max.k)) max.k <- (max(Ks.draw) + 10)
  if (max.k < L) max.k <- L
  if (is.null(extra.iter)) extra.iter <- 3 # used only if mini.batch > 0
  if ((extra.iter == 0) & (mini.batch > 0)) {
    extra.iter <- 1
    warning("Override: extra.iter is set to 1, as it needs to be larger than 0.")
  }

  if (L == 1) {
    # no need to do the full algorithm (we only need the VI-search step):
    if (method == "average" | method == "complete") {
      out <- minVI_hclust(cls.draw_relab, Ks.draw,
        psm, method, max.k, lb,
        L = 1
      )
      part[L, ] <- out$part
      part.evi[L] <- out$evi
    }
    if (method == "greedy") {
      output_minepl <- GreedyEPL::MinimiseEPL(1 + cls.draw_relab,
        par =
          list(Kup = max.k, loss_type = "VI")
      ) # initializing part is randomly sampled
      part[L, ] <- output_minepl$decision
      part.evi[L] <- output_minepl$EPL
    }
    if (method == "salso") {
      output_salso <- salso::salso(x = cls.draw_relab, ...)
      part[L, ] <- as.numeric(output_salso)
      part.evi[L] <- as.numeric(attr(output_salso, "info")[4])
    }
    psm_ret = list(psm)
    if(return_psm & is.null(psm_ret)) {
      psm_ret <- list(mcclust::comp.psm(cls.draw))
    }
    output <- list(
      particles = part, EVI = part.evi, wass.dist = sum(part.evi),
      part.psm = psm_ret, part.weights = 1, draws.assign = NULL
    )
    return(output)
  }

  if (is.null(thin.init)) {
    thin.init <- 10
  }
  if (mini.batch > S) {
    stop("size of mini-batch must be smaller than number of MCMC draws")
  }

  # ------------------------------------------ Initialization ------------------------------------------

  if (method.init == "average" | method.init == "topvi") {
    out <- minVI_hclust(
      cls.draw_relab[seq(1, S, thin.init), ],
      Ks.draw[seq(1, S, thin.init)],
      psm, "average",
      max.k, lb, L
    )
    part <- out$part
    part.evi <- out$evi

    if (method.init == "topvi") {
      part_all <- part
      part.evi_all <- part.evi
    }
    if (suppress.comment == FALSE & method.init != "topvi") {
      cat(paste(
        paste0("Initial particle ", c(1:L)),
        paste0(": number of clusters = ", apply(part, 1, max)),
        paste0(", EVI = ", round(part.evi, 3)), "\n"
      ))
    }
  }

  if (method.init == "complete" | method.init == "topvi") {
    out <- minVI_hclust(
      cls.draw_relab[seq(1, S, thin.init), ],
      Ks.draw[seq(1, S, thin.init)],
      psm, "complete",
      max.k, lb, L
    )
    part <- out$part
    part.evi <- out$evi

    if (method.init == "topvi") {
      part_all <- rbind(part_all, part)
      part.evi_all <- c(part.evi_all, part.evi)
    }
    if (suppress.comment == FALSE & method.init != "topvi") {
      cat(paste(
        paste0("Initial particle ", c(1:L)),
        paste0(": number of clusters = ", apply(part, 1, max)),
        paste0(", EVI = ", round(part.evi, 3)), "\n"
      ))
    }
  }

  if (method.init == "fixed" | (method.init == "topvi" & !is.null(part.init))) {
    if (lb) {
      VI.fxd <- mcclust.ext::VI.lb(part.init, psm)
    } else {
      part.init_relab <- t(apply(part.init, 1, relabel_partition)) - 1
      Ks.fxd <- apply(part.init_relab, 1, function(x) max(x)) + 1
      VI.fxd <- sapply(
        1:L,
        function(i) {
          EVI_Rcpp(
            cls = part.init_relab[i, ], cls.draw = cls.draw_relab[seq(1, S, thin.init), ],
            Ks = Ks.fxd[i], Ks.draw = Ks.draw[seq(1, S, thin.init)]
          )
        }
      )
    }
    part <- part.init
    part.evi <- VI.fxd
    if (method.init == "topvi") {
      part_all <- rbind(part_all, part)
      part.evi_all <- c(part.evi_all, part.evi)
    }

    if (suppress.comment == FALSE & method.init != "topvi") {
      cat(paste(
        paste0("Initial particle ", c(1:L)),
        paste0(": number of clusters = ", apply(part, 1, max)),
        paste0(", EVI = ", round(part.evi, 3)), "\n"
      ))
    }
  }

  if (method.init == "++") {
    part <- matrix(0, L, n)
    cls.draw.thin <- cls.draw_relab[seq(1, S, thin.init), ] # cls.draw.thin will be 0-indexed
    Ks.thin <- Ks.draw[seq(1, S, thin.init)]
    S.thin <- nrow(cls.draw.thin)
    part[1, ] <- cls.draw.thin[sample(1:S.thin, 1, TRUE), ]
    K.part <- max(part[1, ]) + 1
    for (k in 2:L) {
      tmp <- VI_Rcpp(cls.draw.thin, part[1:(k - 1), , drop = FALSE], Ks.thin, K.part)
      tmp[tmp < 0] <- 0 # sometimes numerical errors cause equal particles to have negative (small) VI, so we set them to 0
      vi.init <- apply(tmp, 1, min)
      ik <- which(stats::rmultinom(1, 1, vi.init / sum(vi.init)) == 1)
      # ik = sample_max_jit(vi.init)
      ## this is a potential idea to sample partitions that have higher vi.init values, but I (Ceci) am not sure if it fixes the problem
      # tau = 1; ik = which(stats::rmultinom(1, 1, exp(tau * vi.init)/sum(exp(tau * vi.init)) ) == 1)
      # jit = 0.0001
      # prob = (vi.init/sum(vi.init) + jit)*(vi.init>0)/sum((vi.init/sum(vi.init) + jit)*(vi.init>0))
      # ik = which(stats::rmultinom(1, 1, prob) == 1)
      part[k, ] <- cls.draw.thin[ik, ]
      K.part <- c(K.part, max(part[k, ]) + 1)
    }
    if (lb) {
      VI.pp <- mcclust.ext::VI.lb(part, psm)
    } else {
      Ks.pp <- apply(part, 1, function(x) max(x)) + 1
      VI.pp <- sapply(
        1:L,
        function(i) {
          EVI_Rcpp(
            cls = part[i, ], cls.draw = cls.draw_relab[seq(1, S, thin.init), ],
            Ks = Ks.pp[i], Ks.draw = Ks.draw[seq(1, S, thin.init)]
          )
        }
      )
    }
    part.evi <- VI.pp
    part <- part + 1
    if (suppress.comment == FALSE) {
      cat(paste(
        paste0("Initial particle ", c(1:L)),
        paste0(": number of clusters=", apply(part, 1, max)),
        paste0(", EVI = ", round(part.evi, 3)), "\n"
      ))
    }
  }

  if (method.init == "+++") {
    
    part <- matrix(0, L, n)
    cls.draw.thin <- cls.draw_relab[seq(1, S, thin.init), ]
    Ks.thin <- Ks.draw[seq(1, S, thin.init)]
    S.thin <- nrow(cls.draw.thin)

    # add the partitions from hierarchical clustering
    # if(is.null(max.k)) max.k <- ceiling(sqrt(n)) # ceiling(n/8)
    hclust.comp <- stats::hclust(stats::as.dist(1 - psm), method = "complete")
    cls.comp <- t(apply(matrix(1:max.k), 1, function(x) stats::cutree(hclust.comp, k = x)))
    cls.draw.thin <- rbind(cls.draw.thin, cls.comp - 1)
    Ks.thin <- c(Ks.thin, apply(cls.comp, 1, function(x) max(x)))
    S.thin <- S.thin + max.k

    hclust.comp <- stats::hclust(stats::as.dist(1 - psm), method = "average")
    cls.comp <- t(apply(matrix(1:max.k), 1, function(x) stats::cutree(hclust.comp, k = x)))
    cls.draw.thin <- rbind(cls.draw.thin, cls.comp - 1)
    Ks.thin <- c(Ks.thin, apply(cls.comp, 1, function(x) max(x)))
    S.thin <- S.thin + max.k

    if (!is.null(part.init)) {
      if (dim(part.init)[2] != dim(cls.draw)[2]) {
        warning("part.init cannot be used for +++ initialization, wrong dimension.")
      } else {
        part.init <- t(apply(part.init, 1, relabel_partition))
        cls.draw.thin <- rbind(cls.draw.thin, part.init - 1)
        Ks.thin <- c(Ks.thin, apply(part.init, 1, function(x) max(x)))
        S.thin <- S.thin + dim(part.init)[1]
      }
    }

    part[1, ] <- cls.draw.thin[sample(1:S.thin, 1, TRUE), ]
    K.part <- max(part[1, ]) + 1
    for (k in 2:L) {
      tmp <- VI_Rcpp(cls.draw.thin, part[1:(k - 1), , drop = FALSE], Ks.thin, K.part)
      tmp[tmp < 0] <- 0 # sometimes numerical errors cause equal particles to have negative (small) VI, so we set them to 0
      vi.init <- apply(tmp, 1, min)
      ik <- which(stats::rmultinom(1, 1, vi.init / sum(vi.init)) == 1)
      # ik = sample_max_jit(vi.init)
      ## this is a potential idea to sample partitions that have higher vi.init values, but I (Ceci) am not sure if it fixes the problem
      # tau = 1; ik = which(stats::rmultinom(1, 1, exp(tau * vi.init)/sum(exp(tau * vi.init)) ) == 1)
      # jit = 0.0001
      # prob = (vi.init/sum(vi.init) + jit)*(vi.init>0)/sum((vi.init/sum(vi.init) + jit)*(vi.init>0))
      # ik = which(stats::rmultinom(1, 1, prob) == 1)
      part[k, ] <- cls.draw.thin[ik, ]
      K.part <- c(K.part, max(part[k, ]) + 1)
    }
    if (lb) {
      VI.ppp <- mcclust.ext::VI.lb(part, psm)
    } else {
      Ks.ppp <- apply(part, 1, function(x) max(x)) + 1
      VI.ppp <- sapply(
        1:L,
        function(i) {
          EVI_Rcpp(
            cls = part[i, ], cls.draw = cls.draw_relab[seq(1, S, thin.init), ],
            Ks = Ks.ppp[i], Ks.draw = Ks.draw[seq(1, S, thin.init)]
          )
        }
      )
    }
    part.evi <- VI.ppp
    # let's bring them back to 1-index
    part <- part + 1
    if (suppress.comment == FALSE) {
      cat(paste(
        paste0("Initial particle ", c(1:L)),
        paste0(": number of clusters=", apply(part, 1, max)),
        paste0(", EVI = ", round(part.evi, 3)), "\n"
      ))
    }
  }

  if (method.init == "topvi") {
    # these also avoids repeated partitions
    tmp_reorder <- reorder_find_unique(dim(part_all)[1], dim(part_all)[2], part_all)
    part.evi_all <- part.evi_all[tmp_reorder$config_index]
    part_all <- part_all[tmp_reorder$config_index, ]
    sortvals <- sort(part.evi_all, index.return = T, decreasing = F)
    part.evi <- sortvals$x[1:L]
    part <- part_all[sortvals$ix[1:L], ]

    if (suppress.comment == FALSE) {
      cat(paste(
        paste0("Initial particle ", c(1:L)),
        paste0(": number of clusters=", apply(part, 1, max)),
        paste0(", EVI = ", round(part.evi, 3)), "\n"
      ))
    }
  }

  # ------------------------------------------ End of Initialization ------------------------------------------

  if (max.iter == 0) {
    output <- list(
      particles = part, EVI = part.evi, wass.dist = NA,
      part.psm = NA, part.weights = NA, draws.assign = NA
    )
    return(output)
  }

  part_relab <- t(apply(part, 1, relabel_partition)) - 1
  Ks.part <- apply(part_relab, 1, function(x) max(x)) + 1

  iter <- 1
  part_new <- part
  part_new_relab <- part_relab
  Ks.part_new <- Ks.part
  part.evi_new <- rep(0, L)
  assign.vi <- rep(0, S)
  counts <- rep(0, L)
  diff <- abs(sum(part.evi * counts / S) - sum(part.evi_new * counts / S))
  if (diff == 0) {
    diff <- 1
  }
  # ------------------------------------------ Main optimization ------------------------------------------
  while (iter <= max.iter & diff > eps) {
    # Compute VI for all samples and initial centers
    mb.samp <- c(1:S)
    if (mini.batch > 0) {
      mb.samp <- sample(1:S, mini.batch, replace = FALSE)
    }

    out <- particle_search(
      cls.draw_relab[mb.samp, ], Ks.draw[mb.samp],
      part_relab, Ks.part, part.evi, L, method, swap_countone, max.k, lb, suppress.comment, (iter == 1), ...)

    diff <- abs(sum(part.evi * counts / (mini.batch + (mini.batch == 0) * S)) - sum(out$evi * out$counts / (mini.batch + (mini.batch == 0) * S)))
    if (iter == 1) diff <- 1
    part_relab <- out$part_relab
    part <- part_relab + 1
    part.evi <- out$evi
    Ks.part <- out$Ks.part
    counts <- out$counts

    if (suppress.comment == FALSE) {
      cat(paste("Iteration =", iter, "\n"))
      cat(paste(
        paste0("Particle ", c(1:L)), paste0(": number of clusters=", apply(part, 1, max)),
        paste0(", EVI = ", round(part.evi, 3)), paste0(", sumVI = ", round(part.evi * counts / (mini.batch + (mini.batch == 0) * S), 3)),
        paste0(", w= ", round(counts / (mini.batch + (mini.batch == 0) * S), 3)), "\n"
      ))
      cat(paste("Wasserstein dist =", sum(part.evi * counts / (mini.batch + (mini.batch == 0) * S)), "\n"))
    }
    iter <- iter + 1
  }

  # If minibatch, run as many as another max.iter steps to "converge" and get labels for all
  if (mini.batch > 0) {
    if (suppress.comment == FALSE) cat(paste("*Running full batch after mini-batch*\n"))
    diff <- 1
    while (iter <= (max.iter + extra.iter) & diff > eps) {
      out <- particle_search(
        cls.draw_relab, Ks.draw,
        part_relab, Ks.part, part.evi, L, method, swap_countone, max.k, lb, suppress.comment, FALSE, ...)

      diff <- abs(sum(part.evi * counts / S) - sum(out$evi * out$counts / S))
      part_relab <- out$part_relab
      part <- part_relab + 1
      part.evi <- out$evi
      Ks.part <- out$Ks.part
      counts <- out$counts

      if (suppress.comment == FALSE) {
        cat(paste("Iteration =", iter, "\n"))
        cat(paste(
          paste0("Particle ", c(1:L)), paste0(": number of clusters=", Ks.part),
          paste0(", EVI = ", round(part.evi, 3)), paste0(", sumVI = ", round(part.evi * counts / S, 3)),
          paste0(", w= ", round(counts / S, 3)), "\n"
        ))
        cat(paste("Wasserstein dist =", sum(part.evi * counts / S), "\n"))
      }
      iter <- iter + 1
    } # end while
  }

  assign.vi <- out$assign.vi

  assign.vi_ordered <- match(1:L, order(counts, decreasing = T))[assign.vi]
  tmp_order <- sort(counts, index.return = TRUE, decreasing = TRUE)
  tmp_order <- tmp_order$ix

  part_ordered <- part[tmp_order, ]
  part.evi_ordered <- part.evi[tmp_order]
  counts_ordered <- counts[tmp_order]

  psm_K_ordered <- NULL
  if (return_psm) {
    for (l in 1:L) {
      psm_K_ordered[[l]] <- mcclust::comp.psm(matrix(cls.draw[assign.vi_ordered == l, ], counts_ordered[l], n))
    }
  }

  output <- list(
    particles = part_ordered, EVI = part.evi_ordered, wass.dist = sum(part.evi_ordered * counts_ordered / S),
    part.psm = psm_K_ordered, part.weights = counts_ordered / dim(cls.draw)[1], draws.assign = assign.vi_ordered
  )
  return(output)
}
