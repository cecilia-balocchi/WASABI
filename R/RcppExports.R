# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute Variation of Information (VI) distance between two partitions
#'
#' @param c1 An integer vector of cluster labels for $n$ items. Should be 0-indexed.
#' @param c2 An integer vector of cluster labels for $n$ items. Should be 0-indexed.
#' @param K1 An integer, specifying the number of unique clusters in \code{cl1}.
#' @param K2 An integer, specifying the number of unique clusters in \code{cl2}.
#' @return The VI distance between \code{cl1} and \code{cl2}.
VI_compute_Rcpp <- function(c1, c2, K1, K2) {
    .Call(`_WASABI_VI_compute_Rcpp`, c1, c2, K1, K2)
}

#' Compute Variation of Information (VI) distance between two groups of partitions
#'
#' @param cls1 A matrix of $n$ columns, where each row contains the 0-indexed cluster labels for n items. 
#' @param cls2 A matrix of $n$ columns, where each row contains the 0-indexed cluster labels for n items. 
#' @param K1s An integer vector, specifying the number of clusters in each row of \code{cls1}. Must have length equal to the number of rows in \code{cls1}.
#' @param K2s An integer vector, specifying the number of clusters in each row of \code{cls2}. Must have length equal to the number of rows in \code{cls2}.
#' @return A matrix containing the VI distance between each row of \code{cls1} and \code{cls2}.
VI_Rcpp <- function(cls1, cls2, K1s, K2s) {
    .Call(`_WASABI_VI_Rcpp`, cls1, cls2, K1s, K2s)
}

