## fisher_rao.R
##
## Clark Gaylord
## Computational Biology Institute, The George Washington University
##
## Developed with assistance from Claude (Anthropic).
##
## Reference implementation of Fisher-Rao distance for microbial community
## compositional data.
##
## For the theoretical background and empirical evaluation, see:
##
##   Gaylord C, Pérez-Losada M, Crandall KA.
##   Fisher-Rao distance as a principled dissimilarity for microbial community
##   composition. (in preparation)
##
## The Fisher-Rao distance between two compositions p and q over K taxa is:
##
##   d_FR(p, q) = 2 * arccos( sum_k sqrt(p_k * q_k) )
##
## This is the geodesic distance under the Fisher information metric on the
## probability simplex, appropriate for normalized proportions once sequencing
## depth has been conditioned out. The formula is defined on the closed simplex
## including the boundary; absent taxa contribute zero and require no
## pseudocount or imputation.
##
## Reference: Atkinson & Mitchell (1981); Miyamoto et al. (2024),
##   Information Geometry 7(2):311-354. doi:10.1007/s41884-024-00143-2

## -----------------------------------------------------------------------------
## fisher_rao_dist
##
## Compute pairwise Fisher-Rao distances for a sample-by-taxon matrix.
##
## Arguments:
##   x         Numeric matrix, samples as rows, taxa as columns.
##             Values are raw counts or relative abundances; rows are
##             normalized to sum to 1 internally.
##   clip      Logical (default TRUE). Clip the Bhattacharyya coefficient
##             to [0, 1] before applying arccos, guarding against floating-
##             point values marginally outside this range.
##
## Returns:
##   A 'dist' object of class "dist", compatible with vegan::adonis2,
##   vegan::betadisper, stats::cmdscale, and related functions.
##
## Examples:
##   # From a raw count matrix
##   d <- fisher_rao_dist(count_matrix)
##   vegan::adonis2(d ~ treatment, data = metadata)
##
##   # From a phyloseq object
##   otu <- as(phyloseq::otu_table(ps), "matrix")
##   if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
##   d <- fisher_rao_dist(t(otu))   # transpose: samples as rows
## -----------------------------------------------------------------------------

fisher_rao_dist <- function(x, clip = TRUE) {

  ## --- Input checks ----------------------------------------------------------
  if (!is.matrix(x) && !is.data.frame(x))
    stop("'x' must be a matrix or data frame (samples x taxa).")
  x <- as.matrix(x)
  if (any(x < 0))
    stop("'x' contains negative values; counts or proportions required.")
  if (nrow(x) < 2)
    stop("'x' must have at least two rows (samples).")

  ## --- Normalize to proportions ----------------------------------------------
  ## A small epsilon is added before normalizing, matching the behavior of the
  ## analysis code accompanying this manuscript. This has negligible effect on
  ## non-zero counts at typical sequencing depths and avoids 0/0 arithmetic for
  ## all-zero rows without requiring a separate guard.
  ##
  ## If you prefer exact zeros to remain exact (the theoretically pure behavior
  ## on the closed simplex, where absent taxa contribute zero to the
  ## Bhattacharyya sum), set eps <- 0. The arccos formula is well-defined at
  ## the simplex boundary and requires no imputation.
  eps <- 1e-10
  x <- x + eps
  row_sums <- rowSums(x)
  p <- x / row_sums

  ## --- Pairwise distances ----------------------------------------------------
  n <- nrow(p)
  nms <- rownames(p)
  d <- numeric(n * (n - 1L) / 2L)
  idx <- 1L

  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      bc <- sum(sqrt(p[i, ] * p[j, ]))   # Bhattacharyya coefficient
      if (clip) bc <- min(bc, 1.0)        # guard against fp > 1
      d[idx] <- 2.0 * acos(bc)
      idx <- idx + 1L
    }
  }

  ## --- Return as dist object -------------------------------------------------
  attr(d, "Size")   <- n
  attr(d, "Labels") <- if (!is.null(nms)) nms else as.character(seq_len(n))
  attr(d, "Diag")   <- FALSE
  attr(d, "Upper")  <- FALSE
  attr(d, "method") <- "fisher-rao"
  class(d) <- "dist"
  d
}


## -----------------------------------------------------------------------------
## fisher_rao_pair
##
## Compute Fisher-Rao distance between two composition vectors.
## Useful for scalar comparisons without constructing a full distance matrix.
##
## Arguments:
##   p, q    Numeric vectors of equal length. Normalized internally.
##   clip    Logical (default TRUE). See fisher_rao_dist.
##
## Returns:
##   Scalar numeric distance in [0, pi].
## -----------------------------------------------------------------------------

fisher_rao_pair <- function(p, q, clip = TRUE) {
  if (length(p) != length(q))
    stop("'p' and 'q' must have the same length.")
  if (any(c(p, q) < 0))
    stop("Negative values not permitted.")
  sp <- sum(p + 1e-10); sq <- sum(q + 1e-10)
  p <- (p + 1e-10) / sp
  q <- (q + 1e-10) / sq
  bc <- sum(sqrt(p * q))
  if (clip) bc <- min(bc, 1.0)
  2.0 * acos(bc)
}
