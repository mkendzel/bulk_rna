# Helper to build bubble-size legend limits and breaks from the plotted values
bubble_size_scale <- function(x, n = 3) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(list(limits = NULL, breaks = waiver()))

  r <- range(x)
  d <- if (diff(r) >= 10) 0 else 1
  lim <- c(floor(r[1] * 10^d), ceiling(r[2] * 10^d)) / 10^d

  if (diff(lim) == 0) return(list(limits = lim + c(-1, 1) / 10^d, breaks = lim[1]))

  list(limits = lim,
       breaks = unique(round(seq(lim[1], lim[2], length.out = n), d)))
}
