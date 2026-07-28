#' @export
print.geecor <- function(x, digits = NULL, quote = FALSE, prefix = "", ...)
{
  xg <- x$geese
  if (is.null(digits))
    digits <- options()$digits
  else options(digits = digits)
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(unclass(x$coefficients), digits = digits)
  cat("\nDegrees of Freedom:", length(x$y), "Total (i.e. Null); ",
      x$df.residual, "Residual\n")
  if (!xg$model$scale.fix) {
    cat("\nScale Link:                  ", xg$model$sca.link)
    cat("\nEstimated Scale Parameters:  ")
    print(as.numeric(unclass(xg$gamma)), digits = digits)
  }
  else cat("\nScale is fixed.\n")
  cat("\nCorrelation:  Structure =", x$corstr, " ")
  if (pmatch(xg$model$corstr, "independence", 0) == 0) {
    cat("  Link =", xg$model$cor.link, "\n")
    cat("Estimated Correlation Parameters:\n")
    print(unclass(xg$alpha), digits = digits)
  }
  cat("\nNumber of clusters:  ", length(xg$clusz), "  Maximum cluster size:",
      max(xg$clusz), "\n\n")
  invisible(x)
}

