#' @export
summary.geecor <- function(object, ...)
{
  value <- geepack:::summary.geeglm(object)
  value$corstr <- object$.corstruct
  value$corr <- value$geese$correlation
  colnames(value$corr) <- c("Estimate", "san.se", "Wald", "Pr(>|W|)")

  # value$call <- object$callorig
  value
}
