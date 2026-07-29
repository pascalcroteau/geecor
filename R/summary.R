#' @export
summary.geecor <- function(object, ...)
{
  value <- geepack:::summary.geeglm(object)
  value$corstr <- object$.corstruct
  # value$call <- object$callorig
  value
}
