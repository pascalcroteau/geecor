#' @export
summary.geecor <- function(object, ...)
{
  value <- geepack:::summary.geeglm(object)
  value$corstr <- object$corstr
  value$call <- object$callorig
  value
}
