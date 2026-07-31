nice_collapse <- function(x)
{
  glue::glue_collapse(
    encodeString(
      encodeString(x, quote = "`"),
      quote = "\""),
    sep = ", ", last = " or ")
}
