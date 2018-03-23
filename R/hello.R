# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' A print function
#'
#' Print anything you want after helo
#' @param x What you want to print. Defaults to World.
#' @keywords hello print
#' @export
#' @examples
#' hello("y'all")

hello <- function(x = "World!") {
  print(paste("Hello", x))
}
