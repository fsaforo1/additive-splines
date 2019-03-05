#' @export


# Theming functions for the 'Joy of P-splines'
JOPS_point = function(colour = grey(0.5), size = 1.5) {
  geom_point(colour = grey(0.5), size = 1.5)
}

#' @export
# Overall theme
JOPS_theme = function() {
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
}

#' @export
# Color ramp
JOPS_colors = function(n) rainbow_hcl(n, start = 10, end = 350)
