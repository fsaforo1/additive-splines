#' @export

set_panels = function(rows = 1, cols = 1) {
  # Set number of rows and column for multiple panels
  par(mfrow = c(rows, cols))
  par(mar = c(4, 3, 2, 2))  # Margins around panels
  par(mgp = c(1.6, 0.6, 0))  # Distance of labels and tick marks
  par(tcl = -0.4)  # Tick length
}

