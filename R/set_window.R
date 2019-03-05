#' @export

set_window = function(width = 6, height = 4.5, kill = T) {
  # Set default or custom figure size
  # If kill = T, all graphic windows will first be closed
  # Works only for Windows!
  if (kill)
    graphics.off()
  #dev.new(width, height, noRStudioGD = T)
  windows(width, height)
  set_panels(1, 1)
}

