#' @export
save_PDF = function(fname = "scratch", folder = "../../Graphs") {
  # Save the PDF output in the proper place
  fpdf = paste(fname, ".pdf", sep = "")
  ffull = file.path(folder, fpdf)
  cat(ffull, "\n")
  dev.copy2pdf(file = ffull)
}

