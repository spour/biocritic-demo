## Sample sheet helpers.

read_sample_sheet <- function(path) {
  sheet <- read.delim(path, stringsAsFactors = FALSE)
  sheet$group <- factor(sheet$group, levels = c("Basal", "LP", "ML"))
  sheet$lane <- factor(sheet$lane)
  rownames(sheet) <- sheet$sample
  sheet
}
