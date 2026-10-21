## Mouse mammary epithelium: basal, luminal progenitor (LP) and mature luminal (ML) cells.
## Follows Law et al. (2016) "RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR",
## F1000Research 5:1408 (Bioconductor workflow RNAseq123).

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(jsonlite)
})
source("R/samples.R")
dir.create("results", showWarnings = FALSE)

files <- list.files("data", pattern = "^GSM.*\\.txt\\.gz$")
x <- readDGE(files, path = "data", columns = c(1, 3))
sheet <- read_sample_sheet("data/sample_sheet.tsv")
colnames(x) <- sheet$sample
x$samples$group <- sheet$group
x$samples$lane <- sheet$lane

group <- x$samples$group
lane <- x$samples$lane

keep.exprs <- filterByExpr(x, group = group)
x <- x[keep.exprs, , keep.lib.sizes = FALSE]
x <- calcNormFactors(x, method = "TMM")

design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))

contr.matrix <- makeContrasts(
  BasalvsLP = LP - Basal,
  BasalvsML = Basal - ML,
  LPvsML = LP - ML,
  levels = colnames(design))

v <- voom(x, design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
dt <- decideTests(efit)
print(summary(dt))

tfit <- treat(vfit, lfc = 1)
dt_treat <- decideTests(tfit)

basal.vs.lp <- topTreat(tfit, coef = 1, n = Inf)
write.csv(basal.vs.lp, "results/basal_vs_lp_treat.csv")

metrics <- list(
  n_genes = nrow(x),
  de_BasalvsLP = sum(dt[, "BasalvsLP"] != 0),
  de_BasalvsML = sum(dt[, "BasalvsML"] != 0),
  de_LPvsML = sum(dt[, "LPvsML"] != 0),
  treat_BasalvsLP_up = sum(dt_treat[, "BasalvsLP"] > 0),
  treat_BasalvsLP_down = sum(dt_treat[, "BasalvsLP"] < 0),
  top_treat_gene = rownames(basal.vs.lp)[1],
  top_treat_logFC = basal.vs.lp$logFC[1]
)
write_json(metrics, "metrics.json", auto_unbox = TRUE, digits = 8)
