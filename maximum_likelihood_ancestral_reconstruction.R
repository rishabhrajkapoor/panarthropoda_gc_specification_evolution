# ---------------------------------------------------------------
# SYM, ARD, Γ-rate, 2-rate Hidden-Rates (HRM) models &
# Pagel correlation tests (GC_mech vs Germ_band)
# ---------------------------------------------------------------
# Ancestral-state reconstruction, model comparison, and SVG plotting
# plus inverse & bidirectional Pagel tests with individual LRT P-values
# ---------------------------------------------------------------
# DATA INPUT ----------------------------------------------------
#   Provide a TSV/CSV with:
#       * iq_tree_label – must match phy$tip.label
#       * GC_mech       – "IND" or "GP"
#       * Germ_band     – "long", "short", "intermediate", or <NA>
# ---------------------------------------------------------------

library(phytools)   # loads ape & provides fitMk, fitgammaMk, fitHRM, ancr, fitPagel
library(ape)

# ---------- SETTINGS ----------
outfile         <- "ML_ancestral_reconstruction.pdf"

# Pie-chart colour & border
pie.colors      <- c("grey", "black")   # IND, GP
pie.border.col  <- "black"              # outline colour
pie.border.lwd  <- 1                    # outline thickness

pie.cex.node    <- 0.3
pie.cex.tip     <- 0.2

# ---------- LOAD TREE ----------
phy <- read.tree("mcmctree.tree")

# ---------- LOAD TRAIT DATA ----------
dat <- read.table("phylogenetic_data_with_substitutions.tsv",
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  check.names = FALSE)

label_col <- "iq_tree_label"
if (!label_col %in% names(dat))
  stop("Column 'iq_tree_label' not found.")

# Clean row names ------------------------------------------------
rownames(dat) <- trimws(dat[[label_col]])

# GC_mech --------------------------------------------------------
trait <- factor(trimws(dat$GC_mech),
                levels = c("IND","GP"))
names(trait) <- rownames(dat)
trait <- trait[phy$tip.label]

if (anyNA(trait))
  stop("Some tips are missing GC_mech values.")
if (length(trait) != length(phy$tip.label))
  stop("Tree and trait vector length mismatch.")

# Germ-band (recode to binary) ----------------------------------
if (!"Germ_band" %in% names(dat))
  stop("Column 'Germ_band' not found in the data table.")

gb_raw <- trimws(dat$Germ_band)

recode_fun <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  x_low <- tolower(x)
  if (x_low == "long")                      return("long")
  if (x_low %in% c("short","intermediate")) return("short")
  return(NA_character_)
}

gb_recode <- vapply(gb_raw, recode_fun, character(1))

gb <- factor(gb_recode, levels = c("short","long"))
names(gb) <- rownames(dat)
gb <- gb[phy$tip.label]

if (anyNA(gb))
  warning("Some tips are missing Germ_band values; they will be dropped in Pagel test.")

# ---------- MODEL FITTING (GC_mech) -----------------------------
SYM.args <- list(tree = phy, x = trait, model = "SYM")
ARD.args <- list(tree = phy, x = trait, model = "ARD")

fit_SYM <- do.call(fitMk, SYM.args)
fit_ARD <- do.call(fitMk, ARD.args)

gamma_available <- exists("fitgammaMk", mode = "function")
if (gamma_available) {
  fit_SYM.G <- do.call(fitgammaMk, c(SYM.args, list(k = 4)))
  fit_ARD.G <- do.call(fitgammaMk, c(ARD.args, list(k = 4)))
}

hrm_available <- exists("fitHRM", mode = "function")
if (hrm_available) {
  fit_HRM2 <- fitHRM(tree = phy, x = trait, model = "SYM", rate.cat = 2)
}

# ---------- COLLECT MODELS -------------------------------------
models <- list(SYM = fit_SYM,
               ARD = fit_ARD)
if (gamma_available) {
  models$SYM_G <- fit_SYM.G
  models$ARD_G <- fit_ARD.G
}
if (hrm_available) {
  models$HRM2 <- fit_HRM2
}

# ---------- MODEL COMPARISON TABLE -----------------------------
cat("\n--- MODEL COMPARISON via anova() ---\n")
anova_table <- switch(
  paste(gamma_available, hrm_available),
  "TRUE TRUE"  = anova(fit_SYM, fit_ARD, fit_SYM.G, fit_ARD.G, fit_HRM2),
  "TRUE FALSE" = anova(fit_SYM, fit_ARD, fit_SYM.G, fit_ARD.G),
  "FALSE TRUE" = anova(fit_SYM, fit_ARD, fit_HRM2),
  anova(fit_SYM, fit_ARD)
)
print(anova_table)

# Quick summary --------------------------------------------------
logL <- sapply(models, logLik)
AICc <- sapply(models, AIC)
cat("\n--- SUMMARY (logLik / AIC) ---\n")
print(rbind(logLik = sprintf("%.3f", logL),
            AIC   = sprintf("%.3f", AICc)))

# ---------- SELECT BEST MODEL ----------------------------------
best_lab <- names(which.min(AICc))
best_fit <- models[[best_lab]]
cat(sprintf("\nSelected model for reconstruction (lowest AIC): %s\n", best_lab))

anc <- ancr(best_fit)

## ----- 1. helper --------------------------------------------------
get_probs <- function(label, tips) {
  node <- findMRCA(phy, tips, type = "node")
  vec  <- anc$ace[as.character(node), ]          # 1×k numeric named vector
  data.frame(
    clade      = label,
    node_id    = node,
    t(vec),                       # transpose so each state is a column
    row.names  = NULL,
    check.rows = FALSE,
    check.names = FALSE
  )
}

## ----- 2. build the table ----------------------------------------
tbl <- do.call(
  rbind,
  list(
    get_probs("Diptera",                 c("Drosophila_melanogaster", "Phlebotomus_papatasi")),
    get_probs("Antliophora",             c("Drosophila_melanogaster", "Ctenocephalides_felis")),
    get_probs("Mecoptera+Siphonaptera",  c("Boreidae", "Ctenocephalides_felis")),
    get_probs("Lepidoptera",             c("Hepialidae", "Erebidae")),
    get_probs("Coleoptera",              c("Carabidae", "Meloidae")),
    get_probs("Chrysomelidae+Curculionidae",
              c("Sitophilus_oryzae", "Acanthoscelides_obtectus")),
    get_probs("Hymenoptera (Athalia clade)",
              c("Ichneumonidae", "Athalia_rosae")),
    get_probs("Holometabola",            c("Ichneumonidae", "Drosophila_melanogaster")),
    get_probs("Insecta",                 c("Machilidae", "Drosophila_melanogaster")),
    get_probs("Collembola",              c("Tomoceridae", "Allacma_fusca")),
    get_probs("Poduromorpha+Symphypleona",
              c("Tetrodontophora_bielanensis", "Allacma_fusca")),
    get_probs("Hyalidae+Talitridae",     c("Parhyale_hawaiensis", "Talitridae")),
    get_probs("Decapoda",                c("Macrobrachium_nipponense", "Penaeus_japonicus")),
    get_probs("Pancrustacea",            c("Cytherideidae", "Drosophila_melanogaster")),
    get_probs("Myriapoda",               c("Scolopendra_cingulata", "Glomeridae")),
    get_probs("Chelicerata",             c("Nymphonidae", "Parasteatoda_tepidariorum")),
    get_probs("Arthropoda",              c("Drosophila_melanogaster", "Parasteatoda_tepidariorum")),
    get_probs("Panarthropoda",           c("Drosophila_melanogaster", "Parachela"))
  )
)

## ----- 3. inspect & save -----------------------------------------
print(tbl, digits = 3)                    # quick look in the console
write.csv(tbl, file = "node_marginal_probabilities.csv",
          row.names = FALSE, quote = FALSE)

cat("\n✓ Saved as node_marginal_probabilities.csv\n")
# ---------- PLOT & SAVE PDF ------------------------------------
pdf(outfile, width = 8, height = 10)   # <- pdf() instead of svg()

plot(anc,
     args.plotTree   = list(fsize = 0.5),
     args.nodelabels = list(
       piecol = pie.colors,
       cex    = pie.cex.node,
       border = pie.border.col,
       lwd    = pie.border.lwd
     ),
     args.tiplabels  = list(
       piecol = pie.colors,
       cex    = pie.cex.tip,
       border = pie.border.col,
       lwd    = pie.border.lwd
     ))

title(main = sprintf("Ancestral reconstruction (%s model)", best_lab))

dev.off()
cat(sprintf("\nPDF saved to: %s\n", normalizePath(outfile, mustWork = FALSE)))

