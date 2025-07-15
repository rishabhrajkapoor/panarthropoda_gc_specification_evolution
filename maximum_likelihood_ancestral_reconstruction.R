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

# =================================================================
#                     PAGEL'S CORRELATION TESTS
# =================================================================
cat("\n--- Pagel's discrete correlation tests ---\n")

# Load subset tree ------------------------------------------------
phy_sub <- read.tree("insect_only.tree")

# Align trait vectors to subset tree -----------------------------
GC_sub <- trait[phy_sub$tip.label]
GB_sub <- gb[phy_sub$tip.label]

keep <- !(is.na(GC_sub) | is.na(GB_sub))
phy_sub <- drop.tip(phy_sub, phy_sub$tip.label[!keep])
GC_sub  <- GC_sub[keep]
GB_sub  <- GB_sub[keep]

# Encode as 0/1 numeric for fitPagel -----------------------------
GC_num <- as.numeric(GC_sub == "GP")      # IND = 0, GP = 1
GB_num <- as.numeric(GB_sub == "long")    # short = 0, long = 1

names(GC_num) <- names(GC_sub)
names(GB_num) <- names(GB_sub)

fit.GC_num<-fitPagel(phy_sub,GC_num,GB_num,dep.var="y")
print(fit.GC_num)
print(fit.GC_num$dependent.AIC)
plot(fit.GC_num)

fit.GB_num<-fitPagel(phy_sub,GC_num,GB_num,dep.var="x")
print(fit.GB_num)
print(fit.GB_num$dependent.AIC)
plot(fit.GB_num)

fit.bidirectional<-fitPagel(phy_sub,GC_num,GB_num)
print(fit.bidirectional)
print(fit.bidirectional$independent.AIC)
print(fit.bidirectional$dependent.AIC)
plot(fit.bidirectional)

#----PhyloGLM logistic regression-----------------

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

GC_sub <- trait[phy_sub$tip.label]
GB_sub <- gb[phy_sub$tip.label]

keep <- !(is.na(GC_sub) | is.na(GB_sub))
phy_sub <- drop.tip(phy_sub, phy_sub$tip.label[!keep])
GC_sub  <- GC_sub[keep]
GB_sub  <- GB_sub[keep]

GC_num <- as.numeric(GC_sub == "GP")      # IND = 0, GP = 1
GB_num <- as.numeric(GB_sub == "long")    # short = 0, long = 1

names(GC_num) <- names(GC_sub)
names(GB_num) <- names(GB_sub)

if (anyNA(gb))
  warning("Some tips are missing Germ_band values; they will be dropped in Pagel test.")
df <- data.frame(GC01 = GC_num,
                 GB01 = GB_num,
                 row.names = phy_sub$tip.label)
fit_phyloglm <- phyloglm(GC01 ~ GB01,
                         data = df,
                         phy  = phy_sub,
                         method = "logistic_MPLE")
cat("\n--- phyloglm (after subs filter) ---\n")
print(summary(fit_phyloglm))
cat("Adj. odds ratio (GB01):",
    round(exp(coef(fit_phyloglm)["GB01"]), 3), "\n\n")


dat <- read.table("phylogenetic_data_with_substitutions.tsv",
                  header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("iq_tree_label" %in% names(dat))
rownames(dat) <- trimws(dat$iq_tree_label)

# ---------- ALIGN DATA TO TREE TIPS ----------
tip_order  <- phy_sub$tip.label
dat_align  <- dat[match(tip_order, rownames(dat)), ]  # may introduce NAs

# ---------- BUILD VARIABLES ----------
# GC_mech: response (binary)
trait_raw <- trimws(dat_align$GC_mech)
GC        <- factor(trait_raw, levels = c("IND","GP"))
GC_num    <- as.numeric(GC == "GP")          # IND=0, GP=1

# Germ_band: 3-level predictor
gb_raw <- trimws(dat_align$Germ_band)
GB     <- factor(gb_raw, levels = c("short","intermediate","long"))

# ---------- FILTER OUT MISSING ----------
keep <- complete.cases(GC, GB)
phy_sub <- drop.tip(phy_sub, tip_order[!keep])
GC_num  <- GC_num[keep]
GB       <- GB[keep]

# ---------- BUILD DATA FRAME ----------
df <- data.frame(
  GC01 = GC_num,
  GB   = GB,
  row.names = phy_sub$tip.label
)

# ---------- PHYLOGENETIC LOGISTIC REGRESSION ----------
fit_gp <- phyloglm(
  GC01 ~ GB,
  data   = df,
  phy    = phy_sub,
  method = "logistic_MPLE"
)



