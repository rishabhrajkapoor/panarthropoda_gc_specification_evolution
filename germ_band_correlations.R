# =================================================================
#                     PAGEL'S CORRELATION TESTS
# =================================================================

library(phytools)   # loads ape & provides fitMk, fitgammaMk, fitHRM, ancr, fitPagel
library(ape)

phy <- read.tree("insect_only.tree")

# ---------- LOAD TRAIT DATA ----------
dat <- read.table("insect_germ_band.tsv",
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

# =================================================================
#                     Phylogenetic logistic regression
# =================================================================
phy <- read.tree("insect_only.tree")
dat <- read.table("insect_germ_band.tsv",
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
library(phylolm)
# ---------- PHYLOGENETIC LOGISTIC REGRESSION ----------
fit_gp <- phyloglm(
  GC01 ~ GB,
  data   = df,
  phy    = phy_sub,
  method = "logistic_MPLE"
)

# ---------- OUTPUT RESULTS ----------
print(summary(fit_gp))
