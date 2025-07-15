# ---------- PACKAGES ----------
needed <- c("phytools", "phylolm", "ape")
to_get  <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_get))
  install.packages(to_get, repos = "https://cloud.r-project.org")
library(phytools)
library(phylolm)
library(ape)

# ---------- LOAD TREE ----------
phy_sub <- read.tree("insect_only.tree")

# ---------- LOAD TRAIT DATA ----------
dat <- read.table("oskar_data.tsv",
                  header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)

label_col <- "iq_tree_label"
stopifnot(label_col %in% names(dat))
rownames(dat) <- trimws(dat[[label_col]])
## 1.  Make a data frame keyed by the tree’s tips ------------------------
tip_order <- phy_sub$tip.label                       # the order we need
dat_sub   <- dat[match(tip_order, rownames(dat)), ]  # NA if a tip is missing

# ---------- BUILD VECTORS ----------
trait <- factor(trimws(dat$GC_mech), levels = c("IND", "GP"))
osk   <- factor(trimws(dat$oskar_present), levels = c("True", "False"))
busco <- as.numeric(dat$BUSCO_complete_single_copy)
names(trait) <- names(osk) <- names(busco) <- rownames(dat)

## 3.  One-liner to flag any missing value ------------------------------
keep <- complete.cases(trait, osk, busco)           # TRUE/FALSE for every tipphy_sub$tip.label
phy_sub <- drop.tip(phy_sub, phy_sub$tip.label[!keep])

trait  <- trait [phy_sub$tip.label]
osk    <- osk   [phy_sub$tip.label]
busco  <- busco [phy_sub$tip.label]

# ----------------------------------------------------
# 1.  PHYLOGENETIC LOGISTIC WITH BUSCO COVARIATE
# ----------------------------------------------------
GC_num  <- as.numeric(trait[keep] == "GP")
osk_num <- as.numeric(osk  [keep] == "True")
busco_z <- as.numeric(scale(busco[keep])) 

df <- data.frame(GC01 = GC_num,
                 osk01 = osk_num,
                 busco = busco_z,
                 row.names = phy_sub$tip.label)

# ---------- PHYLOGENETIC LOGISTIC (with BUSCO) ----------
fit_phyloglm <- phyloglm(GC01 ~ osk01 + busco,
                         data = df,
                         phy  = phy_sub,
                         method = "logistic_MPLE")
print(summary(fit_phyloglm))
cat("Adjusted odds ratio for osk presence:",
    round(exp(coef(fit_phyloglm)["osk01"]), 3), "\n")

# REPEAT ANLYSIS EXCLUDING SUBSTITUTIONS

# ---------- LOAD TREE ----------
phy_sub <- read.tree("insect_only.tree")

# ---------- DATA ----------
dat <- read.table("oskar_data.tsv", header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
rownames(dat) <- trimws(dat$iq_tree_label)

## ---- NEW filter: keep rows whose `subs` is empty/NA ------------------
dat <- dat[ is.na(dat$subs) | trimws(dat$subs) == "", ]

# ---------- ALIGN DATA TO TREE TIP ORDER ----------
tip_order <- phy_sub$tip.label
dat_aligned <- dat[match(tip_order, rownames(dat)), ]   # inserts NA where tip missing

# ---------- BUILD VECTORS ----------
trait  <- factor(trimws(dat_aligned$GC_mech),      levels = c("IND","GP"))
osk    <- factor(trimws(dat_aligned$oskar_present),levels = c("True","False"))
busco  <- as.numeric(dat_aligned$BUSCO_complete_single_copy)

# ---------- COMPLETE-CASE FILTER (after subs filter) ----------
keep <- complete.cases(trait, osk, busco)
phy_sub <- drop.tip(phy_sub, tip_order[!keep])
GC_num  <- as.numeric(trait[keep] == "GP")
osk_num <- as.numeric(osk   [keep] == "True")
busco_z <- as.numeric(scale(busco[keep]))    # centred & scaled
phy_sub$tip.label
# ---------- DATA FRAME FOR phyloglm ----------
df <- data.frame(GC01 = GC_num,
                 osk01 = osk_num,
                 busco = busco_z,
                 row.names = phy_sub$tip.label)

# ---------- PHYLOGENETIC LOGISTIC (with BUSCO) ----------
fit_phyloglm <- phyloglm(GC01 ~ osk01 + busco,
                         data = df,
                         phy  = phy_sub,
                         method = "logistic_MPLE")
cat("\n--- phyloglm (after subs filter) ---\n")
print(summary(fit_phyloglm))
cat("Adj. odds ratio (osk):",
    round(exp(coef(fit_phyloglm)["osk01"]), 3), "\n\n")
