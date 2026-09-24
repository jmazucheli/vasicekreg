# =============================================================================
# please_microbiome.R  (data-raw)
# Reconstructs the analytic dataset of Chen & Li (2016) from the public
# repository chvlyl/PLEASE (data from Lewis et al., 2015) and stores it as
# data/please_microbiome.rda for the vasicekreg package.
#
# Model covariates (both components): Baseline + Time(week) + Treat.
# Weeks 1/4/8 are modelled (3 time points); baseline = abundance at week 0.
#
# NB: the filtering thresholds (depth, prevalence, abundance) and the recoding
# follow the ZIBR package conventions, NOT the wording of the paper. The
# stopifnot() traps halt the script if anything diverges from the structure
# reported by Chen & Li (2016).
#
# Run interactively from the package root with:
#   source("data-raw/please_microbiome.R")
# =============================================================================

library(dplyr)
library(tibble)

# -----------------------------------------------------------------------------
# 1. Raw data (same URLs as the ZIBR vignette)
# -----------------------------------------------------------------------------
PLEASE.file <- paste0(
    "https://raw.githubusercontent.com/chvlyl/PLEASE/master/",
    "1_Data/Raw_Data/MetaPhlAn/PLEASE/",
    "G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result.xls"
)
PLEASE.raw <- read.table(
    PLEASE.file,
    sep = "\t", header = TRUE, row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
)
taxa.raw <- t(PLEASE.raw)

human.read.file <- paste0(
    "https://raw.githubusercontent.com/chvlyl/PLEASE/master/",
    "1_Data/Raw_Data/MetaPhlAn/Human_Reads/please_combo_human_reads.xls"
)
human.read <- read.table(
    human.read.file,
    sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE
)

n_samples_raw <- nrow(taxa.raw)
n_taxa_raw    <- ncol(taxa.raw)

# -----------------------------------------------------------------------------
# 2. Quality filters and pre-processing (ZIBR conventions)
# -----------------------------------------------------------------------------
low.depth <- subset(human.read, NonHumanReads < 10000)
taxa.raw  <- taxa.raw[
    !rownames(taxa.raw) %in% rownames(low.depth),
    ,
    drop = FALSE
]

f1 <- apply(taxa.raw, 2, function(x) sum(x > 0) > 0.4 * length(x))
f2 <- apply(taxa.raw, 2, function(x) quantile(x, 0.9) > 1)
taxa.filter <- taxa.raw[, f1 & f2, drop = FALSE]
taxa.filter <- 100 * sweep(taxa.filter, 1, rowSums(taxa.filter), FUN = "/")
taxa.data   <- taxa.filter

cat("Raw samples:", n_samples_raw, "| raw taxa:", n_taxa_raw, "\n")
cat("After depth < 10000 removed:", nrow(taxa.raw), "samples\n")
cat("After genus filter:", ncol(taxa.data), "genera\n")

# -----------------------------------------------------------------------------
# 3. Covariates
# -----------------------------------------------------------------------------
sample.info.file <- paste0(
    "https://raw.githubusercontent.com/chvlyl/PLEASE/master/",
    "1_Data/Processed_Data/Sample_Information/",
    "2015_02_13_Processed_Sample_Information.csv"
)
sample.info <- read.csv(sample.info.file, row.names = 1)

reg.cov <-
    data.frame(Sample = rownames(taxa.data), stringsAsFactors = FALSE) %>%
    left_join(rownames_to_column(sample.info, var = "Sample"), by = "Sample") %>%
    filter(Treatment.Specific != "PEN") %>%
    group_by(Subject) %>% filter(n() == 4) %>% ungroup() %>%
    mutate(
        Treat   = ifelse(Treatment.Specific == "antiTNF", 1, 0),
        Subject = paste0("S", Subject),
        Time    = dplyr::recode(
            as.character(Time),
            "1" = 0, "2" = 1, "3" = 4, "4" = 8,
            .default = NA_real_
        )
    ) %>%
    select(Sample, Subject, Time, Response, Treat) %>%
    as.data.frame()

# ---- AUDIT TRAPS -------------------------------------------------------------
stopifnot(
    "PEN was not removed"           = !any(is.na(reg.cov$Treat)),
    "Time recoding produced NA"     = !any(is.na(reg.cov$Time)),
    "Not every subject has 4 times" = all(table(reg.cov$Subject) == 4),
    "Total rows != 236"             = nrow(reg.cov) == 236,
    "Number of subjects != 59"      = length(unique(reg.cov$Subject)) == 59
)
n_anti <- length(unique(reg.cov$Subject[reg.cov$Treat == 1]))
n_een  <- length(unique(reg.cov$Subject[reg.cov$Treat == 0]))
cat("Subjects: antiTNF =", n_anti, "| EEN =", n_een, "(expected 47 and 12)\n")
stopifnot("Group composition != 47/12" = (n_anti == 47 && n_een == 12))

# -----------------------------------------------------------------------------
# 4. Baseline (t0) as covariate + follow-up samples (t = 1, 4, 8)
# -----------------------------------------------------------------------------
reg.cov.t1 <- subset(reg.cov, Time == 0)
rownames(reg.cov.t1) <- reg.cov.t1$Subject
reg.cov.t234 <- subset(reg.cov, Time != 0)
reg.cov.t234 <- data.frame(
    baseline.sample = reg.cov.t1[reg.cov.t234$Subject, "Sample"],
    reg.cov.t234,
    stringsAsFactors = FALSE
)

stopifnot(
    "baseline.sample missing"       = !any(is.na(reg.cov.t234$baseline.sample)),
    "baseline from another subject" =
        all(reg.cov.t1[reg.cov.t234$Subject, "Subject"] == reg.cov.t234$Subject),
    "Follow-up != 177 rows"         = nrow(reg.cov.t234) == 177
)

spe.all <- colnames(taxa.data)

# -----------------------------------------------------------------------------
# 5. Long format (all genera) -> please_microbiome
# -----------------------------------------------------------------------------
please_microbiome <- do.call(rbind, lapply(spe.all, function(spe) {
    data.frame(
        Genus    = spe,
        Sample   = reg.cov.t234$Sample,
        Subject  = reg.cov.t234$Subject,
        Time     = reg.cov.t234$Time,
        Treat    = reg.cov.t234$Treat,
        Baseline = taxa.data[reg.cov.t234$baseline.sample, spe] / 100,
        Y        = taxa.data[reg.cov.t234$Sample,          spe] / 100,
        stringsAsFactors = FALSE
    )
}))

stopifnot(
    "Y outside [0,1)"        = all(please_microbiome$Y >= 0 &
                                       please_microbiome$Y < 1),
    "Baseline outside [0,1)" = all(please_microbiome$Baseline >= 0 &
                                       please_microbiome$Baseline < 1),
    "Rows != genera * 177"   = nrow(please_microbiome) == length(spe.all) * 177
)

# Treatment as factor; antiTNF is the reference level.
please_microbiome$Treat <- factor(
    please_microbiome$Treat,
    levels = c(1, 0),
    labels = c("antiTNF", "EEN")
)

stopifnot(
    "Treat factor levels wrong" =
        identical(levels(please_microbiome$Treat), c("antiTNF", "EEN"))
)

cat(
    "please_microbiome:", nrow(please_microbiome), "rows |",
    length(spe.all), "genera x 177 obs\n"
)

# -----------------------------------------------------------------------------
# 6. Store as data/please_microbiome.rda
# -----------------------------------------------------------------------------
usethis::use_data(please_microbiome, overwrite = TRUE)
