#' Longitudinal microbiome abundances from the PLEASE study
#'
#' A long-format data frame of genus-level relative abundances and
#' associated clinical covariates from the pediatric study of
#' Lewis et al. (2015). The analytic dataset was reconstructed from the
#' public \code{chvlyl/PLEASE} repository and processed following the
#' filtering and recoding conventions of the \pkg{ZIBR} package
#' (Chen and Li, 2016). Genus-level relative abundances were originally
#' quantified from shotgun metagenomic sequencing using MetaPhlAn 1.7.6
#' (Segata et al., 2012).
#'
#' @format A data frame with 3186 rows and 7 variables. Each row is one
#'   post-baseline observation of a genus in a subject, so the 18 genera
#'   contribute 177 rows each (3 visits \eqn{\times} 59 subjects).
#' \describe{
#'   \item{Genus}{character. Genus name, e.g. \code{"g__Bacteroides"}.}
#'   \item{Sample}{character. Sample identifier.}
#'   \item{Subject}{character. Subject identifier (e.g. \code{"S5001"}).}
#'   \item{Time}{numeric. Follow-up week: \code{1}, \code{4} or \code{8}.
#'     The baseline (week 0) is stored in \code{Baseline} and does not
#'     appear here.}
#'   \item{Treat}{factor. Treatment arm with levels \code{"antiTNF"}
#'     (reference) and \code{"EEN"} (exclusive enteral nutrition).}
#'   \item{Baseline}{numeric. Subject's relative abundance of the same
#'     genus at week 0, on its original proportion scale in \eqn{[0,1]}.}
#'   \item{Y}{numeric. Relative abundance of the genus at the given
#'     post-baseline visit, in \eqn{[0,1]}. May contain zeros.}
#' }
#'
#' @details
#' The dataset covers \eqn{59} subjects (\eqn{47} anti-TNF, \eqn{12} EEN)
#' with observations at all four scheduled visits (baseline and weeks 1, 4
#' and 8) and \eqn{18} bacterial genera, for a total of \eqn{236} samples
#' in the wide format used by Chen and Li (2016). The long format
#' provided here contains the \eqn{177} post-baseline observations per
#' genus, together with the baseline abundance replicated within subject
#' and genus.
#'
#' Processing followed the ZIBR workflow: samples with fewer than
#' \eqn{10{,}000} non-human reads were removed; genera were retained if
#' present in more than \eqn{40\%} of the remaining samples and if their
#' \eqn{90}th percentile of relative abundance, computed including zeros,
#' exceeded \eqn{1\%}; and the abundances of the retained genera were
#' renormalized to sum to one within each sample, so that each modeled
#' abundance is relative to the retained genera rather than to the whole
#' community. These steps preceded the restriction to the anti-TNF and
#' EEN arms (the partial enteral nutrition arm was excluded) and the
#' selection of subjects with observations at all four scheduled visits.
#'
#' In the analysis reported in the companion paper, baseline abundance
#' was used as a subject-level covariate on its original proportion
#' scale, without centering or standardization, and both components of
#' each two-part model included baseline abundance, week, and treatment.
#' One genus, \emph{Bacteroides}, had only three zeros among the
#' post-baseline observations and was excluded from the Beta--Vasicek
#' model comparison, leaving 17 genera; its observations are nonetheless
#' retained here for completeness.
#'
#' @source
#' PLEASE repository:
#' \url{https://github.com/chvlyl/PLEASE}
#'
#' Original study:
#' Lewis, J. D., Chen, E. Z., Baldassano, R. N., et al. (2015).
#' Inflammation, antibiotics, and diet as environmental stressors of the
#' gut microbiome in pediatric Crohn's disease.
#' \emph{Cell Host & Microbe}, \bold{18}(4), 489--500.
#' \doi{10.1016/j.chom.2015.09.008}
#'
#' @references
#' Chen, E. Z. and Li, H. (2016).
#' A two-part mixed-effects model for analyzing longitudinal
#' microbiome compositional data.
#' \emph{Bioinformatics}, \bold{32}(17), 2611--2617.
#' \doi{10.1093/bioinformatics/btw308}
#'
#' Segata, N., Waldron, L., Ballarini, A., et al. (2012).
#' Metagenomic microbial community profiling using unique clade-specific
#' marker genes.
#' \emph{Nature Methods}, \bold{9}(8), 811--814.
#' \doi{10.1038/nmeth.2066}
#'
#' @examples
#' data(please_microbiome)
#' head(please_microbiome)
#'
#' # Dimensions and structure
#' dim(please_microbiome)                     # 3186 x 7
#' length(unique(please_microbiome$Genus))    # 18 genera
#' length(unique(please_microbiome$Subject))  # 59 subjects
#' table(please_microbiome$Treat)             # 47 antiTNF, 12 EEN
#'
#' # Balanced within genus: 177 rows each
#' table(please_microbiome$Genus)
"please_microbiome"