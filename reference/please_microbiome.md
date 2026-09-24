# Longitudinal microbiome abundances from the PLEASE study

A long-format data frame of genus-level relative abundances and
associated clinical covariates from the pediatric study of Lewis et al.
(2015). The analytic dataset was reconstructed from the public
`chvlyl/PLEASE` repository and processed following the filtering and
recoding conventions of the ZIBR package (Chen and Li, 2016).
Genus-level relative abundances were originally quantified from shotgun
metagenomic sequencing using MetaPhlAn 1.7.6 (Segata et al., 2012).

## Usage

``` r
please_microbiome
```

## Format

A data frame with 3186 rows and 7 variables. Each row is one
post-baseline observation of a genus in a subject, so the 18 genera
contribute 177 rows each (3 visits \\\times\\ 59 subjects).

- Genus:

  character. Genus name, e.g. `"g__Bacteroides"`.

- Sample:

  character. Sample identifier.

- Subject:

  character. Subject identifier (e.g. `"S5001"`).

- Time:

  numeric. Follow-up week: `1`, `4` or `8`. The baseline (week 0) is
  stored in `Baseline` and does not appear here.

- Treat:

  factor. Treatment arm with levels `"antiTNF"` (reference) and `"EEN"`
  (exclusive enteral nutrition).

- Baseline:

  numeric. Subject's relative abundance of the same genus at week 0, on
  its original proportion scale in \\\[0,1\]\\.

- Y:

  numeric. Relative abundance of the genus at the given post-baseline
  visit, in \\\[0,1\]\\. May contain zeros.

## Source

PLEASE repository: <https://github.com/chvlyl/PLEASE>

Original study: Lewis, J. D., Chen, E. Z., Baldassano, R. N., et al.
(2015). Inflammation, antibiotics, and diet as environmental stressors
of the gut microbiome in pediatric Crohn's disease. *Cell Host &
Microbe*, **18**(4), 489–500.
[doi:10.1016/j.chom.2015.09.008](https://doi.org/10.1016/j.chom.2015.09.008)

## Details

The dataset covers \\59\\ subjects (\\47\\ anti-TNF, \\12\\ EEN) with
observations at all four scheduled visits (baseline and weeks 1, 4 and
8) and \\18\\ bacterial genera, for a total of \\236\\ samples in the
wide format used by Chen and Li (2016). The long format provided here
contains the \\177\\ post-baseline observations per genus, together with
the baseline abundance replicated within subject and genus.

Processing followed the ZIBR workflow: samples with fewer than
\\10{,}000\\ non-human reads were removed; genera were retained if
present in more than \\40\\\\ of the remaining samples and if their
\\90\\th percentile of relative abundance, computed including zeros,
exceeded \\1\\\\; and the abundances of the retained genera were
renormalized to sum to one within each sample, so that each modeled
abundance is relative to the retained genera rather than to the whole
community. These steps preceded the restriction to the anti-TNF and EEN
arms (the partial enteral nutrition arm was excluded) and the selection
of subjects with observations at all four scheduled visits.

In the analysis reported in the companion paper, baseline abundance was
used as a subject-level covariate on its original proportion scale,
without centering or standardization, and both components of each
two-part model included baseline abundance, week, and treatment. One
genus, *Bacteroides*, had only three zeros among the post-baseline
observations and was excluded from the Beta–Vasicek model comparison,
leaving 17 genera; its observations are nonetheless retained here for
completeness.

## References

Chen, E. Z. and Li, H. (2016). A two-part mixed-effects model for
analyzing longitudinal microbiome compositional data. *Bioinformatics*,
**32**(17), 2611–2617.
[doi:10.1093/bioinformatics/btw308](https://doi.org/10.1093/bioinformatics/btw308)

Segata, N., Waldron, L., Ballarini, A., et al. (2012). Metagenomic
microbial community profiling using unique clade-specific marker genes.
*Nature Methods*, **9**(8), 811–814.
[doi:10.1038/nmeth.2066](https://doi.org/10.1038/nmeth.2066)

## Examples

``` r
data(please_microbiome)
head(please_microbiome)
#>                  Genus  Sample Subject Time   Treat    Baseline            Y
#> 5001-02 g__Bacteroides 5001-02   S5001    1 antiTNF 0.004514608 0.0023451675
#> 5001-03 g__Bacteroides 5001-03   S5001    4 antiTNF 0.004514608 0.0002086198
#> 5001-04 g__Bacteroides 5001-04   S5001    8 antiTNF 0.004514608 0.0000769457
#> 5002-02 g__Bacteroides 5002-02   S5002    1 antiTNF 0.732645555 0.5981502341
#> 5002-03 g__Bacteroides 5002-03   S5002    4 antiTNF 0.732645555 0.0698715533
#> 5002-04 g__Bacteroides 5002-04   S5002    8 antiTNF 0.732645555 0.0460852342

# Dimensions and structure
dim(please_microbiome)                     # 3186 x 7
#> [1] 3186    7
length(unique(please_microbiome$Genus))    # 18 genera
#> [1] 18
length(unique(please_microbiome$Subject))  # 59 subjects
#> [1] 59
table(please_microbiome$Treat)             # 47 antiTNF, 12 EEN
#> 
#> antiTNF     EEN 
#>    2538     648 

# Balanced within genus: 177 rows each
table(please_microbiome$Genus)
#> 
#>        g__Alistipes      g__Bacteroides  g__Bifidobacterium      g__Clostridium 
#>                 177                 177                 177                 177 
#>      g__Collinsella    g__Coprobacillus        g__Dialister            g__Dorea 
#>                 177                 177                 177                 177 
#>      g__Escherichia      g__Eubacterium g__Faecalibacterium      g__Haemophilus 
#>                 177                 177                 177                 177 
#>    g__Lactobacillus  g__Parabacteroides        g__Roseburia     g__Ruminococcus 
#>                 177                 177                 177                 177 
#>    g__Streptococcus      g__Veillonella 
#>                 177                 177 
```
