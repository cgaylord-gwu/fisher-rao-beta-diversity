# Fisher–Rao Distance for Microbial Community Composition

Reference implementation and analysis code for Fisher–Rao distance as a
principled dissimilarity measure for microbial community compositional data.

**Authors:** Clark Gaylord, Marcos Pérez-Losada, Keith A. Crandall
**Affiliation:** Computational Biology Institute and Department of Biostatistics
& Bioinformatics, Milken Institute School of Public Health, The George
Washington University, Washington, DC, USA

---

## Overview

This repository accompanies the manuscript:

> Gaylord, C., Pérez-Losada, M., & Crandall, K.A. Fisher–Rao distance as a
> principled dissimilarity for microbial community composition. (in preparation)

After normalization to relative abundances, microbial community samples are
compositions on the probability simplex. Fisher–Rao distance is the geodesic
distance induced by the Fisher information metric on that simplex and is
uniquely characterized, up to scale, by Čencov's invariance theorem. This
repository provides a reference R implementation of the distance and the
analysis code used to compare it against Bray–Curtis dissimilarity on a
vermicompost microbiome dataset.

---

## Repository contents

```
R/
  fisher_rao.R                  # Reference implementation: fisher_rao_dist()
analysis/
  pairwise_distances_multilevel.R    # CV and pairwise distance comparisons
  permanova_fr_bc_multilevel.R       # PERMANOVA across taxonomic levels
  plot_permanova_multilevel.R        # PERMANOVA heatmap and R² figures
  three_stream_experimental_barplot.R  # Phylum-level composition figure (main text)
  three_stream_barplot_multilevel.R    # Composition figures, all 4 levels (supplement)
data/
  amplicon_native/16S/tables/{seqtab,taxa}.rds   # Native 16S amplicon stream
  megahit_contig/tables/{seqtab,taxa}.rds        # MEGAHIT->Barrnap->cutadapt->BBMap->DADA2 stream
  extracted_read/tables/{seqtab,taxa}.rds        # Read-based (BBDuk) extraction stream
figures/                            # Created by the analysis scripts on first run
```

---

## Reference implementation

The core function is in `R/fisher_rao.R`. It takes a sample-by-taxon
matrix of relative abundances (rows sum to 1) and returns a `dist` object
compatible with `vegan::adonis2`, `vegan::betadisper`, and standard
ordination functions.

```r
source("R/fisher_rao.R")

# p: matrix of relative abundances, samples as rows, taxa as columns
d <- fisher_rao_dist(p)

# Use with vegan
library(vegan)
adonis2(d ~ treatment + time, data = metadata)
```

The distance between samples *i* and *j* is:

$$d_{\mathrm{FR}}(i, j) = 2 \mathrm{arccos} \left(\sum_{k=1}^{K} \sqrt{p_{ik}\ p_{jk}}\right)$$

The formula is defined on the closed simplex, including the boundary; absent
taxa contribute zero to the sum and require no pseudocount or imputation.

---

## Reproducibility note

The processed ASV count tables (`data/`) used by the analysis scripts are
included in this repository, so `analysis/*.R` can be run directly from a
clean checkout — no external data download is required. Each script writes
its outputs to `figures/`, created on first run:

```r
Rscript analysis/pairwise_distances_multilevel.R
Rscript analysis/permanova_fr_bc_multilevel.R
Rscript analysis/plot_permanova_multilevel.R    # after permanova script, reads its CSV
Rscript analysis/three_stream_experimental_barplot.R   # main-text Phylum figure
Rscript analysis/three_stream_barplot_multilevel.R     # supplement Figures S4-S7, all 4 levels
```

These tables were derived from the raw sequencing data in BioProject
PRJNA777435 via an upstream bioinformatics pipeline entailing:
-  read quality control
-  shotgun assembly (MEGAHIT)
-  16S rRNA gene prediction (Barrnap)
-  V4 region extraction (cutadapt)
-  read-level abundance recovery (BBMap), and
-  amplicon sequence variant inference (DADA2 with SILVA 138.1 taxonomy)

Full pipeline documentation, including parameters and workflow code, will
accompany a separate bioinformatics-focused publication.

---

## Data

The vermicompost microbiome dataset used in the analysis is available from
NCBI under BioProject
[PRJNA777435](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA777435)
(Pérez-Losada et al., 2022, *Frontiers in Microbiology*,
doi:[10.3389/fmicb.2022.854423](https://doi.org/10.3389/fmicb.2022.854423)).

---

## Dependencies

R packages: `ggplot2`, `dplyr`, `tidyr`, `vegan`, `phyloseq`, `DADA2`

Assembly and extraction: MEGAHIT 1.2.9, Barrnap 0.9, cutadapt, BBMap

Taxonomy: SILVA 138.1

---

## Citation

If you use this implementation, please cite:

> Gaylord, C., Pérez-Losada, M., & Crandall, K.A. Fisher–Rao distance as a
> principled dissimilarity for microbial community composition. (in preparation)

---

## License

To be determined pending publication.
