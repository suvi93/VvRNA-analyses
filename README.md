# VvRNA-analyses

Scripts used for the RNA-seq analyses in:

> **Jayaprasad, S., Schielzeth, H. & Palacios-Gimenez, O.M. (2026).**
> Evolutionary Dynamics of Dosage Compensation and Sex-biased Gene Expression in Morabine Grasshopper *Vandiemenella viatica*.
> *Genome Biology and Evolution*, 18(2): evag026.
> [https://doi.org/10.1093/gbe/evag026](https://doi.org/10.1093/gbe/evag026)

[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fgbe%2Fevag026-blue)](https://doi.org/10.1093/gbe/evag026)
[![Journal](https://img.shields.io/badge/Journal-Genome%20Biology%20and%20Evolution-orange)](https://academic.oup.com/gbe)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-276DC3?logo=r)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

---

## Overview

This repository contains the scripts used to analyze dosage compensation and sex-biased gene expression in *Vandiemenella viatica* morabine grasshoppers, comparing the ancestral X chromosome (P24X0 race) with derived neo-sex chromosomes (P24XY race) that arose via X–autosome fusion. The pipeline covers RNA-seq read alignment, differential expression, and the generation of all main figures in the paper.

## Repository contents

| Script | Description |
|---|---|
| `alignment.sh` | Aligns all transcriptome samples to the reference genome |
| `deseq2.R` | Differential gene expression analysis (DESeq2) |
| `pca.R` | Principal component analysis of expression data |
| `heatmap.R` | Heatmap of the top 25 up- and down-regulated DEGs per tissue, classified by chromosome |
| `box_plot_fig2.R` | Generates Figure 2 |
| `box_plot_fig3.R` | Generates Figure 3 |
| `box_plot_fig4.R` | Generates Figure 4 |
| `sbg_xavstackedbar_fig5.R` | Generates Figure 5 (sex-biased gene distribution across X/autosome arms) |

## Workflow

```
alignment.sh  →  deseq2.R  →  pca.R / heatmap.R / box_plot_fig*.R / sbg_xavstackedbar_fig5.R
```

1. **`alignment.sh`** aligns raw RNA-seq reads for all samples.
2. **`deseq2.R`** takes the resulting count data and performs differential expression testing between sexes/tissues.
3. The remaining `R` scripts consume the DESeq2 output to produce the exploratory (PCA, heatmap) and publication figures (Figs. 2–5).

## Requirements

- R ≥ 4.0 with `DESeq2`, `ggplot2`, and standard tidyverse packages
- A read aligner as specified in `alignment.sh` (e.g., STAR/HISAT2 — see script for exact tool/version)

## Data availability

- **Genome assemblies** (P24X0 and P24XY races): NCBI BioProject [PRJNA945230](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA945230)
- **RNA-seq reads**: NCBI BioProject [PRJNA668746](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA668746)

## Related repositories

- Genome assembly: Li et al. (2024), *Chromosome-level genome assembly of the morabine grasshopper Vandiemenella viatica*
- Neo-sex chromosome recombination dynamics: Jayaprasad et al. (2024), *Molecular Ecology*, [10.1111/mec.17567](https://doi.org/10.1111/mec.17567)

## Citation

If you use these scripts, please cite:

```bibtex
@article{jayaprasad2026dosage,
  author  = {Jayaprasad, Suvratha and Schielzeth, Holger and Palacios-Gimenez, Octavio M.},
  title   = {Evolutionary Dynamics of Dosage Compensation and Sex-biased Gene Expression in Morabine Grasshopper Vandiemenella viatica},
  journal = {Genome Biology and Evolution},
  volume  = {18},
  number  = {2},
  pages   = {evag026},
  year    = {2026},
  doi     = {10.1093/gbe/evag026}
}
```

## Contact

Suvratha Jayaprasad — [GitHub](https://github.com/suvrathaprasad) · [portfolio](https://suvrathaprasad.github.io)

## License

This project is licensed under the [MIT License](LICENSE).
