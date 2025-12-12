# brca2-repro
Study of reproducibility conducted for BU/BINF-3350.

## Reproducibility

This project uses `renv` for dependency management.

To reproduce the R environment:
```r
install.packages("renv")
renv::restore()
```

R version used: 4.5.2
See `session_info.txt` for details.

### TCGA Data Download

TCGA-BRCA RNA-seq data was obtained from the NCI Genomic Data Commons (GDC) 
using the GDC Data Transfer Tool and the provided manifest file. Cohort 
metadata, clinical annotations, and sample sheets used to construct the 
analysis-ready dataset are already included to ensure exact reproducibility of 
the selected cohort.

To reproduce:

1. Install the GDC Data Transfer Tool:
   https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
   
2. Run:
   ```bash
   bash scripts/00_download_GDC_data.sh
   ```
