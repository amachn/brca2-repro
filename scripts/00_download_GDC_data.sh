#!/bin/bash

# Download TCGA-BRCA RNA-seq data using the GDC Data Transfer Tool
# Requirements:
#   - gdc-client installed and on PATH
#   - gdc_manifest.txt present in project root

gdc-client download \
  --manifest ../data/raw/gdc_data/gdc_manifest.txt \
  --dir ../data/raw/gdc_data/expression
