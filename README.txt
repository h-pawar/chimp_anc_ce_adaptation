# Pawar_et.al._code-HO/

- This directory contains all the code used to make the following figures in Pawar et al. 2022: 2, 4, S2, S3, S4.

## scripts/

- Contains all the scripts used to create the figures.

### Pawar_et.al._figure4-top.R and Pawar_et.al._figure4-bottom.R

- These R scripts make subplots in figure 4 showing selection statistics across the CD4 region. 

### Pawar_et.al._candidate_vcfs.Rmd

- Adds the EPO ancestral alleles to the VCF from de Manuel et al. (2016) and makes VCFs for each 3P-CLR empirical tail. These VCFs are then input into Pawar_et.al._figures_2_S2_S3_S4.ipynb.

### Pawar_et.al._figures_2_S2_S3_S4.ipynb

- Makes the plots showing allele frequency patterns of 3P-CLR candidate windows.

## input/

- Contains the input for the scripts.

### Clean_Callable_Pan_troglodytes_ALL_Concat_1M_HWE_gls_Hetexcess.vcf.gz

- VCF from de Manuel et al. (2016).

### pan_troglodytes_ancestor_CHIMP2.1.4_e86/

- Contains EPO alignment data (downloaded and formatted in Pawar_et.al._candidate_vcfs.Rmd).

### window_subset_0.5percent.txt, window_subset_0.1percent.txt and window_subset_0.05percent.txt

- Coordinates of 3P-CLR candidate windows at each empirical threshold (0.5%. 0.1%, 0.05%).

### *.rds

- All .rds files are input for Pawar_et.al._figure4-top.R and Pawar_et.al._figure4-bottom.R

## data/

- This contains reformatted data from the input/ directory output from Pawar_et.al._candidate_vcfs.Rmd which is the input for Pawar_et.al._figures_2_S2_S3_S4.ipynb.

## output/

- Contains the figures outputted by Pawar_et.al._figures_2_S2_S3_S4.ipyn, Pawar_et.al._figure4-top.R and Pawar_et.al._figure4-bottom.R. 




