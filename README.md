README
================
Emily Lacroix

Last updated: 1/13/2022

This repository (repo) contains the code and data to reproduce all
results and figures in the 2021 JGR: Biogeosciences publication:
**Mineral protection and resource limitations combine to explain
profile-scale carbon persistence**.

`PuuEke_AllData.xlsx` contains all of the geochemical data used in
analyses as well as some microbial data and a README tab explaining data
types and units. `MPC_16S_phyloseq_2019_07_23.rds` is a phyloseq object,
ready for import in R.

The R Markdown files are configured to read directly from the
spreadsheet and phyloseq object on the user’s desktop.

Each RMarkdown file regenerates all results and figures used in the
manuscript. The `.md` files in this repo show formatted code output. The
`.Rmd` files can be downloaded to run and edit on your local machine.

-   `Figures_Geochem.md` shows data processing for geochemical data.
-   `Figures_Microbial.md` shows data processing for microbial data.
-   `MPC_DADA2.R` shows the DADA2 pipeline used to process Illumina
    reads and the creation of the phyloseq object

Please contact the corresponding author(s) with questions or enhanced
access.
