# Examining the neurostructural architecture of intelligence: The Lothian Birth Cohort 1936 Study

## Overview
This repository contains the R and MATLAB code used in "Examining the neurostructural architecture of intelligence: The Lothian Birth Cohort 1936 Study (LBC1936)". The analyses were conducted using 697 participants (~72.5 years of age) from the Lothian Birth Cohort 1936 (LBC1936). Two approaches were employed to model major cognitive domains (processing speed, crystallised ability, memory, and visuospatial ability) as well as general cognitive functioning ('g') using confirmatory factor analysis (CFA) and 13 cognitive tests: 1) a bifactor model that separated g-related variance from cognitive domain scores, providing latent factor scores for each cognitive domain and g; and 2) g as a single-order factor indicated by all 13 test scores. Global analysis investigated the relationship between g and global brain volumes (total brain, grey matter, normal-appearing white matter, and white matter hypointensities), and for white matter microstructure - general factors of fractional anisotropy (gFA) and mean diffusivity (gMD). The vertex-wise analysis investigated the relationship between cognitive domains and cortical volume, measured at 327,684 vertices across the cortical surface.

## Scripts
 - [NeuralArch_RCode.R](NeuralArch_RCode.R) - CFA models and global brain analyses
 - [NeuroArch_surfstat_script_Apr2024.m](NeuroArch_surfstat_script_Apr2024.m) - cortical surface analyses (SurfStat)

## Results
The vertex-wise results are provided in separate files for the four cognitive domains and for 'g'. These files contain t-stats, p-values, and FDR-corrected q-values. 

## Tools Used
Statistical analyses were conducted using R v4.1.0, and CFA models were constructed using R's 'lavaan' package.

Cortical surface computations were performed using SurfStat, a MATLAB toolbox for the statistical analysis of surface and volumetric data. SurfStat is available at [https://www.math.mcgill.ca/keith/surfstat/](https://www.math.mcgill.ca/keith/surfstat/).

## Authors
These scripts were authored by Danielle Page. Minor revisions were made by Joanna Moodie and Colin Buchanan.
