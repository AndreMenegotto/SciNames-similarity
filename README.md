Similarity among scientific names
================

André Menegotto

<br>

## Introduction

<div align="justify">

This repository contains the code used to quantify orthographic
similarity among species names using the Damerau-Levenshtein distance,
which measures the number of edit operations required to transform one
string into another. In the script *`1_NameStats.R`*, orthographic
distances are compared across complete scientific names, as well as
separately for genera among species sharing the same epithet and for
epithets among species sharing the same genus. The frequency and
taxonomic distribution of commonly shared epithets are also examined.

Analyses are based on accepted species names from the World Register of
Marine Species (WoRMS). Raw data are not included in this repository but
can be accessed via the sources referenced within the scripts. All
processed datasets used to generate the results of the study [How
similar are species names and why does this matter for biodiversity
data](https://doi.org/10.3897/arphapreprints.e196971) are provided as
RData files in the *`./Data/tmpRData`* directory.

</div>
