# Filter Faulty sgRNA

### Description
We designed multiple sgRNA for each target site to account for off-target/non-functional sgRNA. This approach compares the distribution of all sgRNA for a given target across
the CM subpopulations from our large screen. sgRNA that have a bias in distribution in comparison to the rest for the same target are excluded from downstream analaysis.

### Basic Overview
* [Part 1: Assigning pvalues for clustering bias](https://github.com/darmen04/Repression-of-CHD-associated-enhancers-delays-human-cardiomyocyte-lineage-commitment/blob/main/Notebooks/sgRNA_Filtering/sgRNA_PVAL_Combos.ipynb)
1. We create a background of control cell distribution using our 5 NC. This is done by using hypergeometric tests to compare distribution of 1 NC to the rest minus 1 across all clusters. The NC which when removed prevents bias of the remaining NCs in all clusters was removed. This creaed the control cell background.
2. For each target, we ran a hypergeometric test comparing the distribution of target sgRNA to the filtered control background. This was done for all possible combinations of guides, across all CM subpopulations, and for both enrichment/depletion.

* [Part 2: Assigning pvalues for clustering bias](https://github.com/darmen04/Repression-of-CHD-associated-enhancers-delays-human-cardiomyocyte-lineage-commitment/blob/main/Notebooks/sgRNA_Filtering/Filter_sgRNA.ipynb)
3. We took each set of pvalues and ordered them from most to least significant. A cutoff was drawn at the halfway mark corresponding to the number of times each guide appears in a combination.
4. We ran a hypergeometric test to determine if a given guide appeared more frequently to the left of the cutoff than the rest. 
5. All sgRNA which were found to be uniquely enriched to the left for a given test were removed.
