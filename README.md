# Species delimitation in an intractable syngameon: Bringing order to the polyphyletic *Heuchera americana* group.

## Morphometrics
* This directory contains scripts, data, and results for morphometric analyses.
- `final_accession_table.xlsx`: Finalized table of specimen accessions and associated data used in the taxonomic treatment.
- `Accession_table_speciespaperversion.docx`: Manuscript version of the accession table formatted for publication.
- `morphometric_data_final.csv`: Final data file used for the morphometric analyses.
- `morphometrics_updated_final.r`: R script used for processing and analyzing the morphometric data.
- `results`: Directory containing the output and figures from the morphometric analyses.

## Phylogenetics
* This directory contains subdirectories pertaining to various genetic analyses. These include: 
- `ASTRAL_pure`: Contains scripts and data files for the inference of an [ASTRAL](https://github.com/smirarab/ASTRAL) phylogeny using individuals assigned as "pure" (taxon assignment > 90%) by [fastSTRUCTURE](https://github.com/rajanil/fastStructure), and a `results` subfolder.
    - `ASTRAL_pure.sh`: Shell script for constructing an [ASTRAL](https://github.com/smirarab/ASTRAL) phylogeny.
    - `all_gene_trees.pure.bs10.tre`: Input file containing individual gene trees (filtered for >10% bootstrap support) used for species tree estimation.
    - `renaming_script_Sept2023.sh`: Script used to batch-rename tips to hypothesized species limits.
    - `new_names_putative_species.txt` and `old_names.txt`: Mapping files used by the renaming script to transition from field IDs/older nomenclature to putative species names.
    - `results`: Subfolder containing the final ASTRAL species tree and figures.
- `bpp_pure`: Contains scripts, control files, and results for Bayesian species delimitation using [BPP](https://github.com/bpp/bpp). 
    - `americana_all_pure.ctl`: Control file specifying parameters for the BPP analysis.
    - `alignments_reduced.phy.zip`: Compressed phylogenomic alignments.
    - `pure_subset.sh` and `summarize_bpp_mcmc.sh`: Shell scripts for subsetting data and summarizing MCMC results.
    - `guide_tree.tre`: The starting guide tree.
	- `lmap.txt` and `lmap_reduced.txt`: Taxon mapping files assigning individual specimens (e.g., A1-1) to hypothesized species groups (e.g., CALYCOSA, BREVIPETALA) for BPP.
    - `pure_individuals.txt` and `mixed_individuals.txt`: Lists of taxa categorized by purity threshold (e.g., "pure" (taxon assignment > 90%) assigned by [fastSTRUCTURE](https://github.com/rajanil/fastStructure)).
    - `outgroup_list.txt`: List of specimen IDs used as outgroups for rooting the phylogeny.
    - `phyx.logfile`: Log file documenting the data processing steps performed with phyx.
    - `best_tree.tre.pdf`, `.tre`, `.png`, and `.svg`: Rendered figures of the BPP species tree results.
- `pop_gen`: Contains the bioinformatic pipeline for variant calling, population structure, and diversity statistics.
    - `1_gatk_variantcall.sh` and `2_gatk_combinefilter.sh`: Initial shell scripts for processing raw reads into a filtered, multi-sample VCF.
    - `3_faststructure.sh`: Script for running [fastSTRUCTURE](https://github.com/rajanil/fastStructure) to infer ancestral clusters.
    - `4_exploratoryPCA.sh`: Principal Component Analysis (PCA) to visualize genetic clusters.
    - `5_popgen_stats.sh`: R script to calculate population-level statistics, including:
        - Pairwise $F_{ST}$ (Weir and Cockerham, 1984; Nei, 1987) with heatmap.
        - Allelic richness ($A_r$), observed heterozygosity ($H_o$), and within-population gene diversity ($H_s$).
        - Wright’s inbreeding coefficient ($F_{IS}$).
    - `americana.SNPall.pruned.vcf.recode.vcf`: The final pruned VCF file used for all downstream population genetic analyses.
    - `americanapairwiseWCfst.csv`: Output matrix of pairwise Weir and Cockerham $F_{ST}$ values.
    - `americana_popfile.txt`, `pca_popfile.txt`, and `sample_list.txt`: Coordinate and mapping files used for assigning individuals to populations and color-coding figures.
    - `new_names_putative_species.txt`: Mapping file to transition from internal MSU specimen IDs to formal species names.
    - `results`: Subfolder containing results.

---
*Markdown documentation for this repository was assisted by Gemini (Google AI).*
