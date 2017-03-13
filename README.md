# divBCR

**divBCR.R** and **divBCR_isotype.R** scripts calculate different biodiversity metrics, as well as metrics that characterize the different aspects of the BCR repertoire. This script depends on the R package “plyr”.

Usage: Rscript divBCR.R inputFile.txt outputFile.txt

The input file is a processed aggregated file that will need to be produced from the V'DJer, so that the columns are in the following order: sample, sequence, cdr3, expected_counts, seq_id, isotype, vregion_identity, aa_cdr3, vgene, jgene, total_count, cluster

Overview:

Step 1) B cell recptor counts are normalized to total RNA-seq counts

Step 2) Calculations of biodiversity and other metrics at the sequence leve

Step 3) Sums clusters

Step 4) Calculations of biodiversity at the cluster level (in case of poor sequence quality)

Step 5) Executing the functions, merging the two data frames and saving the output to a file.


Below are the calculations for each sample (not in column order in final output, column header listed, if by* then calculated both at the sequence and cluster level):

Biodiversity metrics:

1.	Shannon entropy (shannon_entropy_by*): -Σ (count proportion * log(count proportion))

2.	Evenness (evenness_by*): Shannon entropy / log(species richness)

3. True Diversity (true_div_by*): Exp(Shannon entropy)

4. Gini-Simpson (gini_simpson_by*): 1-(sum(count proportion)^2)

BCR repertoire characteristics:

5.	Top clone proportion (top_clone_prop_by*):  max(“species” count) / total counts

6. Second top clone proportion (second_top_clone_by*):  second max (“species” count) / total counts

7. Number of unique sequences or clusters (num_unique*): species richness, number of unique sequences or number of unique clusters

8. Total counts (total_counts): sum(count of each sequence)

9.	 Average CDR3 length (mean_cdr3_aa_len): mean(CDR3 amino acid length)

10.	 Weighted average CDR3 length (w_mean_cdr3_aa_len):  sum(CDR3 length* count of each sequence) / sum(count of each sequence) 

11. Mean V-region identity (mean_v_region_ident): mean(V-region identity)

12.	 Weighted mean V-region identity (w_mean_v_region_ident):  sum(V-region identity* count of each sequence) / sum(count of each sequence)
