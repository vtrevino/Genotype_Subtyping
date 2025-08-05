# Genotype_Subtyping
Genotype Subtyping Approach to Identify Unnoticed Variants in Diseases from GWAS data.

We propose a method referred to as ‘Genotype Subtyping’ that is based on a classical GWAS to identify seed SNPs. Each SNP stratifies cases and controls by its genotype, followed by a sub-GWAS for each genotype, while keeping the number of additional statistical tests low.

Thus, the steps are:

1) Perform a classical GWAS (case-control).
   
2) Obtain Seed-SNP used to stratify:
  - SNPs whose P < 1e-6 (threshold can be changed)
  - A LD clumping among these SNPs is recommended
    
3) Obtain SNPs to be re-tested: 
  - SNPs whose P < 1e-3 (threshold can be changed)
    
4) For each Seed-SNP and each of its 3 genotypes (AA, AB, BB)
  - Obtain the subset of controls whose Seed-SNP == genotype
  - Obtain the subset of cases whose Seed-SNP == genotype
  - Perform a sub-GWAS with only those SNPs from step 3 and the subsets of cases and controls
  - Aggregate results

5) Analyze/highlight those SNPs whose sub-GWAS-P is more significant than the original GWAS-P and whose sub-GWAS-P is < 5e-8.


