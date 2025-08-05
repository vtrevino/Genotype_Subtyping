# Genotype_Subtyping
Genotype Subtyping Approach to Identify Unnoticed Variants in Diseases from GWAS data.

We propose a method referred to as ‘Genotype Subtyping’ that is based on a classical GWAS to identify seed SNPs. Each SNP stratifies cases and controls by its genotype, followed by a sub-GWAS for each genotype. The number of additional statistical tests is kept low to avoid multiplicity test issues.

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


# Estimation of Genotype Subtyping Statistics
The code and description below are based on the assumption that a classical GWAS has been performed.

In R, create a list object named GS. The GS object must contain the following items.
- CASES     : Vector of sample names for cases
- CONTROLS  : Vector of sample names for controls
- GWAS_SEED : Data frame of SEED snps. Columns from GWAS named as : SNP, BETA, P (typically those SNPs whose P < 1e-6 + LD clumping).
- GWAS_TO_GS: Data frame of snps that will be partitioned by GWAS_SEED genotypes (typically those SNPs whose P < 1e-3). Columns from GWAS named as : SNP, BETA, P
- GENOTYPES : Matrix of SNP (rows) and Samples (columns) that contains at least all SNPs from GWAS_SEED and GWAS_TO_GS.
- COVARIATES: Data Frame rownames=Samples, columns=Covariates. NO other columns are allowed (be careful with columns).
- GS_RESULT : Data frame of results.

Using the below R code, the following columns in GS$GS_RESULT will be created:
- SNP_GWAS : The (seed) SNP used to stratify.
- BETA_GWAS : The beta value originally assigned to the corresponding seed SNP.
- P_GWAS : The p value originally assigned to the corresponding seed SNP.
- GENOTYPE_GWAS : The genotype of SNP_GWAS analyzed.
- GS_SNP : The SNP estimated by the sub-GWAS.
- GS_N : The Total number of samples used in the sub-GWAS.
- GS_N_CASES: The total number of cases used in the sub-GWAS.
- GS_N_CONTROLS: The total number of controls used in the sub-GWAS.
- GS_N_CASES_AA: The cases used in the sub-GWAS whose genotype is AA.
- GS_N_CASES_AB: The cases used in the sub-GWAS whose genotype is AB.
- GS_N_CASES_BB: The cases used in the sub-GWAS whose genotype is BB.
- GS_N_CONTROL_AA: The controls used in the sub-GWAS whose genotype is AA.
- GS_N_CONTROL_AB: The controls used in the sub-GWAS whose genotype is AB.
- GS_N_CONTROL_BB: The controls used in the sub-GWAS whose genotype is BB.
- GS_BETA: The estimated beta value for the GS_SNP when SNP_GWAS=GENOTYPE_GWAS.
- GS_P: The p-value for the GS_SNP when SNP_GWAS=GENOTYPE_GWAS.
- GS_CONVERGED: If the glm converged (True/false).
- GS_P_INT: p-value of the interaction of SNP_GWAS and GS_SNP (only estimated when GS_P is < 5e-7)
- snp_GWAS_BETA: The beta value of the GS_SNP from the original GWAS.
- snp_GWAS_P: The P value of the GS_SNP from the original GWAS.
- GS_DELTA_P: The orders of magnitude in p-value gained in the sub-gwas (rare to be positive, those that are interesting).
- GS_CASES_Ns: A short representation of cases within the sub-gwas.
- GS_CONTROLS_Ns: A short representation of controls within the sub-gwas.
- GS_CASES_Rs: A short representation of the fraction of cases within the sub-gwas.
- GS_CONTROLS_Rs: A short representation of the fraction of controls within the sub-gwas.
