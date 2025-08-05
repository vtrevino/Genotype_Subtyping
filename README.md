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


## RUN
The GS_RESULT within the GS list object will contain all combinations independently of whether the sub-GWAS was performed or not. Therefore, only those rows whose GS_P is not NA have been estimated.

Code in R for estimation:
``````
   GS$GS_RESULT <- c()
   MinSamplesPerGT <- 500
   ntot <- length(GS$CONTROLS) + length(GS$CASES)
   ncalc <- 0
   zeroMatTab <- matrix(0,nrow=2,ncol=3,dimnames=list(c("0","1"),c("0","1","2")))

   for (i in 1:nrow(GS$GWAS_SEED)) {
       SNPi <- GS$GWAS_SEED$SNP[i]
       geneResDF <- c()
       for (gt in 0:2) {
           WGi <- which(GS$GENOTYPES[SNPi,] == gt)
           nGTi <- length(WGi) #sum(GS$GENOTYPES[SNPi,] == gt, na.rm=TRUE)
           colsGTi <- names(WGi)
           for (j in 1:nrow(GS$GWAS_TO_GS)) {
               SNPj <- GS$GWAS_TO_GS$SNP[j]
               if (SNPi != SNPj) {
                   if (j %% ifelse(nGTi > 5000, 20, 100) == 1) {
                       cat("i=",i," of ",nrow(GS$GWAS_SEED),", SNP GWAS:",SNPi,", GT:",gt," (",nGTi,"), j=",j,", SNP GS:",SNPj,", nCalc:",ncalc,"\n",sep="")
                       flush.console()
                   }
                   gtMatrix <- cbind(
                       Phenotype=ifelse(colsGTi %in% GS$CONTROLS, 0, 1), # 0=Control, 1=Case
                       Genotype=GS$GENOTYPES[SNPj, colsGTi])
                   beta <- p <- fit <- intp <- NA
                   if (nGTi >= MinSamplesPerGT && ntot-nGTi >= MinSamplesPerGT) {
                       ncalc = ncalc + 1
                       gtMatrix <- cbind(gtMatrix,
                           Co=GS$COVARIATES[colsGTi, ])
                       model <- glm(Phenotype ~ ., family = "binomial", data = gtMatrix, maxit = 500)
                       beta <- model$coefficients["Genotype"]
                       fit <- model$converged
                       if (!is.na(beta)) {
                           p <- summary(model)$coefficients["Genotype",4]
   
                           if (p < 5e-7) {
                               # Test for Interaction (in all samples)
                               cols <- c(GS$CONTROLS, GS$CASES)
                               intMatrix <- cbind(
                                   Phenotype=ifelse(cols %in% GS$CONTROLS, 0, 1), # 0=Control, 1=Case
                                   GenotypeI=GS$GENOTYPES[SNPi, cols], 
                                   GenotypeJ=GS$GENOTYPES[SNPj, cols], 
                                   Co=GS$COVARIATES[cols, ])
                               intModel <- glm(Phenotype ~ . + GenotypeI:GenotypeJ, family = "binomial", data = intMatrix, maxit = 500)
                               intBeta <- intModel$coefficients["GenotypeI:GenotypeJ"]
                               if (!is.na(intBeta)) {
                                   intp <- summary(intModel)$coefficients["GenotypeI:GenotypeJ",4]
                               }
                           }
                       }
                   }
                   # Store results
                   xT <- table(gtMatrix[,"Phenotype"], gtMatrix[,"Genotype"])
                   if (nrow(xT) != 2 || ncol(xT) != 3) {
                       xMT <- zeroMatTab
                       if (nGTi > 0) for (lei in 1:nrow(xT)) xMT[rownames(xT)[lei],colnames(xT)] <- xT[lei,] 
                       xT <- xMT
                   }
                   geneResDF <- rbind(geneResDF,
                     data.frame(I=i, J=j, SNP_GWAS = SNPi, BETA_GWAS = GS$GWAS_SEED$BETA[i],
                       P_GWAS = GS$GWAS_SEED$P[i], GENOTYPES_GWAS = paste0("G",gt), N_GT_GWAS = nGTi,
                       GS_SNP = SNPj, GS_N = nGTi, 
                       GS_N_CASES = sum(xT["1",],na.rm=TRUE),          #sum(Xpt == 1, na.rm=TRUE),
                       GS_N_CONTROLS = sum(xT["0",],na.rm=TRUE),       #sum(Xpt == 0, na.rm=TRUE),
                       GS_N_CASES_AA_0 = sum(xT["1","0"],na.rm=TRUE),  #sum(Xpt == 1 & Xgt == 0, na.rm=TRUE),
                       GS_N_CASES_AB_1 = sum(xT["1","1"],na.rm=TRUE),  #sum(Xpt == 1 & Xgt == 1, na.rm=TRUE),
                       GS_N_CASES_BB_2 = sum(xT["1","2"],na.rm=TRUE),  #sum(Xpt == 1 & Xgt == 2, na.rm=TRUE),
                       GS_N_CONTROLS_AA_0 = sum(xT["0","0"],na.rm=TRUE),  #sum(Xpt == 0 & Xgt == 0, na.rm=TRUE),
                       GS_N_CONTROLS_AB_1 = sum(xT["0","1"],na.rm=TRUE),  #sum(Xpt == 0 & Xgt == 1, na.rm=TRUE),
                       GS_N_CONTROLS_BB_2 = sum(xT["0","2"],na.rm=TRUE),  #sum(Xpt == 0 & Xgt == 2, na.rm=TRUE),
                       GS_BETA = beta, GS_P = p, GS_CONVERGED = fit,
                       GS_P_INT = intp,
                       snp_GWAS_BETA = GS$GWAS_TO_GS$BETA[j], snp_GWAS_P = GS$GWAS_TO_GS$P[j],
                       stringsAsFactors=FALSE
                       ))
                 }
           } 
       }
       GS$GS_RESULT <- rbind(GS$GS_RESULT, geneResDF)
   }

   GS$GS_RESULT$DELTA_P <- log10(GS$GS_RESULT$snp_GWAS_P)-log10(GS$GS_RESULT$GS_P)
   GS$GS_RESULT$GS_CASES_Ns <- with(GS$GS_RESULT, paste(GS_N_CASES,GS_N_CASES_AA_0,GS_N_CASES_AB_1,GS_N_CASES_BB_2,sep="|"))
   GS$GS_RESULT$GS_CONTROLS_Ns <- with(GS$GS_RESULT, paste(GS_N_CONTROLS,GS_N_CONTROLS_AA_0,GS_N_CONTROLS_AB_1,GS_N_CONTROLS_BB_2,sep="|"))
   GS$GS_RESULT$GS_CASES_Rs <- with(GS$GS_RESULT, paste(GS_N_CASES,round(GS_N_CASES_AA_0/GS_N_CASES,3),round(GS_N_CASES_AB_1/GS_N_CASES,3),round(GS_N_CASES_BB_2/GS_N_CASES,3),sep="|"))
   GS$GS_RESULT$GS_CONTROLS_Rs <- with(GS$GS_RESULT, paste(GS_N_CONTROLS,round(GS_N_CONTROLS_AA_0/GS_N_CONTROLS,3),round(GS_N_CONTROLS_AB_1/GS_N_CONTROLS,3),round(GS_N_CONTROLS_BB_2/GS_N_CONTROLS,3),sep="|"))
   
   save(GS, file="<your-file>.Rdata")
   
   subset(GS$GS_RESULT, GS_P < 5e-7 & DELTA_P > 1 & snp_GWAS_P > 5e-8)
``````
