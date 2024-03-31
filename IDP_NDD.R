##############   load packages ########################
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(MRPRESSO)
library(dplyr)
library(mr.raps)
##############   Extracting significant SNPs associated with exposure ########################

exposure_name <- ''#need modification

exp_data <- fread(exposure_name,header=T)
expo_data <- subset(exp_data,exp_dat$pval < 5e-8)
expo_clump_data <- clump_data(expo_data,clump_r2=0.001,clump_kb=10000)

filename <- ''#need modification
write.csv(expo_clump_data, file = filename,row.names =FALSE) 

##############   Extracting information of IVs in outcome #######################

exposure_filename <- ''#need modification
outcome_name <- ''#need modification

outcome_GWAS <- fread(outcome_name,header =TRUE,sep = " ")
exposure <- read.csv(exposure_filename,header =TRUE)
exposure_outcome <- format_data(dat=outcome_GWAS,
                                type = "outcome",#need modification
                                snps = exposure$SNP,#need modification
                                header = TRUE,
                                phenotype_col = "phenotype",#need modification
                                snp_col = "SNP",#need modification
                                beta_col = "BETA",#need modification
                                se_col = "SE",#need modification
                                eaf_col = "EAF",#need modification
                                effect_allele_col = "EA",#need modification
                                other_allele_col = "OA",#need modification
                                pval_col = "P",#need modification
                                chr_col = "CHR",#need modification
                                pos_col = "POS",#need modification
                                samplesize_col = "N"#need modification
                                )

filename <- ''#need modification
write.csv(exposure_outcome, file = filename,row.names =FALSE) 

##############  MR analyses  ##########################

exposure_filename <- ''#need modification
exposure_outcome_filename <- ''#need modification
outSNP_name <- ''#need modification
confounder_name <- ''#need modification

confounder <- read.table(confounder_name,header = T)
outSNP <- read.csv(outSNP_name,header = T)

exposure <- read_exposure_data(filename =exposure_filename,clump = FALSE,sep= ",",
                                    snp_col = "SNP",
                                    beta_col = "beta.exposure",
                                    se_col = "se.exposure",
                                    effect_allele_col ="effect_allele.exposure",
                                    other_allele_col = "other_allele.exposure",
                                    eaf_col = "eaf.exposure",
                                    pval_col = "pval.exposure",
                                    phenotype_col = "exposure",
                                    samplesize_col = "samplesize.exposure",
                                    chr_col = "chr.exposure",
                                    pos_col = "pos.exposure")

exposure_outcome <- read.csv(exposure_outcome_filename,header =T)

mydata <- harmonise_data(exposure_dat = exposure,outcome_dat = exposure_outcome,action = 2)
  
  
if (length(intersect(mydata$SNP,outSNP$SNP)) >0){
   mydata <- anti_join(mydata,outSNP,by='SNP')
   }
        
if (length(intersect(mydata$SNP,confounder$SNP)) >0){
    mydata <- anti_join(mydata,confounder,by='SNP')
    }

res <- mr(mydata,method_list=c( 'mr_wald_ratio','mr_two_sample_ml','mr_egger_regression','mr_weighted_median','mr_ivw'))

mr_sen_raps <- data.frame(mr.raps(mydata$beta.exposure,mydata$beta.outcome,mydata$se.exposure,mydata$se.outcome))
          
Inter <- mr_egger_regression(mydata$beta.exposure, mydata$beta.outcome, mydata$se.exposure, mydata$se.outcome)
         
het <- mr_heterogeneity(mydata)

pleio <- mr_pleiotropy_test(mydata)
  
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                                OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                                data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
            
