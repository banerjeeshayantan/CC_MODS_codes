#THERE ARE TWO SETS OF FEATURES TESTED ACROSS ALL TRAIN-VALIDATION-TEST DATASETS - mods_predictors and combined_mortality_predictors. Users can change the code to accomodate the changes in the analyses accordingly.

#..................................PREPARING GSE66099.........................................
#gse66099gpl570GEM.RDS, gse66099gpl570KEY.RDS, mortality_predictors.rds are available from https://www.synapse.org/#!Synapse:syn5612563
#preparing GSE66099 data with mortality labels
gse66099_processed_data = readRDS("gse66099gpl570GEM.RDS")
gse66099_gene_key = readRDS("gse66099gpl570KEY.RDS")
gse66099_gene_key = as.data.frame(gse66099_gene_key)
colnames(gse66099_gene_key) = "Symbols"
combined_mortality_predictors = readRDS("mortality_predictors.rds")
gse66099_processed_data$Genes = gse66099_gene_key[match(rownames(gse66099_processed_data), rownames(gse66099_gene_key)), "Symbols"]
setdiff(combined_mortality_predictors, gse66099_processed_data$Genes) #OR52R1 not present in the gene list (~22000 genes) belonging to the GPL70 platform for GSE66099 (train set). hence omitting this.

combined_mortality_predictors = combined_mortality_predictors[-c(which(combined_mortality_predictors == "OR52R1"))]
length(combined_mortality_predictors)
gse66099_processed_data_mortality_genes = gse66099_processed_data[which(gse66099_processed_data$Genes %in% combined_mortality_predictors),]
dim(gse66099_processed_data_mortality_genes) #110 genes by 363 patients


#preparing GSE66099 data with MODS labels
gse66099_mods_data = read.delim("DGE_RMA_MODS_day3day7.txt")
mods_predictors = read.delim("Top_feats_Day3Day7MODS_without_propensity.txt")
#subsetting gse66099 using samples from DGE_RMA
#MODS_pheno should be in your environment
gse66099_clinical = MODS_pheno
gse66099_clinical = subset(gse66099_clinical, gse66099_clinical$Day.3.MODS.Trajectory_2_Day.7.MODS..1_No.Day.3.MODs.and.Yes.Day.7.MODS...0_No.Day.7.MODS != 1)
gse66099_clinical$sample_id = unlist(lapply(gse66099_clinical$Microarray.ID, function(x) strsplit(strsplit(x,"_")[[1]][2], ".CEL")))
newrma_aggregrated=aggregate(gse66099_mods_data,
                             by = list(Gene = mods_predictors$X),
                             FUN = mean,
                             na.rm = TRUE)
rownames(newrma_aggregrated) = newrma_aggregrated$Gene
newrma_aggregrated$Gene = NULL
newrma_aggregrated = as.data.frame(t(newrma_aggregrated))
newrma_aggregrated$Labels = gse66099_clinical[match(rownames(newrma_aggregrated), gse66099_clinical$sample_id), "Day.3.MODS.Trajectory_2_Day.7.MODS..1_No.Day.3.MODs.and.Yes.Day.7.MODS...0_No.Day.7.MODS"]
newrma_aggregrated = newrma_aggregrated[-c(which(is.na(newrma_aggregrated$Labels))), ]


#..................................PREPARING E_MTAB-5882.........................................
newrma_aggregrated_5882 = read.delim("External_Validation_Data_Day3Day7_5882.txt")
rownames(newrma_aggregrated_5882) = paste("X", rownames(newrma_aggregrated_5882), sep="")
emtab_5882 = read.delim("AE_normalised_counts_Cabrera_et_al_5_12_13.txt")
emtab_5882_84_patients = emtab_5882[, which(colnames(emtab_5882) %in% paste("X", rownames(newrma_aggregrated_5882), sep=""))]
emtab_5882_84_patients$ProbeID = emtab_5882$Probe_Id
emtab_5882_probeids = read.delim("A-GEOD-10558.adf.txt")
emtab_5882_84_patients$Gene_name = emtab_5882_probeids[match(emtab_5882_84_patients$ProbeID, emtab_5882_probeids$Reporter.Name), 6]
setdiff(mods_predictors$X, emtab_5882_84_patients$Gene_name) #use combined_mortality_predictors instead of mods_predictors while comparing model performances using mods-based features vs. mortality-based features
emtab_5882_84_patients[which(emtab_5882_84_patients$Gene_name == "C11ORF74"), "Gene_name"] <- "C11orf74"
setdiff(mods_predictors$X, emtab_5882_84_patients$Gene_name)
mods_predictors_1 = mods_predictors$X[-c(which(mods_predictors$X %in% setdiff(mods_predictors$X, emtab_5882_84_patients$Gene_name)
))]
emtab_5882_84_patients = emtab_5882_84_patients[which(emtab_5882_84_patients$Gene_name %in% mods_predictors_1), ]
emtab_5882_84_patients$ProbeID = NULL
newrma_aggregrated_5882_new=aggregate(emtab_5882_84_patients[, -c(85)],
                                      by = list(Gene = emtab_5882_84_patients$Gene_name),
                                      FUN = mean,
                                      na.rm = TRUE)
str(newrma_aggregrated_5882_new)
rownames(newrma_aggregrated_5882_new) = newrma_aggregrated_5882_new$Gene
newrma_aggregrated_5882_new$Gene = NULL
newrma_aggregrated_5882_new_t = as.data.frame(t(newrma_aggregrated_5882_new))
newrma_aggregrated_5882_new_t$Label = newrma_aggregrated_5882[match(rownames(newrma_aggregrated_5882_new_t), rownames(newrma_aggregrated_5882)), "Label"]


#..................................PREPARING E-MTAB-10938.........................................
emtab_10938 = readRDS("EMTAB_10938.RDS")
setdiff(mods_predictors$X, emtab_10938$ID)
ext_vali_data_10938_Proulx = read.delim("ext_vali_data_10938_Proulx.txt")
samples = rownames(ext_vali_data_10938_Proulx)
colnames(emtab_10938)[grep("Norm_72.01", colnames(emtab_10938))] = "Norm_72.01"
emtab_10938_32_patients = emtab_10938[, which(colnames(emtab_10938) %in% samples)]
emtab_10938_32_patients$Gene_name = emtab_10938$ID
emtab_10938_32_patients_52_genes = emtab_10938_32_patients[which(emtab_10938_32_patients$Gene_name %in% mods_predictors$X), ]
rownames(emtab_10938_32_patients_52_genes) = emtab_10938_32_patients_52_genes$Gene_name
emtab_10938_32_patients_52_genes$Gene_name = NULL
newrma_aggregrated_10938 = as.data.frame(t(emtab_10938_32_patients_52_genes))
dim(newrma_aggregrated_10938)
newrma_aggregrated_10938_new = newrma_aggregrated_10938[samples,]
newrma_aggregrated_10938_new$Label = ext_vali_data_10938_Proulx$Label


#..................................PREPARING E-MTAB-1548.........................................
emtab_1548 = readRDS("EMTAB1548GEM.RDS")
emtab_1548_key = readRDS("EMTAB1548KEY.RDS")
emtab_mods = read.delim("E-MTAB-1548-all_genes.txt")
dim(emtab_mods)
emtab_1548_key = as.data.frame(tibble::enframe(emtab_1548_key))
colnames(emtab_1548_key) = c("ProbeID", "Gene")
emtab_1548 = as.data.frame(emtab_1548)
emtab_1548$Probe = rownames(emtab_1548)
emtab_1548$Gene = emtab_1548_key[match(emtab_1548$Probe, emtab_1548_key$ProbeID), "Gene"]
dim(emtab_1548)
emtab_1548_new = emtab_1548[which(emtab_1548$Gene %in% mods_predictors$X), ]
newrma_aggregrated_1548_new=aggregate(emtab_1548_new[, -c(75,76)],
                                      by = list(Gene = emtab_1548_new$Gene),
                                      FUN = mean,
                                      na.rm = TRUE)
rownames(newrma_aggregrated_1548_new) = newrma_aggregrated_1548_new$Gene
newrma_aggregrated_1548_new$Gene = NULL

dim(newrma_aggregrated_1548_new_t)
newrma_aggregrated_1548_new_t = as.data.frame(t(newrma_aggregrated_1548_new))
newrma_aggregrated_1548_new_t$Label = emtab_mods$Mortality
write.table(newrma_aggregrated_1548_new_t, "E-MTAB_1548_MODS_data.txt", col.names = T, row.names = T, sep="\t", quote = F)


#..................................PREPARING GSE144406.........................................
GSE144406_expression_data_0.96 <- read.delim("GSE144406_expression_data_0.96.txt")
gse144406 = read.csv("GSE144406_All_samples_TPM_count.csv", header = T)
library(org.Hs.eg.db)
annots <- select(org.Hs.eg.db, keys=gse144406$X, 
                 columns="SYMBOL", keytype="ENSEMBL")
gse144406$gene_symbol = annots$SYMBOL
annots_sub = annots[which(annots$SYMBOL %in% mods_predictors),]
gse144406_new = gse144406[which(gse144406$X %in% annots_sub$ENSEMBL), ]
dim(gse144406_new)
gse144406_new$gene_symbol = annots_sub$SYMBOL
gse144406_new$X = NULL
rownames(gse144406_new) = gse144406_new$gene_symbol
gse144406_new$gene_symbol = NULL
newrma_aggregrated_144406 = as.data.frame(t(gse144406_new))
newrma_aggregrated_144406$Label = GSE144406_expression_data_0.96$Label