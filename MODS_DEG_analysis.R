library(GEOquery)
library(affy)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(sva)
library(limma)
library(RColorBrewer)
library(annotate)
library(hgu133plus2.db)
library(calibrate)
library(affyio)
#setwd to the folder titled data containing the .CEL files. This large folder has ~200 files
setwd("~/CC_MODS_codes_github/data")
#subset the new phenotypic information stored in the variable MODS_pheno. This step assumes that you already have the phenotypic information for the MODS data
tab=data.frame(Samples=MODS_pheno$Microarray.ID, 
               Outcome=MODS_pheno$Day.3.MODS.Trajectory_2_Day.7.MODS..1_No.Day.3.MODs.and.Yes.Day.7.MODS...0_No.Day.7.MODS)
rownames(tab)=tab$Samples
tab[,1] %in% fns
fns <- list.celfiles()
fns
#use the get.celfile.dates() from the affyio package to get the dates and 
#extract all the dates for 201 patients
#adjusting for batch effects arising due to differences in collection times
dates=unlist(lapply(get.celfile.dates(fns),function(x) substr(x,1,4)))
batch_effect_file=data.frame(Samples=fns,Batch_Variable=dates)
pd1=pd[which(pd$Samples %in% tab$Samples),]
head(pd1)
pd2=tab[match(pd1$Samples, tab$Samples),c("Outcome","Samples")]
pd2$Batch=pd1$Batch_Variable

#reading CEL files and normalizing using rma method
rawData=ReadAffy(verbose=TRUE, filenames=tab[,1], cdfname="hgu133plus2",phenoData=tab)
print("rawData completed")
par(mfrow=c(2,1))
boxplot(rawData,col="red",main="Before Normalization",names = FALSE)
normData <- expresso(rawData, normalize.method="quantiles", bgcorrect.method="rma", summary.method="medianpolish", pmcorrect.method="pmonly")
rma=exprs(normData)
boxplot(rma,col="blue", main="After Normalization",names = FALSE)
rma=rma[which(!grepl("AFFX", rownames(rma))),]
print("AFFX")
#Format values to 5 decimal places
rma=format(rma, digits=5)
print("format rma")
#Map probe set identifiers to Entrez gene symbols and IDs and then combine with raw data.
#Map probe sets to gene symbols or other annotations
#To see all available mappings for this platform
ls("package:hgu133plus2.db") #Annotations at the exon probeset level

probes=row.names(rma)
#probes
#Extract probe ids, entrez symbols, and entrez ids
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
#Symbols
print("Symbols")
#Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))
#hgu133a2SYMBOL
#Combine gene annotations with raw data
rma=data.frame(probes,Symbols,rma)


#Write RMA-normalized, mapped data to file
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#BiocManager::install("limma")
library(limma)

table <-read.table("rma.txt",header=TRUE,sep="\t")
#head(table)
#comb2 <- subset(table, select = -c(Symbols))
comb2 <- subset(rma, select = -c(probes)) #PUT IT BACK FOR GENE LEVEL  b bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
comb2 <- comb2[complete.cases(comb2[ , 1]),]
symbols=comb2$Symbols
comb2$Symbols=NULL
#rownames(comb2)=comb2$Symbols
str(comb2)
print("comb2_1")
comb2 <- data.frame(lapply(comb2, as.character), stringsAsFactors = FALSE)
comb2 <- data.frame(lapply(comb2, as.numeric), stringsAsFactors = FALSE)
#comb2 <- comb2[complete.cases(comb2[ , 1]),]
#head(comb2)
#print("comb2_2")
comb2[1:10,1:10]
comb2$Symbols=symbols
out <- aggregate(. ~ Symbols, comb2, mean) #PUT IT BACK FOR GENE LEVEL bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
#out <- comb2
symb=out$Symbols
head(out)
# Format values to 5 decimal places
out = format(out, digits = 5)
out[1:3, 1:3]
rownames(out)=out$Symbols
out$Symbols=NULL
write.table(out, file = "DGE_RMA_MODS_day3day7.txt", quote = FALSE, sep = "\t", row.names = FALSE)
out <-read.table("DGE.RMA.txt",header=TRUE,sep="\t")
dim(comb2)
dim(out)
out <- lapply(out, function(x) as.numeric(as.character(x)))
out=as.data.frame(out)
rownames(out)=symb


#Batch correction and DEG analysis

genes <- out$Symbols
Batch <- pd2$Batch
mod = model.matrix(~as.factor(pd2$Outcome) + as.factor(Batch),data=pd2)
rownames(out)=out$Symbols
#out <- out[,-1]
out <- data.frame(lapply(out, as.numeric), stringsAsFactors = FALSE)
fit = lm.fit(mod, t(out))
hist(fit$coefficients[2,],col=2,breaks=100)
dev.off()
table(pd2$Outcome,Batch)
modcombat = model.matrix(~1, data=pd2)
modstatus = model.matrix(~pd2$Outcome, data=pd2)
combat_out = ComBat(dat=as.matrix(out), batch=Batch,mod=modcombat, par.prior=TRUE, prior.plot=FALSE)
combat_fit = lm.fit(modstatus,t(combat_out))
dim(t(combat_out))
combat_fit1 = lmFit(combat_out,modstatus)
hist(combat_fit$coefficients[2,],col=2,breaks=100)
plot(fit$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="Linear Model",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)
mod = model.matrix(~pd2$Outcome, data=pd2)
mod0 = model.matrix(~1, data=pd2)
sva1 = sva(as.matrix(out),mod,mod0,n.sv=2)
#head(sva1)
summary(lm(sva1$sv[,1] ~ pd2$Batch))
boxplot(sva1$sv[,2] ~ pd2$Batch, xlab="Batch", ylab="Surrogate Variable")
points(sva1$sv[,2] ~ jitter(as.numeric(pd2$Batch)),col=as.numeric(pd2$Batch))
modsv = cbind(mod,sva1$sv)
options(digits=3)
fit3 <- eBayes(combat_fit1)
writefile = topTable(fit3, n=Inf, genelist=genes, adjust.method = "fdr", sort.by="logFC")
writefile = topTable(fit3, number=Inf, genelist=genes, adjust.method = "BH", sort.by="logFC")



#read the phenotype file MODS_pheno.txt (formerly complex_course_pheno.txt) and the expression matrix out.txt from "/data/sayantan/comp_course_extension/Comp_course_extension"
arr=numeric(201)

MODS_pheno_subset=subset(MODS_pheno, MODS_pheno$Day.3.MODS.Trajectory_2_Day.7.MODS..1_No.Day.3.MODs.and.Yes.Day.7.MODS...0_No.Day.7.MODS!=1)
pheno=MODS_pheno_subset$Day.3.MODS.Trajectory_2_Day.7.MODS..1_No.Day.3.MODs.and.Yes.Day.7.MODS...0_No.Day.7.MODS
out=out[,MODS_pheno_subset$Microarray.ID]
#You can change the above variable to accomodate specific classes, for instance
#to exclude SIRS patients, complex_course_pheno_subset=subset(complex_course_pheno, complex_course_pheno$Classification!="SIRS")
pheno[pheno==0]="NO"
pheno[pheno==2]="YES"
arr[which(pheno=="YES")]="YES"
arr[which(pheno=="NO")]="NO"
#tail(clinical_data_without_controls)
#arr[which(arr=="0")]<-"X"
arr[which(arr=="YES")]<-"1"
arr[which(arr=="NO")]<-"0"
gsms=paste(arr, collapse = "")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sel <- which(sml != "X")
sml <- sml[sel]
sml
#x=out[,sel]
x=out
x[1:10,1:10]
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
g=as.data.frame(t(x))
g$description <- fl
arr=arr[which(arr=="0" | arr=="1")]
arr
design <- model.matrix(~0+factor(c(arr)))
colnames(design) <- levels(fl)
#x[] <- lapply(x, function(x) as.numeric(as.character(x)))
x1=as.data.frame(x)
fit <- lmFit(x1, design)
contrast.matrix <- makeContrasts(diff=G1-G0, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
genes=rownames(x1)
tT = topTable(fit2, number=Inf, genelist=genes, adjust.method = "BH", sort.by="logFC")
head(tT,15)
tT_pvalue=tT[order(tT$adj.P.Val,decreasing = FALSE),]
tT_pvalue_0.05=tT_pvalue[which(tT_pvalue$adj.P.Val<0.05),]
tT_signi_up_down=subset(tT_pvalue, tT_pvalue$logFC > 0.5 | tT_pvalue$logFC < -0.5) #18 genes
ordered_adj_pval=tT_signi_up_down[order(tT_signi_up_down$adj.P.Val, decreasing = FALSE),]
#plot heatmap
df=data.frame(l=arr)
colnames(df)="MODS Labels"
df[which(df$`MODS Labels`==0),"MODS Labels"]="Resolving MODS"
df[which(df$`MODS Labels`==1),"MODS Labels"]="Evolving MODS"

rownames(df)=colnames(x)
genes=c(rownames(tT_signi_up_down))
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
sub=x[c(which(rownames(x) %in% genes)),]
data_subset_norm <- t(apply(sub, 1, cal_z_score))
#plotting heatmap
pheatmap(data_subset_norm,annotation_col = df,annotation_row = t, cluster_cols = TRUE,show_colnames = FALSE,show_rownames = FALSE,clustering_method = "complete",cluster_rows = T,scale="row",color=colorRampPalette(c("green", "white", "red"))(50))
#plotting volcano plot
volcano_names <- subset(tT_pvalue_0.05, abs(logFC)>=0.5)
#volcano_names <- data.frame(ifelse(abs(writefile$logFC)>=1, writefile$ID, NA))
head(volcano_names)
class(volcano_names)

print("volcano_names")

volcanoplot(fit2, style = "p-value", highlight = 10, names = rownames(volcano_names),
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

with(tT_pvalue_0.05, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot", xlim=c(-2.5,4)))
# Add colored points: red if adj.P.Val<0.05, orange of log2FC>0.5, green if both)
with(subset(tT_pvalue_0.05, adj.P.Val<.05 ), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
with(subset(tT_pvalue_0.05, abs(logFC)>0.5), points(logFC, -log10(adj.P.Val), pch=20, col="orange"))
with(subset(tT_pvalue_0.05, adj.P.Val<.05 & abs(logFC)>0.5), points(logFC, -log10(adj.P.Val), pch=20, col="green"))
with(subset(tT_pvalue_0.05, adj.P.Val<.05 & logFC>0.5), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
with(subset(tT_pvalue_0.05, adj.P.Val<.05 & logFC < -0.5), points(logFC, -log10(adj.P.Val), pch=20, col="black"))

#Do you want to label the genes then use the following the code. Might be messy for logFC 0.5 cutoff - lots of genes
#library(calibrate)
with(subset(tT_pvalue_0.05, adj.P.Val<.05 & logFC>2), textxy(logFC, -log10(adj.P.Val), labs=ID, cex=.3))
with(subset(tT_pvalue_0.05, adj.P.Val<.05 & logFC < -1), textxy(logFC, -log10(adj.P.Val), labs=ID, cex=.3))

#with(subset(writefile, abs(logFC)>0.5), textxy(logFC, -log10(adj.P.Val), labs=ID, cex=.3))

legend(3,2.5, legend=c("Upregulated", "Downregulated","Not Significant"),bty = "n",
       col=c("blue", "black","red"), cex=0.75,pch=16)




#MA plot
number_genes <- dim(out)[1]

# False discovery rate cut off set to 0.05.

FDR_cutoff <- 0.05
p_values <- fit2$P.value
adjusted_p_values <- tT_pvalue_0.05$adj.P.Val

# Identify significant genes :

significant_genes <- tT_pvalue_0.05[tT_pvalue_0.05$adj.P.Val <=
                                      FDR_cutoff,]
gene_index <- rownames(significant_genes)

# MA-plot displaying the log fold change between diseased and healthy
#samples as a function of the average expression level across all
#samples.


status <- character (length=number_genes)
status <- rep ("not changing", number_genes)
names (status) <- seq (1,number_genes,1)
status [gene_index] <- "significant"

pdf(file="saving_plot8.1.pdf")
geneplotter::plotMA (tT_pvalue_0.05, status=status, col=c("red", "blue"))
text(x=11, y=3, labels=paste("FDR = ", FDR_cutoff), col="black", font=2)
abline(h=c(0.5,-0.5), col="green")
dev.off()



#GO analysis
Sign_genes_up = subset(tT_signi_up_down, logFC > 0.5)
Sign_genes_down = subset(tT_signi_up_down, logFC < -0.5)
de=tT_signi_up_down
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$logFC > 0.5 & de$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$logFC < -0.5 & de$adj.P.Val < 0.05] <- "DOWN"
p <- ggplot(data=de, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point() + theme_minimal()
# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

#GO analysis
all_genes_human = unique(rownames(out))
GO_up <- enrichGO(gene         = rownames(up),
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  universe      = all_genes_human,
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
GO_down <- enrichGO(gene         = rownames(down),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    universe      = all_genes_human,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
dotplot(GO_up, x="count", showCategory=30, color = 'p.adjust', title  = "DEG_UP COVID vs Normal")
dotplot(GO_down, x="count", showCategory=30, color = 'p.adjust', title  = "DEG_DOWN COVID vs Normal")
#plotting heatmap using pheatmap
df=data.frame(l=labels)



