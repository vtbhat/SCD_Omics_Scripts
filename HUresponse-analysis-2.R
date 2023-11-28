library(limma)
library(edgeR)
library(DESeq2)
library(snm)
library(statmod) 
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Homo.sapiens)
library(tibble)
library(stringr)
library(glmnet)
library(fgsea)
library(msigdbr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(caTools)
library(GSVA)
library(ComplexHeatmap)
library(GenomicFeatures)
library(PCAtools)
library(pvca)
library(ggforce)

####Function to remove Ensembl version ID and annotate gene symbols
GeneSymbolAnnot <- function(input_matrix, justver = F)
{ 
  ensemblgenes <- rownames(input_matrix)
  ensembl_nover <- c()
  for(i in 1:nrow(input_matrix))
  {
    ensembl_nover <- append(ensembl_nover, 
                            str_replace(ensemblgenes[i], "\\.\\d+$", ""))
  }
  rownames(input_matrix) <- ensembl_nover
  if(justver == F) {
  genes <- as.vector(rownames(input_matrix))
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  mappedgenes <- getBM(values=genes, attributes = c("ensembl_gene_id", "hgnc_symbol"), mart=mart)
  #Convert rownames to first column in count matrix
  input_matrix <- rownames_to_column(as.data.frame(input_matrix), "ensembl_gene_id")
  #Add HGNC gene symbols to count matrix
  input_matrix <- merge(x=input_matrix, y=mappedgenes, by="ensembl_gene_id", all.x=TRUE)
  }
  return(input_matrix)
}




####Load HTSeq files using DESeq2's function
####WILL NOT USE DESeq2 for analysis
directory <- ("C:/Users/varsh/Desktop/GT Research/HU_Response_Omics_2023/HTSeqCounts_nowashout")
sampleTable <- read.csv("C:/Users/varsh/Desktop/GT Research/HU_Response_Omics_2023/sampleTable_HU_nowashout.csv", header=TRUE)
#directory <- ("C:/Users/varsh/Desktop/GT Research/HU_Response_Omics_2023/HTSeqCounts_rnaseqenv")
#sampleTable <- read.csv("C:/Users/varsh/Desktop/GT Research/HU_Response_Omics_2023/sampleTable_HU.csv", header=TRUE)
sampleTable$condition <- as.factor(sampleTable$condition)
sampleTable$condition <- relevel(sampleTable$condition, ref="pre_HU")
sampleTable$sex <- as.factor(sampleTable$sex)
sampleTable$agecategory <- as.factor(sampleTable$agecategory)
#sampleTable$Drug.status <- as.factor(sampleTable$Drug.status)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                       directory=directory, design=~condition)
counts <- counts(ddsHTSeq)

####Filter out genes with zero expression in less than 10% of samples
counts_reduced <- rowSums(counts >= 1) >= 5
counts <- counts[counts_reduced, ]

####Build TMM-normalized matrix
dge <- DGEList(counts=counts)
#Calculate normalization factors
dge <- calcNormFactors(dge, method="TMM")
#Get counts normalized by TMM method - logCPM
#Do log2(cpm)
counts_TMM <- cpm(dge, log=T, prior.count = 3)



####PCA including washout samples
rownames(sampleTable) <- sampleTable$sampleName
rownames(sampleTable) <- gsub("_counts.txt", "",rownames(sampleTable))
colnames(snm_counts) <- gsub("_counts.txt", "",colnames(snm_counts))
p <- pca(snm_counts, metadata = sampleTable, removeVar = 0.1)
p$yvars <- gsub("_counts.txt", "", p$yvars)
png("Fig2.png", units="in", width=8, height=8, res=600)
biplot(p, lab=NULL, colby = c('Drug.status'), shape=('age.category'), 
       pointSize = 4,title = ('PCA biplot: Age + condition'),
       shapekey = c('Above5' = 15, 'Below5' = 17),
       xlim = c(-150,100), ylim = c(-80, 100),
       legendPosition="right")+ geom_point() + geom_line(aes(group=sampleTable$DNA_ID), linetype = 1) + theme_minimal() +
ggforce::geom_mark_ellipse(aes(fill = sampleTable$Drug.status, color = sampleTable$Drug.status)) +
  theme(panel.grid.minor = element_blank())
dev.off()
#xlim = c(-200,100), ylim = c(-100, 100),


 

####Principal Variance Components Analysis
sampleTable_pdata <- subset(sampleTable, select=c("agecategory", "sex", "condition"))
rownames(sampleTable_pdata) <- sampleTable$sampleName
#Construct expression set object
phenodata <- new("AnnotatedDataFrame", data = sampleTable_pdata)
pvca_input <- ExpressionSet(assayData = counts_TMM, phenoData = phenodata)
pct_threshold <- 0.6
batch.factors <- c("agecategory", "sex", "condition")
pvcaObj <- pvcaBatchAssess (pvca_input, batch.factors, pct_threshold)
#Plot PVCA results
bp <- barplot(pvcaObj$dat, xlab = "Effects",
            ylab = "Weighted average proportion variance",
            ylim= c(0,1.1),col = c("blue"),
            main="PVCA estimation bar chart", las=2)
axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, cex.names = 1.5, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)
#Plot PVCA results with better resoltuion
pvcadf <- t(as.data.frame(pvcaObj$dat*100))
pvcadf <- as.data.frame(pvcadf)
pvcadf$V1 <- pvcadf[order(pvcadf$V1),]
pvcadf$labels <- c("Sex:Condition", "Age:Sex", "Age:Condition", "Sex",
                   "Condition", "Age", "Residual")
bp <- barplot(pvcaObj$dat*100,  xlab= "Effects",names=pvcadf$labels,
        ylab = "Weighted average proportion variance",
        ylim= c(0,100),col = c("blue"),
        main="PVCA estimation bar chart", las=2) 
colnames(pvcadf) <- c("WAPV", 
                      "Effects")
pvcadf$WAPV <- round(pvcadf$WAPV, 2)
ggplot(data=pvcadf, aes(x=Effects, y=WAPV)) +
  geom_bar(stat="identity", fill="steelblue")+
  ylab("Weighted average proportion variance") + ylim(0, 100)+
  geom_text(aes(label=paste(WAPV, "%")), vjust=-0.3, size=3.5) + 
  theme_classic()+scale_x_discrete(limits = pvcadf$Effects)


####Sample-sample correlation matrix
cormat <- counts_TMM
cormat <- as.data.frame(cor(cormat))
cormat$sample <- rownames(cormat)
cormat <- cormat[match(sampleTable$sampleName, cormat$sample),]
Condition = sampleTable$condition
annot_sscor = HeatmapAnnotation(Condition = Condition)
cormat$sample <- NULL
colnames(cormat) <- gsub("_counts.txt", "", colnames(cormat))
rownames(cormat) <- gsub("_counts.txt", "", rownames(cormat))
cormat <- as.matrix(cormat)
Heatmap(cormat, cluster_rows = T, cluster_columns = T, 
        top_annotation =  annot_sscor, 
        heatmap_legend_param=list(title="Correlation Scores"),)

####Scatterplot for age effects of genes in nowashout matrix
agedf <- GeneSymbolAnnot(tpmmat)
agedf$ensembl_gene_id <- NULL
agedf <- t(agedf)
colnames(agedf) <- agedf[47,]
agedf <- agedf[-47,]
agedf <- as.data.frame(agedf)

agedf[] <- lapply(agedf, as.numeric)
agedf <- log2(agedf + 1)

agedf$agecategory <- sampleTable$agecategory
agedf$condition<- sampleTable$condition
agedf <- agedf[!duplicated(as.list(agedf))]
ggplot(agedf, aes(condition, ZBTB7B)) +
  geom_point(aes(color=agecategory), size=3)


####Run SNM for sex and age as fixed effects
bio.var <- model.matrix(~condition, data = sampleTable)
adj.var <- model.matrix(~sex + agecategory, data = sampleTable)
snm_obj <- snm(counts_TMM, bio.var, adj.var, rm.adj=T)
snm_counts <- snm_obj$norm.dat
snm_counts <- as.data.frame(snm_counts)
colnames(snm_counts) <- colnames(counts_TMM)

####Run Limma for DGE -- simple paired analysis
design <- model.matrix(~ condition, sampleTable)
corfit <- duplicateCorrelation(snm_counts,design, block = sampleTable$subject)
fit <- lmFit(snm_counts, design, block = sampleTable$subject,
             correlation=corfit$consensus)
fit <- eBayes(fit, trend=TRUE)
res_top <- topTable(fit, coef=ncol(design), genelist=fit$genes,number=Inf,adjust="BH", p.value =0.05)
res <- topTable(fit, coef=ncol(design), genelist=fit$genes,number=Inf, adjust="BH")

####Annotate with gene symbols
#Remove version number
res <- GeneSymbolAnnot(res, justver = T)

####Enrichment analysis
resSig <- res
msigdbr_df <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)
resSig$stat <- -log10(resSig$P.Value)*resSig$logFC
resSig$ensembl_gene_id <- rownames(resSig)
res2 <- resSig %>% 
  dplyr::select(ensembl_gene_id, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)
set.seed(42)
fgseaRes <- fgsea(pathways=pathwaysH, stats=ranks)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
#Plotting it
fgseaResTidy$pathway <- gsub("^.{0,9}", "", fgseaResTidy$pathway)
fgseaResTidy$pathway <- gsub("_", " ", fgseaResTidy$pathway)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA (Limma matrix, age and sex with SNM, subject with duplicateCorrelation())") + 
  scale_fill_manual(
    values = c("#d7e3f5","#89ABE3")) +theme_classic()

####Calculating TPMs
#Calculating transcript length from GTF file
#Importing the GTF file also used for HTSeq
txdb <- makeTxDbFromGFF("C:/Users/varsh/Desktop/GT Research/HU_Response_Omics_2023/gencode.v42.annotation.gtf",format="gtf")
#Collect exons by gene ID
exons.list.per.gene <- exonsBy(txdb,by="gene")
#For each gene, reduce all the exons to a set of non overlapping exons, 
#calculate their lengths (widths) and sum them
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
#Calculate TPMs
tpmmat <- counts / exonic.gene.sizes
tpmmat <- t( t(tpmmat) * 1e6 / colSums(tpmmat) )
tpmmat <- tpmmat[rowSums(tpmmat[])>0,]
#Duplicates: "ENSG00000230417" "ENSG00000276085"

##########################################################
####Correlation test between all genes and corrected.HbF
metadata_washout <- read.csv("C:/Users/varsh/Downloads/HUresponse_washout.csv", 
                             header=TRUE)
##To remove washout rows and row with negative deltaHbF
metadata_washout <- subset(metadata_washout, !(metadata_washout$RNA.ID %in% c("HU11", "HU12","HU21","HU22","HU35","HU36")))
#Number of correlation coefficients = number of genes
transpose_counts <- subset(counts_TMM, select = -c(HU11_counts.txt, HU12_counts.txt))
bio.var <- model.matrix(~corrected.HbF, data = metadata_washout)
adj.var <- model.matrix(~agecategory, data = metadata_washout)
snm_obj <- snm(transpose_counts, bio.var, adj.var, rm.adj=T)
snm_counts <- snm_obj$norm.dat
snm_counts <- as.data.frame(snm_counts)
colnames(snm_counts) <- colnames(transpose_counts)
transpose_counts <- as.data.frame(t(transpose_counts))
sampleTable$RNA.ID <- gsub("_counts.txt", "", sampleTable$sampleName)
snm_counts_HbF <- t(snm_counts)
library(ppcor)
corrgenes <- data.frame(matrix(NA, nrow = 23164, ncol = 3))
corrgenes[,1] <- colnames(snm_counts_HbF)
colnames(corrgenes) <- c("gene", "coeff", "pval")
for(i in 1:23164)
{
  cortestres <- cor.test(snm_counts_HbF[,i], as.vector(log(metadata_washout$corrected.HbF)), method="pearson")
  #cortestres <- pcor.test(transpose_counts[,i], as.vector(log(metadata_washout$corrected.HbF)), as.vector(metadata_washout$age))
  corrgenes[i,2] <- cortestres$estimate
  corrgenes[i,3] <- cortestres$p.value
}
corrgenes$padj <- p.adjust(corrgenes$pval, method = "BH")
#corrgenes$padj <- qvalue(corrgenes$pval)$qvalues
strongcor <- subset(corrgenes, corrgenes$padj < 0.05)
strongcor <- subset(strongcor, strongcor$coeff > 0.5 | strongcor$coeff < -0.5)
ggplot(metadata_washout, aes(y=log(corrected.HbF), x=age, color=Drug.status)) 
+geom_point() + geom_text(aes(label=RNA.ID))

#Annotate the strongcor matrix with the gene symbol function
rownames(strongcor) <- strongcor$gene
strongcor <- GeneSymbolAnnot(strongcor, justver=F)
#write.csv(strongcor, "C:/Users/varsh/Desktop/GT Research/HU_Response_Omics_2023/HU_and_log_HbF_SNMforAge_r6.csv", q=F, row.names=F)
genes_for_model <- rownames(strongcor)

####ReactomePA
library(clusterProfiler)
library(gson)
reactome <- read.gmt("C:/Users/varsh/Downloads/ReactomePathways.gmt")
snm_genes <- snm_counts
snm_genes <- GeneSymbolAnnot(snm_genes, justver=F)

snm_genes$symbol = mapIds(org.Hs.eg.db,
                          keys=snm_genes$ensembl_gene_id, #Column containing Ensembl gene ids
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
background <- as.vector(unique(snm_genes$hgnc_symbol))
background <- setdiff(background, c("CCL3L1"))
strongcor_neg <- subset(strongcor, strongcor$coeff<0)
strongcor_pos <- subset(strongcor, strongcor$coeff>0)
y <- as.vector(strongcor_neg$hgnc_symbol)
x <- enricher(gene=y, universe=background,
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = reactome)
res_enricher <- head(x)

####Overrepresentation analysis with the correlated genes
sign1 <- as.vector(strongcor$ensembl_gene_id)
msigdbr_df <- msigdbr(species = "human")
pathwaysH = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)
fisherTestRes <- as.data.frame(names(pathwaysH))
fisherTestRes$pvals <- NA
fisherTestRes$padj <- NA
fisherTestRes$genes <- NA
#Fisher's exact test
for(j in 1:length(pathwaysH))
{ 
  sign2 <- as.vector(unlist(pathwaysH[j]))
  contintable <- data.frame(matrix(NA, nrow = 2, ncol =2))
  contintable[1, 1] <- length(intersect(sign1, sign2))
  contintable[2,1] <- length(setdiff(sign1, sign2))
  contintable[1, 2] <- length(setdiff(sign2, sign1))
  contintable[2, 2] <- 22000 - (contintable[1,1] + contintable[2,1] + contintable[1,2])
  test <- fisher.test(contintable)
  fisherTestRes$pvals[j] <- test$p.value
  fisherTestRes$genes[j] <- paste(sign2, collapse =" ")
}

fisherTestRes$padj <- p.adjust(fisherTestRes$pvals, method = "BH")
hbfpathways <- subset(fisherTestRes, fisherTestRes$pvals < 0.05)
genes1 <- unlist(strsplit(gsub("\\.","",hbfpathways$genes[1])," "))
genes2 <- unlist(strsplit(gsub("\\.","",hbfpathways$genes[2])," "))
genes3 <- unlist(strsplit(gsub("\\.","",hbfpathways$genes[3])," "))
genes4 <- unlist(strsplit(gsub("\\.","",hbfpathways$genes[4])," "))
commongenes <- intersect(genes3, genes4) 
commongenes <- intersect(commongenes, genes1)
commongenes <- intersect(commongenes, genes2)


####Single-sample GSEA
counts_ssgsea <- counts_TMM
#counts_ssgsea <- tpmmat - for running enrichment analysis with TPM values
counts_ssgsea <- GeneSymbolAnnot(counts_ssgsea, justver = T)
ssgsea_scores <- gsva(counts_ssgsea, gset.idx.list = pathwaysH, method="gsva")
ordered <- c(seq(1,46,2), seq(2,47,2))
ssgsea_scores <- ssgsea_scores[, c(ordered)]
#Fetch significant pathways
#sig_pathways <- fgseaResTidy$pathway[fgseaResTidy$padj < 0.05 & 
#                       fgseaResTidy$NES > 1.5 | fgseaResTidy$NES < -1.5]
sig_pathways1 <- fgseaResTidy$pathway[fgseaResTidy$padj < 0.05 & 
                                       fgseaResTidy$NES > 1.5]
#sig_pathways1 <- paste0(sig_pathways1, " (UP)")
sig_pathways2 <- fgseaResTidy$pathway[fgseaResTidy$padj < 0.05 & 
                                        fgseaResTidy$NES < -1.5]
#sig_pathways2 <- paste0(sig_pathways2, " (DOWN)")
sig_pathways <- append(sig_pathways1, sig_pathways2)
ssgsea_scores <- subset(ssgsea_scores, rownames(ssgsea_scores) %in% sig_pathways)
#Construct annotation dataframe
heatmap_ann <- read.csv("C:/Users/varsh/Desktop/GT Research/HU_Response_Omics_2023/GBT_sampleTable_HU_noWashout.csv",
                        header= TRUE)
heatmap_ann <- subset(heatmap_ann, select=c(sampleName, condition, age))
#heatmap_ann <- subset(sampleTable, select=c(sampleName, condition, agecategory))
#Sort dataframe
ssgsea_scores_t <- as.data.frame(t(ssgsea_scores))
ssgsea_scores_t$sample <- rownames(ssgsea_scores_t)
heatmap_ann <- heatmap_ann[order(heatmap_ann$condition,heatmap_ann$age ),]
ssgsea_scores_t <- ssgsea_scores_t[match(heatmap_ann$sampleName, rownames(ssgsea_scores_t)),]
heatmap_ann$sampleName <- NULL
colnames(heatmap_ann) <- c("Condition", "Age")
library(circlize)
col_fun = colorRamp2(c(0, 19.5), c("white", "#10b59c"))
ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))
ha = HeatmapAnnotation(df = heatmap_ann, 
     col = list(Age=col_fun,
       Condition = c("pre_HU" = "#ADD8E6", "HU_MTD" = "#00008B")))
                                                    
ssgsea_scores_t$sample <- NULL
ssgsea_scores_t <-  t(ssgsea_scores_t)
ssgsea_scores_t <- ssgsea_scores_t[match(sig_pathways, rownames(ssgsea_scores_t)),]
colnames(ssgsea_scores_t) <- gsub("_counts.txt", "", colnames(ssgsea_scores_t))
rownames(ssgsea_scores_t) <- gsub("^.{0,9}", "", rownames(ssgsea_scores_t))
rownames(ssgsea_scores_t) <- gsub("_", " ", rownames(ssgsea_scores_t))
#colnames(ssgsea_scores_t) <- NULL #For plotting only
Heatmap(ssgsea_scores_t, cluster_rows = F, cluster_columns = F,
             bottom_annotation=NULL,
        top_annotation =  ha, heatmap_legend_param=list(title="GSVA Scores"))

#Test correlation between GSVA scores and condition age
gsva_main <- data.frame(matrix(nrow=0, ncol=4))
for(k in 6:9)
{
HU_age <- heatmap_ann$Age[24:46]
gsva <- ssgsea_scores_t[k, 24:46]
gsva_scat1 <- as.data.frame(cbind(HU_age, gsva))
gsva_scat1$condition = "pre-HU"
cf <- coef(lm(gsva_scat1$gsva~gsva_scat1$HU_age))
print(paste("Slope for pathway pre-HU", rownames(ssgsea_scores_t)[k],"is", cf[2]))
HU_age <- heatmap_ann$Age[1:23]
gsva <- ssgsea_scores_t[k, 1:23]
gsva_scat2 <- as.data.frame(cbind(HU_age, gsva))
gsva_scat2$condition = "HU MTD"
cf <- coef(lm(gsva_scat2$gsva~gsva_scat2$HU_age))
print(paste("Slope for pathway HU MTD", rownames(ssgsea_scores_t)[k],"is", cf[2]))
gsva_scat <- rbind(gsva_scat1, gsva_scat2)
gsva_scat$pathway <- rownames(ssgsea_scores_t)[k]
gsva_main <- rbind(gsva_main, gsva_scat)
}

ggplot(gsva_main, aes(x = HU_age, y = gsva)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)  +
  theme_light() + 
  ggtitle("GSVA scores vs Age - pre HU and HU MTD (Pathways downregulated at HU MTD)") +
  xlab("Age (yrs)") + 
  ylab("GSVA scores")+
  facet_grid(condition ~ pathway) +
  theme(strip.text = element_text(size = 11, 
                                    color="black", face="bold"))






####Axes of variation
tcell <- c("BIN1", "ANXA2R", "DDX24", "NOP53", "IMP3", "LAGE3", "LIME1", "OCIAD2",
           "SAE1", "SNRPD2")
reticulocyte <- c("EPB42", "GMPR", "IFIT1B", "OR2W3", "PBX1", "SELENBP1", "SLC4A1", "SLC6A10P", 
                  "SNCA", "TNS1")
neutrophil <- c("CXCR1", "C5AR1", "NUP214", "AQP9", "PHC2", "SIRPA", "TSEN34", 
                "MBOAT7", "HCK", "JAML")
bcell <- c("AFF3", "BLK", "CD19", "CD72", "CD79A", "EBF1", "NIBAN3", "FCRLA",
           "POU2AF1", "VPREB3")
interferon <- c("IFIT2", "HERC5", "RSAD2", "EPSTI1", "OAS3", "IRF7", "SAMD9L",
                "SERPING1", "MX1", "DDX58")
axisg <- c("BCLAF1", "DYRK1A", "HNRNPK", "NPTN", "TENT2", "SLK", "SRP54",
           "TRIM33", "WIPF1", "ZFAND5")
counts_anno <- GeneSymbolAnnot(counts_TMM)
subset_tcell <- counts_anno[counts_anno$hgnc_symbol %in% reticulocyte,]
subset_tcell$ensembl_gene_id <- NULL
subset_tcell$hgnc_symbol <- NULL
#Pre-HU dataframe
col_odd <- seq_len(ncol(subset_tcell)) %% 2
tcell_pre_HU <- subset_tcell[ , col_odd == 1]
results <- prcomp(t(tcell_pre_HU), scale=TRUE)$x
axis_tcell1 <- as.data.frame(results[,1])
colnames(axis_tcell1) <- c("pc1")
#HU MTD dataframe
tcell_HU_MTD <- subset_tcell[ , col_odd == 0]
results <- prcomp(t(tcell_HU_MTD), scale=TRUE)$x
axis_tcell2 <- as.data.frame(results[,1])
colnames(axis_tcell2) <- c("pc1")
summary(results) #Proportion of variance should be high
newdf <- rbind(axis_tcell1, axis_tcell2)
newdf$condition <- c(rep("pre_HU", 23), rep("HU_MTD", 23))
ggplot(newdf, aes(x = condition, y = pc1, color = condition)) +
  geom_point()
t.test(axis_tcell1$pc1, axis_tcell2$pc1, paired=F)
ggplot()

####Find genes in those 4 signatures
rownames(counts_TMM) <- str_replace(rownames(counts_TMM), "\\.\\d+$", "")
subset_comm <- counts_TMM[rownames(counts_TMM) %in% commongenes,]
subset_comm <- subset(subset_comm, select = -c(HU11_counts.txt, HU12_counts.txt))
results <- prcomp(t(subset_comm), scale=FALSE)$x
axis_subset_corr$comm <- results[,1]
ggplot(axis_subset_corr, aes(x = comm, y = PC1))  +
  geom_point()
cor.test(axis_subset_corr$PC1, axis_subset_corr$comm)

####Volcano plot
ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()
Heatmap(snm_counts)

####nonHbF Hu benefits
gmtfile <- gmtPathways("C:/Users/varsh/Downloads/gene_set_library_crisp.gmt")
resSig <- res
pathwaysH = gmtfile
resSig$stat <- -log10(resSig$P.Value)*resSig$logFC
resSig$ensembl_gene_id <- rownames(resSig)
res2 <- resSig %>% 
  dplyr::select(hgnc_symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(hgnc_symbol) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)
set.seed(42)
fgseaRes <- fgsea(pathways=pathwaysH, stats=ranks)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
#Plotting it
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA (Limma matrix, age and sex with SNM, subject with duplicateCorrelation())") + 
  theme_minimal()


####################################################
library(celldex)
ref <- BlueprintEncodeData()

rbc_genes <- gmtfile$'erythrocyte'
fetal_hemo <- c("HBA1")
rif <- read.csv("C:/Users/varsh/Downloads/rbc_rif_harmonizome.csv", header= F)
rif_genes <- rif$V3
cd22 <-read.csv("C:/Users/varsh/Desktop/GT Research/HU_Response_Omics_2023/PC1_ToppFun.csv", h =T)
cd22_path <- cd22$heme_scavenging
cd22_path <- setdiff(cd22_path, c("", "JCHAIN"))
####Calculate PC1 scores for RBC genes and see if it correlates with Hbf levels
counts_anno <- GeneSymbolAnnot(snm_counts)
#counts_anno <- subset(counts_anno, select  = -c(HU11_counts.txt, HU12_counts.txt))
subset_tcell <- counts_anno[counts_anno$hgnc_symbol %in% cd22_path,]
subset_tcell$ensembl_gene_id <- NULL
subset_tcell$hgnc_symbol <- NULL#Pre-HU dataframe
results <- prcomp(t(subset_tcell),  scale=T)$x
axis_tcell1 <- as.data.frame(results[,1])
colnames(axis_tcell1) <- c("pc1")
metadata_washout <- subset(metadata_washout, !(metadata_washout$RNA.ID %in% c("HU11","HU12","HU21","HU22","HU35","HU36")))
#metadata_washout <- read.csv("C:/Users/varsh/Downloads/HUresponse_washout.csv", 
#                             header=TRUE)
cor.test(axis_tcell1$pc1, log(metadata_washout$corrected.HbF))
ggplot(metadata_washout, aes(x=(corrected.HbF))) + geom_density()
ggplot(metadata_washout, aes(x=log(corrected.HbF))) + geom_density()


####Scatterplot of HBG1/HBG2 vs corrected.HbF
hbgplot <- subset(metadata_washout, select=c(corrected.HbF))
hbgplot$HBG1 <- as.numeric(counts_TMM['ENSG00000213934.9',])
hbgplot$HBG2 <- as.numeric(counts_TMM['ENSG00000196565.15',])
library(ggplot2)
ggplot(hbgplot, aes(x=corrected.HbF, y=HBG1, label=RNA.ID), size=15)  +
  geom_point(aes(color=Drug.status)) + xlab("Corrected HbF%") + ylab("HBG1 TMM/log")+ geom_text(hjust=0, vjust=0) 
hbg <- data.frame(exp=unlist(hbgplot, use.names = FALSE))
hbg$gene <- c(rep("HBG1", 50), rep("HBG2", 50))
ggplot(hbg, aes(x=gene, y=exp)) + geom_violin()+ geom_boxplot(width=0.1) + theme_classic()
#For HBG1/G2 split plot
hbgplot$corrected.HbF <- NULL
ordered <- c(seq(1,50,2), seq(2,51,2))
hbgplot <- hbgplot[c(ordered),]
hbg <- data.frame(exp=unlist(hbgplot, use.names = FALSE))
hbg$gene <- c(rep("HBG1: pre-HU", 25), rep("HBG1: HU MTD", 25),
              rep("HBG2: pre-HU", 25), rep("HBG2: HU MTD", 25))
hbg$condition <- c(rep("pre-HU", 25), rep("HU MTD", 25),
                   rep("pre-HU", 25), rep("HU MTD", 25))
hbg$reorder <- as.vector(1:nrow(hbg))
ggplot(hbg, aes(x=fct_reorder(gene, reorder, sum), y=exp, fill=condition)) + 
  geom_violin()+ geom_boxplot(width=0.1) + theme_classic() +
  xlab("Gene") + ylab("Expression (TMM-normalized and log-transformed)")+
  ggtitle("HBG1/G2 Expression in pre-HU and HU MTD")
 
####Linear model
transpose_counts1 <- subset(transpose_counts, select=as.vector(strongcor$gene))
corrgenes <- data.frame(matrix(NA, nrow = 23164, ncol = 4))
corrgenes[,1] <- colnames(transpose_counts)
colnames(corrgenes) <- c("gene", "coeff", "pval", "coeffage")
for(i in 1:23164)
{
  lmres <- lm(metadata_washout$corrected.HbF ~ metadata_washout$age + transpose_counts1[,i])
  lmres_est <- summary(lmres)
  corrgenes[i,2] <- lmres_est$coefficients[3,1]
  corrgenes[i,3] <- lmres_est$coefficients[3,4]
  corrgenes[i,4] <- lmres_est$coefficients[2,1]
}
corrgenes$padj <- p.adjust(corrgenes$pval, method = "BH")
x <- lm(metadata_washout$corrected.HbF ~ metadata_washout$agebinary + transpose_counts[,2])
lmres <- summary(x)


####Correlation between ANC and log fold change
corrgenes <- data.frame(matrix(NA, nrow = 23164, ncol = 3))
corrgenes[,1] <- colnames(res_topx)
colnames(corrgenes) <- c("gene", "coeff", "pval")
for(i in 1:23164)
{
  cortestres <- cor.test(transpose_counts[,i], as.vector(log(metadata_washout$deltaANC)))
  corrgenes[i,2] <- cortestres$estimate
  corrgenes[i,3] <- cortestres$p.value
}
corrgenes$padj <- p.adjust(corrgenes$pval, method = "BH")

###ANC correlation calculation
corrgenes <- data.frame(matrix(NA, nrow = 50, ncol = 3))
corrgenes[,1] <- colnames(ssgsea_scores)
colnames(corrgenes) <- c("pathway", "coeff", "pval")
metadata_washout2 <- subset(metadata_washout, Drug.status=="HU MTD")
metadata_washout2 <- subset(metadata_washout2, !(metadata_washout2$RNA.ID %in% c("HU21","HU22","HU35","HU36")))
for(i in 1:50)
{
  cortestres <- cor.test(ssgsea_scores[,i], as.vector(metadata_washout2$ANC))
  corrgenes[i,2] <- cortestres$estimate
  corrgenes[i,3] <- cortestres$p.value
}
corrgenes$padj <- p.adjust(corrgenes$pval, method = "BH")
####Cell type deconvolution with dtangle
library(celldex)
ref <- BlueprintEncodeData(ensembl = T)
celltyperef <- assays(ref)
celltyperef <- t(celltyperef[[1]])
counts_TMM <- GeneSymbolAnnot(counts_TMM, justver =T)
counts_TMM_t <- t(counts_TMM)
newvec <- as.vector(colnames(celltyperef))
ints <- intersect(newvec, colnames(counts_TMM_t))
counts_TMM_t1 <- subset(counts_TMM_t, select = ints)
celltyperef <- subset(celltyperef, select = ints)
counts_TMM_t1[counts_TMM_t1 < 0] <- 0  
counts_TMM_t1 <- as.data.frame(counts_TMM_t1)
celltyperef <- as.data.frame(celltyperef)
counts_TMM_t1 <- counts_TMM_t1[ints]
celltyperef <- celltyperef[ints]
celltyperef <- as.matrix(celltyperef)
counts_TMM_t1 <- as.matrix(counts_TMM_t1)
library(dtangle)
dtn <- dtangle(counts_TMM_t1,references=celltyperef, data_type = "rna-seq")
res1 <- dtn$estimates