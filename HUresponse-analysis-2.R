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




####Run SNM for sex and age as fixed effects
bio.var <- model.matrix(~condition, data = sampleTable)
adj.var <- model.matrix(~sex + agecategory + ARC + ANC + WBC, data = sampleTable)
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




