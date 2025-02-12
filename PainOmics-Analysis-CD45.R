library(limma)
library(edgeR)
library(DESeq2)
library(snm)
library(statmod) 
library(biomaRt)
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
###For entrez IDS: entrezgene_id



####Load HTSeq files using DESeq2's function
####WILL NOT USE DESeq2 for analysis
directory <- ("C:/Users/varsh/Desktop/GT Research/PainOmics_2023/HTSeq_PainOmics")
sampleTable <- read.csv("C:/Users/varsh/Desktop/GT Research/PainOmics_2023/metadata_CD45.csv", header=TRUE)
#sampleTable <- read.csv("C:/Users/varsh/Desktop/GT Research/PainOmics_2023/metadata_CD71.csv", header=TRUE)
sampleTable <- na.omit(sampleTable)
sampleTable <- subset(sampleTable, sampleTable$Timepoint!="")
#sampleTable$sampleName <- paste0(sampleTable$CD71_NWGC_ID, ".txt")
sampleTable <- subset(sampleTable, Timepoint!="Baseline")
sampleTable <- subset(sampleTable, Timepoint!="Metformin_Baseline")
sampleTable$Timepoint <- as.factor(sampleTable$Timepoint)
sampleTable$Unique.ID <- as.factor(sampleTable$Unique.ID)
sampleTable$Chronic.Pain <- as.factor(sampleTable$Chronic.Pain)
sampleTable$Sex <- as.factor(sampleTable$Sex)
sampleTable$Age <- round(sampleTable$Age,0)
sampleTable$CD45_libprep_batch <- as.factor(sampleTable$CD45_libprep_batch)
#sampleTable$CD71_libprep_batch <- as.factor(sampleTable$CD71_libprep_batch)
sampleTable$HU_Therapy <- as.factor(sampleTable$HU_Therapy)
#sampleTable$sampleName <- sampleTable$CD45_NWGC_ID
#sampleTable$fileName <- paste0(sampleTable$CD45_NWGC_ID,".txt")
sampleTable$Status <- "Steady_State"
for(i in 1:nrow(sampleTable))
{
  if(sampleTable$Timepoint[i]=="Inpatient_VOC")
  {
    sampleTable$Status[i] <- "Inpatient_VOC"
  }
  else if(sampleTable$Timepoint[i]=="Inpatient_FU")
  {
    sampleTable$Status[i] <- "Inpatient_FU"
  }
}
sampleTable$Met_therapy <- "No"
for(i in 1:nrow(sampleTable))
{
  if(sampleTable$Timepoint[i]=="Metformin_MTD")
  {
    sampleTable$Met_therapy[i] <- "Yes"
  }
  else if(sampleTable$Timepoint[i]=="Inpatient_VOC" && sampleTable$Therapy[i]=="HU + Metformin")
  {
    sampleTable$Met_therapy[i] <- "Yes"
  }
  else if(sampleTable$Timepoint[i]=="Inpatient_FU" && sampleTable$Therapy[i]=="HU + Metformin")
  {
    sampleTable$Met_therapy[i] <- "Yes"
  }
}
sampleTable$AgeCategory <- "TenandBelow"
for(i in 1:nrow(sampleTable))
{
  if(sampleTable$Age[i]>10)
  {
    sampleTable$AgeCategory[i] <- "AboveTen"
  }
}
sampleTable$Met_therapy <- as.factor(sampleTable$Met_therapy)
sampleTable$AgeCategory <- as.factor(sampleTable$AgeCategory)
###To keep a subset of patients with chronic pain only
sampleTable <- subset(sampleTable, sampleTable$Chronic.Pain == "Yes")

#sampleTable_limma$Status <- as.factor(sampleTable_limma$Status)
sampleTable$Status <- as.factor(sampleTable$Status)
sampleTable$Status <- relevel(sampleTable$Status, ref="Steady_State")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                       directory=directory, design=~Status)
counts <- counts(ddsHTSeq)

####Filter out genes with zero expression in less than 10% of samples (34 for CD71 and 31 for CD45)
#counts_reduced <- rowSums(counts >= 1) >= 22 ##CD45
counts_reduced <- rowSums(counts >= 1) >= 4 ##CD45 chronic pain = 4, no CP =18
counts <- counts[counts_reduced, ]

####Build TMM-normalized matrix
dge <- DGEList(counts=counts)
#Calculate normalization factors
dge <- calcNormFactors(dge, method="TMM")
#Get counts normalized by TMM method - logCPM
#Do log2(cpm)
counts_TMM <- cpm(dge, log=T, prior.count=3)


counts_anno <- GeneSymbolAnnot(snm_counts)
tcell <- c("BIN1", "ANXA2R", "DDX24", "NOP53", "IMP3", "LAGE3", "LIME1", "OCIAD2",
           "SAE1", "SNRPD2")
subset_tcell <- counts_anno[counts_anno$hgnc_symbol %in% tcell,]
subset_tcell$ensembl_gene_id <- NULL
subset_tcell$hgnc_symbol <- NULL
trans <- t(subset_tcell)
results <- princomp(scale(trans))$scores
axis_tcell1 <- as.data.frame(results[,1])
colnames(axis_tcell1) <- c("pc1")
bcell <- c("AFF3", "BLK", "CD19", "CD72", "CD79A", "EBF1", "NIBAN3", "FCRLA",
           "POU2AF1", "VPREB3")
subset_bcell <- counts_anno[counts_anno$hgnc_symbol %in% bcell,]
subset_bcell$ensembl_gene_id <- NULL
subset_bcell$hgnc_symbol <- NULL
trans <- as.data.frame(t(subset_bcell))
set.seed(42)
results <- princomp(scale(trans))$scores
axis_bcell1 <- as.data.frame(results[,1])
colnames(axis_bcell1) <- c("pc1")


sampleTable$pc1_tcell <- axis_tcell1$pc1
sampleTable$pc1_bcell <- axis_bcell1$pc1
sampleTable$pc1_neut <- axis_neut$pc1
####SNM Analysis
##Sex and Libprep_batch as factors. 
##Not adjusting for age as longitudinal samples so the effect may be masked
bio.var <- model.matrix(~Status, data = sampleTable)
#adj.var <- model.matrix(~Sex + CD45_libprep_batch + pc1_tcell + pc1_bcell + AgeCategory + Met_therapy + Chronic.Pain, data = sampleTable)
#adj.var <- model.matrix(~Sex + CD45_libprep_batch + AgeCategory + Met_therapy + Chronic.Pain, data = sampleTable)
####Forchronic pain pc1_bcell + pc1_tcell
adj.var <- model.matrix(~Sex + pc1_bcell + pc1_tcell + CD45_libprep_batch, data = sampleTable)
####For no chronic pain
#adj.var <- model.matrix(~Sex + pc1_bcell + pc1_tcell + CD45_libprep_batch +Met_therapy + AgeCategory, data = sampleTable)
#adj.var <- model.matrix(~Sex + CD45_libprep_batch +Met_therapy + AgeCategory, data = sampleTable)
snm_obj <- snm(counts_TMM, bio.var, adj.var, rm.adj=T)
snm_counts <- snm_obj$norm.dat
snm_counts <- as.data.frame(snm_counts)
colnames(snm_counts) <- colnames(counts_TMM)


####Run Limma for DGE Analysis
design <- model.matrix(~0 + Status, sampleTable)
colnames(design) <- c("Steady_State", "Inpatient_FU", "Inpatient_VOC")
corfit <- duplicateCorrelation(snm_counts,design, block = sampleTable$Unique.ID)
fit <- lmFit(snm_counts, design, block = sampleTable$Unique.ID,
             correlation=corfit$consensus)
#colnames(design) <- c("Steady_State", "Inpatient_FU", "Inpatient_VOC")
contrast.matrix <- makeContrasts(Inpatient_FU-Steady_State, Inpatient_VOC-Steady_State, Inpatient_VOC-Inpatient_FU, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
#VOCvsteadyState
##1: FU vs SS 2: VOC vs SS, 3: VOC vs FU
res_top <- topTable(fit2, coef=2, genelist=fit$genes,number=Inf,adjust="BH", p.value =0.05)
#res_top <- topTable(fit, coef=ncol(design), genelist=fit$genes,number=Inf,adjust="BH", p.value =0.05)
res <- topTable(fit2, coef=2, genelist=fit$genes,number=Inf, adjust="BH")

####Annotate with gene symbols
#Remove version number
res <- GeneSymbolAnnot(res, justver = T)
#write.csv(res_x, "C:/Users/varsh/Desktop/GT Research/PainOmics_2023/CD71_DEGs_VOCvsSS_ChronicPain.csv", row.names=FALSE)

####Enrichment analysis
resSig <- res
#msigdbr_df <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
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
#pdf(file="C:/Users/varsh/Desktop/GT Research/PainOmics_2023/ASH_2023_Poster/Fig1.pdf", width = 12)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA - CD45+ cells (Limma matrix, sex, batch, age, CP and Met status with SNM, VOC vs SS)") + 
  theme_classic()  +
  scale_fill_manual(
    values = c("#f0e1f5","#333D79FF"))

fgseaResSubset <- subset(fgseaResTidy, fgseaResTidy$padj<0.05)
ggplot(fgseaResSubset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme(axis.text.x = element_text(size = 12)) +
  theme_classic(base_size = 12) +
  scale_fill_gradient(high = "#C33764", low = "#1D2761") +
  labs(fill="Adjusted p-val")
#dev.off()

ggplot(fgseaResSubset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme(axis.text.x = element_text(size = 12)) +
  theme_classic(base_size = 12) +
  scale_fill_gradient(high = "#EC2F4B", low = "#009FFF") +
  labs(fill="Adjusted p-val")

ggplot(fgseaResSubset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme(axis.text.x = element_text(size = 12)) +
  theme_classic(base_size = 12) +
  scale_fill_gradient(high = "#4B194D", low = "#273469") +
  labs(fill="Adjusted p-val")


snm_countsorig <- subset(snm_counts, select=sampleTable_noFU$sampleName)
snm_countsorig <- subset(snm_countsorig, rownames(snm_counts) %in% rownames(res_top))
snm_countsorig <- as.matrix(snm_countsorig)
Heatmap(snm_countsorig, cluster_columns=TRUE, cluster_rows=FALSE)
res_top2 <- res_top[order(-abs(res_top$logFC)),]
res_top2 <- res_top2[1:20,]
snm_countsorig2 <- subset(snm_counts, rownames(snm_counts) %in% rownames(res_top))
pheatmap(snm_countsorig, annotation_col=rowann, cluster_cols=FALSE)
rowann <- subset(sampleTable_noFU, select = c(sampleName, Status))
rowann <- rowann[order(rowann$Status),]
colnames(rowann) <- "Status"
snm_countsorig2 <- snm_countsorig2[,select=rownames(rowann)]
rowann$sampleName <- NULL

#sampleTablelong <- subset(sampleTable_noFU, sampleTable_noFU$Unique.ID %in% samplevec_ssvoc_CD45)
snm_countsorig <- subset(snm_counts, select=sampleTable_noFU$sampleName)
snm_countsorig <- subset(snm_countsorig, rownames(snm_counts) %in% rownames(res_top))
snm_countsorig <- as.matrix(snm_countsorig)


sampleTablelong <- subset(sampleTable_noFU, sampleTable_noFU$Unique.ID %in% samplevec_ssvoc_CD45)
snm_countsorig <- subset(snm_counts, select=sampleTablelong$sampleName)
snm_countsorig <- subset(snm_countsorig, rownames(snm_counts) %in% rownames(res_top))
snm_countsorig <- as.matrix(snm_countsorig)

snm_counts <- GeneSymbolAnnot(snm_counts)
genes <- c("MS4A4A", "SERPINB2", "IL1B", "FAM20A")
genes <- c("CD177", "SOCS3", "ANXA3")
subset_tcell <- snm_counts[snm_counts$hgnc_symbol %in% genes,]
subset_tcell$ensembl_gene_id <- NULL
subset_tcell$hgnc_symbol <- NULL
trans <- t(subset_tcell)
results <- princomp(scale(trans))$scores
axis_tcell1 <- as.data.frame(results[,1])
colnames(axis_tcell1) <- c("pc1")
summary(results) #Proportion of variance should be high
newdf <- axis_tcell1
newdf$Status <- sampleTable$Status
newdf$Status <- as.character(newdf$Status)
newdf <- subset(newdf, newdf$Status!="Inpatient_FU")
for(i in 1:nrow(newdf)){
  if( newdf$Status[i] == "Inpatient_VOC") { newdf$Status[i]="VOE"}}
ggplot(newdf, aes(x = Status, y = pc1, fill = Status)) +
  geom_boxplot() + geom_jitter() + ggtitle("PC1 Scores: Proposed VOE Biomarkers")+
  theme_minimal() + scale_fill_manual(values = c("#4B878BFF", "#D01C1FFF")) 


library(dplyr)

ss <- sampleTablehmm %>% distinct(Unique.ID, .keep_all=TRUE)


larger <- read.csv("C:/Users/varsh/Desktop/GT Research/PainOmics_2023/painomics_updated_metadata.csv", header=TRUE)
larger <- subset(larger, select=c("Date.of.collection_CD71", "CD71_NWGC_ID"))
larger <- subset(larger, larger$CD71_NWGC_ID %in% sampleTable$CD71_NWGC_ID)

crp <- c("PTX3", "FAM20A", "SERPINB2")
subset_plasma <- counts_anno[counts_anno$hgnc_symbol %in% crp,]
subset_plasma$ensembl_gene_id <- NULL
rownames(subset_plasma) <- subset_plasma$hgnc_symbol
subset_plasma$hgnc_symbol <- NULL
subset_plasma <- t(subset_plasma)
fam20a <- as.data.frame(subset_plasma[,1])
fam20a$gene <- "FAM20A"
colnames(fam20a) <- c("Expression", "Gene")
nr3c1 <- as.data.frame(subset_plasma[,2])
serpinb2 <- as.data.frame(subset_plasma[,3])
nr3c1$gene <- "CYP17"
serpinb2$gene <- "SERPINB2"
colnames(serpinb2) <- c("Expression", "Gene")
colnames(nr3c1) <- c("Expression", "Gene")

newdf <- rbind(fam20a,serpinb2)
newdf <- rbind(newdf, nr3c1)
newdf$Sample <- rownames(newdf)
ggplot(newdf, aes(x=Sample, y=Expression)) +
  geom_line(aes(color=factor(Gene), group=factor(Gene))) +
  geom_point(aes(color=factor(Gene))) + theme_classic()


##GSVA: Dot-Plot
counts_ssgsea <- snm_counts
#counts_ssgsea <- tpmmat - for running enrichment analysis with TPM values
#counts_ssgsea1 <- GeneSymbolAnnot(counts_ssgsea, justver = T)
counts_ssgsea <- as.matrix(counts_ssgsea1)
ssgsea_scores <- gsva(counts_ssgsea, gset.idx.list = pathwaysH, method="gsva")
sig_pathways <- fgseaResTidy$pathway[fgseaResTidy$padj < 0.05 &
                                       fgseaResTidy$NES > 1.5]
fgseaResTidy1 <- fgseaResTidy[order(fgseaResTidy$padj),]
sig_pathways <- fgseaResTidy1$pathway[fgseaResTidy1$padj < 0.05]
sig_pathways1 <- sig_pathways[1:6]
#3 columns: pathway name, sample name, enrichment score
# Load the tidyverse package

#Construct annotation dataframe
sampleTablepaired <- subset(sampleTable, sampleTable$Unique.ID %in% samplevec_ssvoc_CD45)
sampleTablepaired <- subset(sampleTablepaired, sampleTablepaired$Status!="Inpatient_FU")
STunique <- as.data.frame(sampleTable %>% distinct(Unique.ID,Status, .keep_all = TRUE))
sampleTablepaired <- subset(sampleTablepaired, sampleTablepaired$sampleName %in% STunique$sampleName)
#sampleTablepaired <- subset(sampleTablepaired, sampleTablepaired$Chronic.Pain=="No")
heatmap_ann <- sampleTablepaired
heatmap_ann$Chronic.Pain <- unfactor(heatmap_ann$Chronic.Pain)
heatmap_ann$Status <- unfactor(heatmap_ann$Status)
heatmap_ann <- subset(heatmap_ann, select=c(sampleName, Chronic.Pain, Status))
for(i in 1:nrow(heatmap_ann))
{
  if(heatmap_ann[i,]$Chronic.Pain=="No") 
  {
    heatmap_ann[i,]$Chronic.Pain="No Chronic Pain"
    }
 else {heatmap_ann[i,]$Chronic.Pain="Chronic Pain"}
}
#heatmap_ann <- subset(sampleTable, select=c(sampleName, condition, agecategory))
#Sort dataframe
ssgsea_scores <- subset(ssgsea_scores, rownames(ssgsea_scores) %in% sig_pathways1)
ssgsea_scores <- subset(ssgsea_scores, select = sampleTablepaired$fileName)
ssgsea_scores_t <- as.data.frame(t(ssgsea_scores))
ssgsea_scores_t$sample <- rownames(ssgsea_scores_t)
heatmap_ann <- heatmap_ann[order(heatmap_ann$Status, heatmap_ann$Chronic.Pain ),]
for (i in 1:nrow(heatmap_ann))
{
  if(heatmap_ann$Status[i]=="Steady_State")
  {heatmap_ann$Status[i]="Steady State"}
}
for (i in 1:nrow(heatmap_ann))
{
  if(heatmap_ann$Status[i]=="Inpatient_VOC")
  {heatmap_ann$Status[i]="VOE"}
}
heatmap_ann <- arrange(heatmap_ann, Status)
heatmap_ann <- arrange(heatmap_ann, desc(Chronic.Pain))
ssgsea_scores_t <- ssgsea_scores_t[match(heatmap_ann$sampleName, rownames(ssgsea_scores_t)),]
heatmap_ann$sampleName <- NULL
colnames(heatmap_ann) <- c("Pain Status", "Timepoint")

library(circlize)
col_fun = colorRamp2(c(0, 19.5), c("white", "#10b59c"))
ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))
ha = HeatmapAnnotation(df = heatmap_ann, 
                       col = list(Chronic.Pain=col_fun,
                                  Status = c("Steady State" = "#ADD8E6", "VOE" = "#00008B")))

#Format the SSGSEA scores matrix
ssgsea_scores_t$sample <- NULL
ssgsea_scores_t <-  t(ssgsea_scores_t)
ssgsea_scores_t <- ssgsea_scores_t[match(sig_pathways1, rownames(ssgsea_scores_t)),]
colnames(ssgsea_scores_t) <- gsub("_counts.txt", "", colnames(ssgsea_scores_t))
rownames(ssgsea_scores_t) <- gsub("^.{0,9}", "", rownames(ssgsea_scores_t))
rownames(ssgsea_scores_t) <- gsub("_", " ", rownames(ssgsea_scores_t))
totalpescore <- colSums(as.data.frame(ssgsea_scores_t))
z_scores <- (totalpescore-mean(totalpescore))/sd(totalpescore)
colnames(ssgsea_scores_t) <- round(z_scores,2) #For plotting only
Heatmap(ssgsea_scores_t, cluster_rows = F, cluster_columns = F,
        bottom_annotation=NULL,column_names_rot = 65, column_names_gp = gpar(fontsize = 14),
        top_annotation =  ha, heatmap_legend_param=list(title="GSVA Scores"))

sig_pathways1 <- sig_pathways[1:6]
ssgsea_scores
