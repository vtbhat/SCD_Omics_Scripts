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
sampleTable <- read.csv("C:/Users/varsh/Desktop/GT Research/PainOmics_2023/metadata_CD71.csv", header=TRUE)
#sampleTable <- read.csv("C:/Users/varsh/Desktop/GT Research/Papers_Proposals/Pain_Txp_Manuscript/metadataCD71_SRA.csv", header=TRUE)
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
sampleTable$CD71_libprep_batch <- as.factor(sampleTable$CD71_libprep_batch)
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
#sampleTable <- subset(sampleTable, sampleTable$Chronic.Pain == "No")

#sampleTable_limma$Status <- as.factor(sampleTable_limma$Status)
sampleTable$Status <- as.factor(sampleTable$Status)
sampleTable$Status <- relevel(sampleTable$Status, ref="Steady_State")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                       directory=directory, design=~Status)
counts <- counts(ddsHTSeq)

####Filter out genes with zero expression in less than 10% of samples (34 for CD71 and 31 for CD45)
counts_reduced <- rowSums(counts >= 1) >= 20 ##20 CD71 #16 for noCP
counts <- counts[counts_reduced, ]

####Build TMM-normalized matrix
dge <- DGEList(counts=counts)
#Calculate normalization factors
dge <- calcNormFactors(dge, method="TMM")
#Get counts normalized by TMM method - logCPM
#Do log2(cpm)
counts_TMM <- cpm(dge, log=T, prior.count=3)


##Sex and Libprep_batch as factors. 
##Not adjusting for age as longitudinal samples so the effect may be masked
bio.var <- model.matrix(~Status, data = sampleTable)
#adj.var <- model.matrix(~Sex + CD71_libprep_batch + AgeCategory + Met_therapy + Chronic.Pain, data = sampleTable)
####Forchronic apin
#adj.var <- model.matrix(~Sex  + CD71_libprep_batch, data = sampleTable)
####For no chronic pain
adj.var <- model.matrix(~Sex + CD71_libprep_batch +Met_therapy + AgeCategory, data = sampleTable)
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
fgseaResSubset <- subset(fgseaResTidy, fgseaResTidy$padj<0.05)
ggplot(fgseaResSubset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme(axis.text.x = element_text(size = 12)) +
  theme_classic(base_size = 12) +
  scale_fill_gradient(high = "#45ba87", low = "#023D54") +
  labs(fill="Adjusted p-val")
#dev.off()

##For acute pain
ggplot(fgseaResSubset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme(axis.text.x = element_text(size = 12)) +
  theme_classic(base_size = 12) +
  scale_fill_gradient(high = "#EC2F4B", low = "#009FFF") +
  labs(fill="Adjusted p-val")
