## Analysis of Klingelhutz keratinocyte RNA-seq data
## Date: 03.25.19
## Author: Michael Chimenti
## Organism: hg38 / human 
## Aligners: hisat2 / salmon
## Design: Two cells lines +/- treatments with replicates
## Reps: 4

##########
## Imports
##########


#negative binomial GLM and related
library('DESeq2')
library('calibrate')
library('tximport')
library('readr')
#annotation
library('biomaRt')
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library('tidyverse')
library('pcaExplorer')
#pathway and gene clusters
library('DEGreport')
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

setwd("~/iihg/RNA_seq/klingelhutz/project_klingelhutz_rnaseq_032119/")

###########
##Function Defs
###########

get_res <- function(dds, meta_col, cond1, cond2, anno, pval = 0.05){
  res <- results(dds, contrast = c(meta_col,cond1,cond2))
  res <- na.omit(res)
  res_sig <- res[res$padj < pval & res$baseMean > 5.0,]
  res_ord <- res_sig[order(res_sig$padj),]
  res_ord$ext_gene <- anno[row.names(res_ord), "gene_name"]
  return(res_ord)
}

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#######################################
## tximport > DESeq2 > PCAExplorer
#######################################


samples <- read.table("samples.csv", sep=',', header=TRUE)
rownames(samples) <- samples$sample
#samples$batch <- as.factor(samples$batch)

files <- file.path(getwd(), samples$sample, 'salmon', 'quant.sf')
names(files) <- samples$sample

tx2gene <- read.csv(file.path(getwd(), "tx2gene.csv"), header = FALSE, as.is = c(1:2)) 

tx2gene$V1 <- tx2gene$V1 %>% 
  strsplit(split = '.', fixed = TRUE) %>%
  sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

samples$group <- paste(samples$line, samples$cond, sep='_')

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ group + rep)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi <- DESeq(ddsTxi)

##---------------launch PCA Explorer on dds object 
anno <- get_annotation(ddsTxi, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)
rldTxi <- rlog(ddsTxi, blind=FALSE)

#rld_tab <- as.data.frame(assay(rldTxi))
#rld_tab$gene <- anno[row.names(rld_tab), "gene_name"]
#write.csv(x=rld_tab, file = "rlog_transform_counts_genes.csv")

pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)
pcaExplorer::pcaplot(rldTxi, intgroup = "group", ellipse = FALSE, text_labels = TRUE)

plotDispEsts(ddsTxi)
plotMA(object = ddsTxi, alpha = 0.05, ylim = c(-10,10))

## CHecking CXCL8 expression b/c of qPCR data
cxcl8 <- plotCounts(ddsTxi, gene = "ENSG00000169429", intgroup="group", returnData = TRUE)
p <- ggplot(cxcl8, aes(x=as.factor(group), y = count)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.05)
p <- p + ggtitle("CXCL8 normalized expression in human keratinocytes")
p

#### SAMPLE SWAP
## Samples 1807TSS2 and 1807NT3 appear swapped by PCA / distance heatmap
## Samples 1818TSS1 and 1818NT2 appear swapped by PCA and distance heatmap

# samples_SWAPFIX <- read.table("samples_SWAPFIX.csv", sep=',', header=TRUE)
# rownames(samples_SWAPFIX) <- samples_SWAPFIX$sample
# #samples$batch <- as.factor(samples$batch)
# 
# files <- file.path(getwd(), samples_SWAPFIX$sample, 'salmon', 'quant.sf')
# names(files) <- samples_SWAPFIX$sample
# 
# tx2gene <- read.csv(file.path(getwd(), "tx2gene.csv"), header = FALSE, as.is = c(1:2)) 
# 
# tx2gene$V1 <- tx2gene$V1 %>% 
#   strsplit(split = '.', fixed = TRUE) %>%
#   sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists
# 
# txi <- tximport(files, type="salmon", tx2gene=tx2gene)
# 
# samples_SWAPFIX$group <- paste(samples_SWAPFIX$line, samples_SWAPFIX$cond, sep='_')
# 
# ddsTxi_FIX <- DESeqDataSetFromTximport(txi,
#                                    colData = samples_SWAPFIX,
#                                    design = ~ group + rep)
# 
# ddsTxi_FIX <- ddsTxi_FIX[ rowSums(counts(ddsTxi_FIX)) > 5, ]
# ddsTxi_FIX <- DESeq(ddsTxi_FIX)
# rldTxi_FIX <- rlog(ddsTxi_FIX, blind = FALSE)
# 
# pcaExplorer(dds = ddsTxi_FIX, rlt = rldTxi_FIX, annotation = anno)
# pcaExplorer::pcaplot(rldTxi_FIX, intgroup = "group", ellipse = FALSE, text_labels = TRUE)


##### DROP SAMPLES 1807TSS2, 1807NT3, 1818TSS1, 1818NT2
ddsTxi_drop <- ddsTxi[,!colData(ddsTxi)$sample %in% c("1807TSS2","1807NT3","1818TSS1","1818NT2")]
colData(ddsTxi_drop)

ddsTxi_drop <- DESeq(ddsTxi_drop)
rldTxi_drop <- rlog(ddsTxi_drop, blind=FALSE)

#rld_tab <- as.data.frame(assay(rldTxi))
#rld_tab$gene <- anno[row.names(rld_tab), "gene_name"]
#write.csv(x=rld_tab, file = "rlog_transform_counts_genes.csv")

pcaExplorer(dds=ddsTxi_drop,annotation=anno,rlt=rldTxi_drop)
pcaExplorer::pcaplot(rldTxi_drop, intgroup = "group", ellipse = FALSE, text_labels = TRUE)


##############
## DE analysis
##############

res_1807TSS_v_1807NT <- get_res(dds = ddsTxi_drop, meta_col = 'group', cond1 = 'euro_1807_TSS', cond2= 'euro_1807_NT', anno=anno)
res_1818TSS_v_1818NT <- get_res(dds = ddsTxi_drop, meta_col = 'group', cond1 = 'afro_1818_TSS', cond2 ='afro_1818_NT', anno=anno)
res_1818TSS_v_1807TSS <- get_res(dds = ddsTxi_drop, meta_col = 'group', cond1 = 'afro_1818_TSS', cond2 ='euro_1807_TSS', anno=anno)

res_1807SEB_v_1807NT <- get_res(dds=ddsTxi_drop, meta_col = 'group', cond1='euro_1807_SEB', cond2='euro_1807_NT', anno=anno)
res_1818SEB_v_1818NT <- get_res(dds=ddsTxi_drop, meta_col = 'group', cond1='afro_1818_SEB', cond2='afro_1818_NT', anno=anno)
res_1818SEB_v_1807SEB <- get_res(dds=ddsTxi_drop, meta_col = 'group', cond1='afro_1818_SEB', cond2='euro_1807_SEB', anno=anno)

res_1807SEB_v_1807TSS <- get_res(dds=ddsTxi_drop, meta_col = 'group', cond1='euro_1807_SEB',cond2='euro_1807_TSS',anno=anno)
res_1818SEB_v_1818TSS <- get_res(dds=ddsTxi_drop, meta_col = 'group', cond1='afro_1818_SEB',cond2='afro_1818_TSS',anno=anno)

## these analyses require a different model 

## All TSS vs All NT
## All SEB vs All NT
## All SEB vs All TSS
## 1807NT vs 1818NT

ddsTxi2 <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ line + cond + rep)  #assuming no interaction term...

ddsTxi2_drop <- ddsTxi2[,!colData(ddsTxi2)$sample %in% c("1807TSS2","1807NT3","1818TSS1","1818NT2")]

ddsTxi2_drop <- DESeq(ddsTxi2_drop)
res_TSS_v_NT <- get_res(dds = ddsTxi2_drop, meta_col = 'cond', cond1 = 'TSS', cond2 ='NT', anno=anno)
res_SEB_v_NT <- get_res(dds = ddsTxi2_drop, meta_col = 'cond', cond1 = 'SEB', cond2 = 'NT', anno=anno)
res_SEB_v_TSS <- get_res(dds=ddsTxi2_drop, meta_col = 'cond', cond1 = 'SEB', cond2 = 'TSS', anno=anno)
res_1807NT_v_1818NT <- get_res(dds=ddsTxi2_drop, meta_col = 'line', cond1 = 'euro_1807', cond2 = 'afro_1818', anno=anno)

## write gene lists
mycols <- c("baseMean","log2FoldChange","padj","ext_gene")

write.csv(res_1807TSS_v_1807NT[,mycols], file = "DE_genes_1807TSS_v_1807NT_padj_0p05.csv")
write.csv(res_1818TSS_v_1818NT[,mycols], file = "DE_genes_1818TSS_v_1818NT_padj_0p05.csv")
write.csv(res_1818TSS_v_1807TSS[,mycols], file = "DE_genes_1818TSS_v_1807TSS_padj_0p05.csv")

write.csv(res_1807SEB_v_1807NT[,mycols], file="DE_genes_1807SEB_v_1807NT_padj_0p05.csv")
write.csv(res_1818SEB_v_1818NT[,mycols], file="DE_genes_1818SEB_v_1818NT_padj_0p05.csv")
write.csv(res_1818SEB_v_1807SEB[,mycols], file="DE_genes_1818SEB_v_1807SEB_padj_0p05.csv")

write.csv(res_1807SEB_v_1807TSS[,mycols], file="DE_genes_1807SEB_v_1807TSS_padj_0p05.csv")
write.csv(res_1818SEB_v_1818TSS[,mycols], file="DE_genes_1818SEB_v_1818TSS_padj_0p05.csv")

write.csv(res_TSS_v_NT[,mycols], file="DE_genes_TSS_v_NT_padj_0p05.csv")
write.csv(res_SEB_v_NT[,mycols], file="DE_genes_SEB_v_NT_padj_0p05.csv")
write.csv(res_SEB_v_TSS[,mycols], file="DE_genes_SEB_v_TSS_padj_0p05.csv")
write.csv(res_1807NT_v_1818NT[,mycols], file="DE_genes_1807NT_v_1818NT_padj_0p05.csv")

########
png("tss_v_nt_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_TSS_v_NT, main = "Volcano Plot: DE genes TSS v NT", lfcthresh=1, sigthresh=0.05, textcx=.35, xlim=c(-5, 5), ylim = c(3,20))
dev.off()

png("seb_v_nt_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_SEB_v_NT, main = "Volcano Plot: DE genes SEB v NT", lfcthresh=1, sigthresh=0.05, textcx=.35, xlim=c(-5, 5), ylim = c(3,20))
dev.off()

png("seb_v_tss_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_SEB_v_TSS, main = "Volcano Plot: DE genes SEB v TSS", lfcthresh=1, sigthresh=0.05, textcx=.35, xlim=c(-5, 5), ylim = c(3,10))
dev.off()

png("1807TSS_v_1807NT_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_1807TSS_v_1807NT, main = "Volcano Plot: DE genes 1807TSS v. 1807NT", lfcthresh=1, sigthresh=0.05, textcx=.35, xlim=c(-5, 5), ylim = c(3,100))
dev.off()

png("1807SEB_v_1807NT_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_1807SEB_v_1807NT, main = "Volcano Plot: DE genes 1807SEB v. 1807NT", lfcthresh=1, sigthresh=0.05, textcx=.35, xlim=c(-5, 5), ylim = c(3,100))
dev.off()

## full exp datatable

rld_tab <- as.data.frame(assay(rldTxi_drop))
rld_tab$gene <- anno[row.names(rld_tab), "gene_name"]
write.csv(x=rld_tab, file = "rlog_transform_counts_genes.csv")

##################
## Heatmaps 

library("pheatmap")
library("RColorBrewer")
library("viridis")

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
################
## virus response 
virus_resp <- read.csv('virus_response.csv', header = FALSE, as.is = TRUE)
names(virus_resp) <- "gene_id"

virus_resp <- left_join(virus_resp, anno, on = 'gene_id')
virus_resp <- unique(virus_resp[,c("gene_id","gene_name")])


df <- as.data.frame(colData(ddsTxi_drop)[,"group"])


virus_tib <- as.data.frame(assay(rldTxi_drop)) %>%
  rownames_to_column(var = "gene_id") %>%
  as.tibble() %>%
  filter(gene_id %in% virus_resp$gene_id)

virus_tib <- left_join(virus_tib, virus_resp, by = "gene_id")
virus_mat <- as.matrix(virus_tib[,2:21])
row.names(virus_mat) <- virus_tib$gene_name
rownames(df) <- colnames(virus_mat)
colnames(df) <- "group"

#ordering <- c(2,10,15,5,8,13,3,7,11,4,6,12,1,9,14)
mat_breaks <- quantile_breaks(virus_mat, n = 11)

p <- pheatmap(virus_mat, 
              breaks = mat_breaks,
              cluster_rows=TRUE, 
              show_rownames=TRUE,
              cluster_cols=FALSE, 
              annotation_col=df,
              color = viridis(10),
              fontsize = 6,
              drop_levels = TRUE,
              show_colnames = FALSE,
              cutree_rows = 5,
              treeheight_row = 20,
              treeheight_col = 20,
              fontsize_row = 6,
              display_numbers = TRUE
)

png("virus_response_heatmap.png", 1700, 1100, pointsize=20, res=200)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

#######################
## virus membrane entry 
virus_resp <- read.csv('virus_membrane.csv', header = FALSE, as.is = TRUE)
names(virus_resp) <- "gene_id"

virus_resp <- left_join(virus_resp, anno, on = 'gene_id')
virus_resp <- unique(virus_resp[,c("gene_id","gene_name")])

library("pheatmap")
library("RColorBrewer")
library("viridis")

df <- as.data.frame(colData(ddsTxi_drop)[,"group"])


virus_tib <- as.data.frame(assay(rldTxi_drop)) %>%
  rownames_to_column(var = "gene_id") %>%
  as.tibble() %>%
  filter(gene_id %in% virus_resp$gene_id)

virus_tib <- left_join(virus_tib, virus_resp, by = "gene_id")
virus_mat <- as.matrix(virus_tib[,2:21])
row.names(virus_mat) <- virus_tib$gene_name
rownames(df) <- colnames(virus_mat)
colnames(df) <- "group"

#ordering <- c(2,10,15,5,8,13,3,7,11,4,6,12,1,9,14)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(virus_mat, n = 11)

p <- pheatmap(virus_mat, 
              breaks = mat_breaks,
              cluster_rows=TRUE, 
              show_rownames=TRUE,
              cluster_cols=FALSE, 
              annotation_col=df,
              color = viridis(10),
              fontsize = 6,
              drop_levels = TRUE,
              show_colnames = FALSE,
              cutree_rows = 5,
              treeheight_row = 20,
              treeheight_col = 20,
              fontsize_row = 6,
              display_numbers = TRUE
)

png("virus_membrane_response_heatmap.png", 1700, 1100, pointsize=20, res=200)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

###########
## GO annotation
library(clusterProfiler)
library(org.Hs.eg.db)

virus_resp_GO <- c("GO:0039502","GO:0039644","GO:0039657","GO:0039656","GO:0039652","GO:0051607","GO:0098586")
virus_memb_GO <- c("GO:0005534","GO:0005537","GO:0075509","GO:0016020","GO:0016021","GO:0001618","GO:0039707","GO:0046718")

TSS_NT_1807_genes <- bitr(row.names(res_1807TSS_v_1807NT), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]
TSS_NT_1807_ego<- enrichGO(gene = TSS_NT_1807_genes, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.csv(TSS_NT_1807_ego@result[,1:7], file = "GO_enrich_TSS_NT_1807.csv")

head(TSS_NT_1807_ego@result[,1:7], 40)

virus_resp_GO %in% TSS_NT_1807_ego@result$ID
virus_memb_GO %in% TSS_NT_1807_ego@result$ID

SEB_NT_1807_genes <- bitr(row.names(res_1807SEB_v_1807NT), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]
SEB_NT_1807_ego <- enrichGO(gene = SEB_NT_1807_genes, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
head(SEB_NT_1807_ego@result[,1:7], 40)
write.csv(SEB_NT_1807_ego@result[,1:7], file = "GO_enrich_SEB_NT_1807.csv")

TSS_NT_1818_genes <- bitr(row.names(res_1818TSS_v_1818NT), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]
TSS_NT_1818_ego <- enrichGO(gene = TSS_NT_1818_genes, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.csv(TSS_NT_1818_ego@result[,1:7], file = "GO_enrich_TSS_NT_1818.csv")

SEB_NT_1818_genes <- bitr(row.names(res_1818SEB_v_1818NT), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")[,2]
SEB_NT_1818_ego <- enrichGO(gene = SEB_NT_1818_genes, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
write.csv(SEB_NT_1818_ego@result[,1:7], file = "GO_enrich_SEB_NT_1818.csv")


