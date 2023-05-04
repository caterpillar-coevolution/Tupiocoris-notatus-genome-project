#JKGs DESeq2 analysis of Tnot gene expression

library("pheatmap")
library("tidyverse")
library("DEGreport")
library( "DESeq2" )
library( "EnhancedVolcano" )
library(dplyr)
library(tibble)
library(ggplot2)

#######################################################################
# load the datasets for the first analysis (datura vs. tobacco)
metaData <- read.csv("Gene_Metadata.csv")
geneData <- read.csv("DEQ_dataset_022223.csv")

#set significance cutoff
padj.cutoff <- 0.05

# construct dataset from genecount table
dds <- DESeqDataSetFromMatrix(countData=geneData, 
                              colData=metaData, 
                              design=~plant, tidy = TRUE)
# define the baseline/control tratment
dds$plant <- relevel(dds$plant, ref = "tobacco")

# run analysis and save results
dds <- DESeq(dds)
res <- results(dds)

#order output by pvalue
res <- res[order(res$padj),]
head(res)

#make a subset of the data with only significant genes
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# print results as csvs
write.csv(as.data.frame(res),
          file="Tnot_datura_vs_tobacco_deseq_output.csv")
write.csv(as.data.frame(sig_res),
          file="Tnot_datura_vs_tobacco_deseq_output_sig_only.csv")

# make a PCA plot of the results by sample (colored according to host plant)
vsdata <- vst(dds, blind=FALSE)
PCAplot <- plotPCA(vsdata, intgroup="plant")
PCAplot + theme_classic()

# make a volcano plot of the results by individual genes
EnhancedVolcano(res,
                lab = NA,
                labSize = 2.0,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Datura wrightii fed vs. Nicotiana attenuata-fed',
                pointSize = 1.5,
                pCutoff = 0.01452,
                FCcutoff = 2,
                col=c('black', 'black', 'red3', 'red3'),
                colAlpha = 1.)
#########################################################################
# now do the second analysis
# Nicotiana attenuata fed individuals only
# essentially a re-analysis of Cristina's data

# load the datasets for the second analysis (Na-EV vs. Na-irAOC)
metaDataNAonly <- read.csv("Gene_Metadata_NAonly.csv")
geneDataNAonly <- read.csv("DEQ_dataset_022323_NAonly.csv")

# construct dataset from genecount table
ddsNAonly <- DESeqDataSetFromMatrix(countData=geneDataNAonly, 
                              colData=metaDataNAonly, 
                              design=~line, tidy = TRUE)
# define the baseline/control tratment
ddsNAonly$line <- relevel(ddsNAonly$line, ref = "NaEV")

# run analysis and save results
ddsNAonly <- DESeq(ddsNAonly)
resNAonly <- results(ddsNAonly)

# sort results by p-value (adjusted for multiple testing) and print them as a new csv
resNAonly <- resNAonly[order(resNAonly$padj),]
head(resNAonly)

#make a subset of the data with only significant genes
sig_resNAonly <- resNAonly %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

write.csv(as.data.frame(resNAonly),
          file="Tnot_NaOnly_deseq_output.csv")
write.csv(as.data.frame(sig_resNAonly),
          file="Tnot_NaOnly_deseq_output_sig_only.csv")


# make a PCA plot of the results by sample (colored according to host plant)
vsdataNAonly <- vst(ddsNAonly, blind=FALSE)
PCAplotNAonly <- plotPCA(vsdataNAonly, intgroup="line")
PCAplotNAonly + theme_classic()

# make a volcano plot of the results by individual genes
EnhancedVolcano(resNAonly,
                lab = NA,
                labSize = 2.0,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Na EV vs. NA irAOC',
                pointSize = 1.5,
                pCutoff = 0.00014,
                FCcutoff = 1,
                col=c('black', 'black', 'red3', 'red3'),
                colAlpha = 1.)



