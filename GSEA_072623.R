#jay trying GSEA for Tnotatus data - Gene enriched when feeding on Datura, relative to Tobacco

#install the package
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("DBI")

#load the package
require("clusterProfiler")
require("tidyverse")

#load the data. 3x pre-sorted csv files
upreg <- read.csv("genelist_upreg.csv")
baseline <- read.csv("genelist_baseline_expression.csv")

#load the TERM2GENE file (connects gene IDs to their corresponding Go terms)
bp_terms <- read.csv("gene2bp.csv")
mf_terms <- read.csv("gene2mf.csv")
cc_terms <- read.csv("gene2cc.csv")

#load the GO2NAME files
bp_names <- read.csv("bp2name.csv")
mf_names <- read.csv("mf2name.csv")
cc_names <- read.csv("cc2name.csv")


#subset upregulated genes to only include genes with functional annotations
bp_joined <- inner_join(upreg, bp_terms, by = "gene_ID") 
mf_joined <- inner_join(upreg, mf_terms, by = "gene_ID") 
cc_joined <- inner_join(upreg, cc_terms, by = "gene_ID") 

#subset baseline genelist like above
bp_joined_base <- inner_join(baseline, bp_terms, by = "gene_ID") 
mf_joined_base <- inner_join(baseline, mf_terms, by = "gene_ID") 
cc_joined_base <- inner_join(baseline, cc_terms, by = "gene_ID") 

#write joined datasets to csvs for later reference
write.csv(as.data.frame(bp_joined),
          file="bp_joined.csv")
write.csv(as.data.frame(mf_joined),
          file="mf_joined.csv")
write.csv(as.data.frame(cc_joined),
          file="cc_joined.csv")

#do the same for baseline datasets
write.csv(as.data.frame(bp_joined_base),
          file="bp_joined_baseline.csv")
write.csv(as.data.frame(mf_joined_base),
          file="mf_joined_baseline.csv")
write.csv(as.data.frame(cc_joined_base),
          file="cc_joined_baseline.csv")

#QC BP dataset
 bp_genes = bp_joined[,2]
 names(bp_genes) = as.character(bp_joined[,1])
 bp_genes = sort(bp_genes, decreasing = TRUE)

# run analysis and make dotplot (BP, upreg)
y <- GSEA(bp_genes, TERM2GENE = bp_terms, TERM2NAME = bp_names)
write.csv(y,
          file="bp_results.csv")
dotplot(y, showCategory=30)

#QC MF dataset
mf_genes = mf_joined[,2]
names(mf_genes) = as.character(mf_joined[,1])
mf_genes = sort(mf_genes, decreasing = TRUE)

# run analysis and make dotplot (BP, upreg)
x <- GSEA(mf_genes, TERM2GENE = mf_terms, TERM2NAME = mf_names)
write.csv(x,
          file="mf_results.csv")
dotplot(x, showCategory=30)

#QC MF dataset
cc_genes = cc_joined[,2]
names(cc_genes) = as.character(cc_joined[,1])
cc_genes = sort(cc_genes, decreasing = TRUE)

# run analysis and make dotplot (BP, upreg)
z <- GSEA(cc_genes, TERM2GENE = cc_terms, TERM2NAME = cc_names)
write.csv(z,
          file="cc_results.csv")
dotplot(z, showCategory=30)

#QC BP baseline dataset
bp_based = bp_joined_base[,2]
names(bp_based) = as.character(bp_joined_base[,1])
bp_based = sort(bp_based, decreasing = TRUE)

# run analysis and make dotplot (BP, baseline)
a <- GSEA(bp_based, TERM2GENE = bp_terms, TERM2NAME = bp_names)
head(a)
dotplot(a, showCategory=30)

#QC MF baseline dataset
mf_based = mf_joined_base[,2]
names(mf_based) = as.character(mf_joined_base[,1])
mf_based = sort(mf_based, decreasing = TRUE)

# run analysis and make dotplot (BP, baseline)
b <- GSEA(mf_based, TERM2GENE = mf_terms, TERM2NAME = mf_names)
head(b)
dotplot(b, showCategory=30)

#QC CC baseline dataset
cc_based = cc_joined_base[,2]
names(cc_based) = as.character(cc_joined_base[,1])
cc_based = sort(cc_based, decreasing = TRUE)

# run analysis and make dotplot (BP, baseline)
c <- GSEA(cc_based, TERM2GENE = cc_terms, TERM2NAME = cc_names)
head(c)
dotplot(c, showCategory=30)


