## ----library, include = FALSE-----------------------------------------------------------------------------------------
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(tibble)
library(DT)

## ----import_quantification--------------------------------------------------------------------------------------------
stringTie_quant <- importIsoformExpression(
  parentDir = "./Source/bladder",
  addIsofomIdAsColumn = FALSE,
  readLength = 150
)

(head(stringTie_quant$abundance, 10))


## ----design_matrix-----------------------------------------------------------------------------------------------------
design <- data.frame(
  sampleID = colnames(stringTie_quant$abundance),
  condition = gsub(".{1}$", "", colnames(stringTie_quant$abundance)) # Remove number
)

(design)


## ----merge_file, warning = FALSE, message = FALSE---------------------------------------------------------------------
switch_analyze_Rlist <- importRdata(
  isoformCountMatrix   = stringTie_quant$counts,
  isoformRepExpression = stringTie_quant$abundance,
  designMatrix         = design,
  isoformExonAnnoation = "./Source/bladder/BCaMerge.gtf",
)

names(switch_analyze_Rlist)


## ----gene_count-------------------------------------------------------------------------------------------------------
gene_count <- extractGeneExpression(
  switch_analyze_Rlist,
  extractCounts = TRUE # set to FALSE for abundances
)

(head(gene_count,10))


## ----import_VGF-------------------------------------------------------------------------------------------------------
V38.gtf <- rtracklayer::import("./Source/gencode.v38.annotation.gtf")
unique(V38.gtf$gene_type)


## ----filter_lncRNA----------------------------------------------------------------------------------------------------
lncRNA_subset <- V38.gtf$gene_type == "lncRNA"
lncRNA <- V38.gtf[lncRNA_subset]$gene_name

gene_count_no_lncRNA <- gene_count %>% filter(!(gene_name %in% unique(lncRNA)))
(head(gene_count_no_lncRNA, 10))

## Get only count matrix
count_matrix <- gene_count_no_lncRNA %>%
  distinct(gene_name, .keep_all = TRUE) %>% # Remove duplicate gene_name 
  column_to_rownames("gene_name") %>% 
  select(-gene_id) %>% 
  as.matrix

(head(count_matrix, 10))


## ----boxplot_example--------------------------------------------------------------------------------------------------
count_matrix %>% 
  cpm(log=TRUE) %>% 
  head(20) %>% 
  t %>% 
  as.data.frame() %>% 
  rownames_to_column("type") %>% 
  tidyr::pivot_longer(-type) %>%
  mutate(type = gsub("\\d", "", type)) %>% 
  ggplot(aes(x = name, y =value, fill = type)) + geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  


## ----differential_expression------------------------------------------------------------------------------------------
library(edgeR)

colnames(count_matrix)
group <- design$condition
diff_design <- model.matrix(~0+group)
diff_design


## ----normalization----------------------------------------------------------------------------------------------------
dge <- DGEList(counts=count_matrix, group=group) 

keep <- filterByExpr(dge, group=group,min.count=2, min.prob=0.5)

dge <- dge[keep,]

# Calculate normalization factor
genexp <- calcNormFactors(dge)

# GLM Common dispersion
genexp <- estimateGLMCommonDisp(dge, diff_design)

# Estimate GLM trended dispersions
genexp <- estimateGLMTrendedDisp(genexp, diff_design)

# Tagwise dispersion of each gene
genexp <- estimateGLMTagwiseDisp(genexp, diff_design)

plotBCV(genexp)

(head(genexp$counts))


## ----fit_model--------------------------------------------------------------------------------------------------------
fit <- glmQLFit(genexp, diff_design)

genediff <- glmQLFTest(fit, contrast=c(-1,1))

# All genes
all_gene <- topTags(genediff, n = Inf, p.value = 1, adjust.method = "fdr")

all_gene$table %>% 
  rownames_to_column("gene_name") %>% 
readr::write_csv("all_gene.csv")
# Only significant value
sig_gene <- topTags(genediff, n = Inf, p.value = 0.05, 
                    adjust.method = "fdr", sort.by = "logFC")

# Total diffenetiated gene
summary(decideTests(genediff))

# Summary table
(sig_gene$table)


## ----pca_data---------------------------------------------------------------------------------------------------------
library(PCAtools)

# Calculate log-counts-per-million
logcpm <- cpm(dge, prior.count = 2, log = TRUE)

# Create a metadata table
metadata <- data.frame(row.names = colnames(logcpm),
                       group = c(rep(1,3), rep(2,3)))

(metadata)


## ----pca_analysis-----------------------------------------------------------------------------------------------------
# Perform PCA analysis
pc <- pca(logcpm, metadata = metadata, removeVar = 0.1)

# Create PCA plot
biplot(pc, colby = "group")


## ----pairsplot, fig.height=8------------------------------------------------------------------------------------------
pairsplot(pc)


## ----screeplot--------------------------------------------------------------------------------------------------------
screeplot(pc)


## ----heatmap, fig.height=7, fig.width=4-------------------------------------------------------------------------------
library(ComplexHeatmap)

# Filter only significant DEGs gene (from EdgeR)
DEGGene <- logcpm[rownames(sig_gene$table),]

Heatmap(DEGGene, row_km = 6, name = "RNAseq\nWorkshop",
        row_names_gp = gpar(fontsize = 3))  


## ----volcano_plot, fig.height=7---------------------------------------------------------------------------------------
library(EnhancedVolcano)

# Create Volcano plot
EnhancedVolcano(all_gene$table, lab = rownames(all_gene$table), x = "logFC", y = "PValue", 
                xlim = c(-12, 12), labSize = 2.0, pCutoff = 0.05,
                title = "RNAseq workshop", max.overlaps = 50)
