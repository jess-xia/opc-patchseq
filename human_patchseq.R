setwd("C:/Tripathy Lab/opc_ephys")

library(data.table)
library(dplyr)
library(textshape)
library(tidyverse)


# read human patch-seq data (from https://github.com/keon-arbabi/patch-seq-microglia/blob/main/1-calculate-qc-metrics.Rmd)

Samp_datp_h = fread(
  file = "20200625_patchseq_metadata_human.csv",
  data.table = FALSE
)
Expr_datp_h = fread(
  file = "20200512_Human_PatchSeq_Release_cpm.csv",
  data.table = FALSE
) %>%
  dplyr::rename(gene = 1) %>%
  # Removes duplicate rows in a data frame, keeps all columns in input dataframe
  # There's only one duplicate, messed up because gene name are dates, randomly remove one
  dplyr::distinct(gene, .keep_all = TRUE) %>% 
  column_to_rownames(var = "gene") %>%
  # Swap row and columns so that each row is a sample and each column is a gene
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample_id")
# Order meta data to match order of expression data
Samp_datp_h = Samp_datp_h[match(Expr_datp_h$sample_id, Samp_datp_h$transcriptomics_sample_id), ]

saveRDS(Samp_datp_h, "Samp_datp_h.rds")

Samp_datp_h <- readRDS("~/opc-ephys/Samp_datp_h.rds")
Expr_datp_h <- readRDS("~/opc-ephys/Expr_datp_h.rds")

kpSampP = 1:dim(Samp_datp_h)[1] # use all cells

annoPat_all = Samp_datp_h[kpSampP, ]
annoPat_all$dendcluster_color = annoPat_all$cluster_color
datPat_all = as.matrix(Expr_datp_h[kpSampP, names(Expr_datp_h) != "sample_id"])
rownames(datPat_all) = annoPat_all$transcriptomics_sample_id
datPat_all = t(datPat_all)





# My code, trying things out
before_filter <- fread(file = "20200512_Human_PatchSeq_Release_cpm.csv",
                                   data.table = FALSE) %>%
  dplyr::rename(gene = 1) %>%
  filter(gene == "01-Mar") %>%
  column_to_rownames(var = "gene")

after_filter <- fread(file = "20200512_Human_PatchSeq_Release_cpm.csv",
                                  data.table = FALSE) %>%
  dplyr::rename(gene = 1) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  filter(gene == "01-Mar")


filter(Expr_datp_h, gene == "01-Mar")

subset(Expr_datp_h, select = c("PDGFRA"))


Expr_datp_h %>%
  ggplot(aes(x=sample_id, y=log2(MOG+1))) + 
  geom_point()




Expr_datp_h %>%
  ggplot(aes(x=sample_id, y=log2(PDGFRA+1))) + 
  geom_point()


#MOG, PDGFRA, C1QC










