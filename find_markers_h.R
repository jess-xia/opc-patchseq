devtools::install_github("AllenInstitute/patchseqtools")
devtools::install_github('PavlidisLab/patchSeqQC')
devtools::install_github("AllenInstitute/VENcelltypes")

library(dplyr)
library(data.table)
library(textshape)
library(tidyverse)

library(magrittr)
library(Seurat)
library(edgeR)
#library(here)
library(data.table)
library(ggpubr)
library(patchseqtools)
library(patchSeqQC)
library(VENcelltypes)
library(feather)
library(matrixStats)
library(ggplotify)
library(cowplot)
library(ggpubr)
library(grid)
library(svglite)
library(RColorBrewer)
library(ggbeeswarm)


# load and format human FACs data
Samp.dat.h = fread(file = "/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/metadata.csv", data.table = F) %>% 
  filter(region_label == "MTG") %>% #middle temporal gyrus, where the patchseq data was obtained
  filter(!(subclass_label == "" | is.na(subclass_label))) %>% #only keep cells with celltype labels
  relocate(subclass_label, .after = sample_name) #reorganize column order
rownames(Samp.dat.h) = Samp.dat.h$sample_name #use sample names as rownames

table(Samp.dat.h$subclass_label) #counts for each subclass
leaveout = Samp.dat.h %>% 
  dplyr::count(subclass_label, sort = TRUE) %>% 
  filter(n < 20) %>% 
  pull(subclass_label) #identify the subclasses that have less than 20 counts
kpSamp2 = Samp.dat.h %>% 
  filter(!(subclass_label %in% leaveout)) %>% 
  pull(sample_name) #identify names of samples that are in subclasses with greater than 20 counts



Expr.dat = fread(file = "/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/matrix.csv", data.table = F) %>%
  column_to_rownames(var = "sample_name") #count matrix, same as metadata, the sample names are set as the row names

annoFACs2.h = Samp.dat.h[kpSamp2,] #subset metadata of cells that belong to subclasses with greater than 20 counts
datFACs2.h = Expr.dat[match(annoFACs2.h$sample_name, rownames(Expr.dat)),] #Reorder count matrix such that the rows are in the same order as that in the metadata
datFACs2.h = edgeR::cpm(t(datFACs2.h)) #computes cpm, cpm typically have decimal points?



# A character vector of same length as names with each changed to a syntactically valid name
annoFACs2.h$subclass_label = make.names(annoFACs2.h$subclass_label)
# adjust subclass labels to match human MTG markers from Lee et al. (https://elifesciences.org/articles/65482#s2)
annoFACs2.h %<>%
  mutate(subclass_label_new = case_when(
    subclass_label %in% c("IT", 'L4.IT') ~ "Superficial.Layers",
    subclass_label %in% c("L5.6.IT.Car3",  "L5.6.NP", "L5.ET",  "L6.CT", "L6b") ~ "Deep.Layers",
    subclass_label %in% c("LAMP5",  "PAX6") ~ "LAMP5.PAX6.Other",
    #subclass_label == "Oligodendrocyte" ~ "Oligo.OPC",
    subclass_label == "OPC" ~ "Oligo.OPC",
    subclass_label == "Microglia" ~ "Microglia",
    TRUE ~ subclass_label
  )) %>%
  relocate(subclass_label_new, .after = subclass_label) #change column positions

# Define class labels
classBr.h = annoFACs2.h$subclass_label_new
# Based on original subclass labels, rename the neuronal subtypes to gabaergic or glutamatergic
classBr.h[annoFACs2.h$class_label!="Non-neuronal"] = annoFACs2.h$class_label[annoFACs2.h$class_label!="Non-neuronal"]
classBr.h = factor(classBr.h)
clustersF.h = factor(annoFACs2.h$subclass_label_new)

gc()

tmp = cbind(annoFACs2.h, datFACs2.h['PVALB',])
names(tmp)[ncol(tmp)] = "gene"
tmp %>% ggplot(aes(x = gene, y = subclass_label_new)) + geom_jitter()



# find markers and format
# markers = defineClassMarkers(datFACs2.h, clustersF.h, classBr.h, numMarkers = 50)
# allMarkers = unique(unlist(markers))
# markerTable = NULL
# for (i in 1:length(markers))
#   markerTable = cbind(markerTable, markers[[i]])
# colnames(markerTable) = names(markers)
# write.csv(markers,
#           "/external/rprshnas01/kcni/jxia/opc-ephys/human_MTG_markers_calculated.csv",
#           row.names = FALSE)

markers = read_csv("opc-ephys/human_MTG_markers_calculated.csv")
markers = as.list(markers)
allMarkers = unique(unlist(markers))


# Format the data for FACS and patch-seq
rownames(datPat_all)[rownames(datPat_all) == '03-Mar'] = "MARCH4" # fix error
tmp = datPat_all
rownames(tmp) = make.names(rownames(tmp))
pat_df = as.data.frame(t(tmp[allMarkers, annoPat_all$transcriptomics_sample_id])+1)
pat_df$sample_id = rownames(pat_df)
rownames(datFACs2.h)[rownames(datFACs2.h) == '03-Mar'] = "MARCH4" # fix error
tmp = datFACs2.h
rownames(tmp) = make.names(rownames(tmp))
facs_df = as.data.frame(t(tmp[allMarkers,])+1)
facs_df$sample_id = rownames(facs_df)
facs_df$major_type = as.character(classBr.h)
facs_df$contam_type = as.character(clustersF.h)


# define which subclass each patch-seq cell is assigned to, based on maximal marker expression
nm  = names(markers)
isOn = substr(nm, nchar(nm) - 2, nchar(nm)) == "_on"
useThese = nm[isOn & (!is.element(nm, paste0(nm, "_on")))]
useThese = setdiff(useThese, c("CR_on", "Meis2_on")) # These types are rare and unlikely to be actually patched.
subclassDat = calcContamAllTypes(pat_df, markers[useThese])  # Identify subclass based on marker gene expression
subclass = colnames(subclassDat)[subclassDat %>% apply(1, which.max)]
subclass = gsub("_on", "", subclass)
pat_df$contam_type = subclass

tmp2 = match(pat_df$contam_type, annoFACs2.h$subclass_label_new)

pat_df$major_type  = as.character(classBr.h)[tmp2]
pat_df$contam_type = paste0(pat_df$contam_type, "_on")
# check
tmp = annoPat_all %>%
  dplyr::rename(sample_id = transcriptomics_sample_id) %>%
  dplyr::select(sample_id, corresponding_AIT2.3.1_alias)
pat_df = merge(pat_df, tmp, by = "sample_id", all.y = FALSE)


# calculate the QC metrics 

range01 = function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# calculate
qcMetrics = calculatePatchSeqQCMetrics2(pat_df, facs_df, markers)
# NMS score of 0.4 as a pass/fail call
qcMetrics$QC_pass = c(TRUE, FALSE)[(qcMetrics$marker_sum_norm < 0.40) +
                                     1]

qcMetrics$Microglia = range01(qcMetrics$Microglia)

write.csv(qcMetrics, file = "/external/rprshnas01/kcni/jxia/opc-ephys/qcMetrics.csv")


# plot marker expression heatmap
## human

if(spcs == "human"){
  num_markers = 20
  plot_cell_types = c('Astrocyte','Microglia')
  plot_marker_list = c(markers[plot_cell_types])
  # don't include unknown marker genes in plot
  rm_markers = str_detect(allMarkers, "^LOC|^LINC|^KIAA|^SLC|^MT-|^RP[0-9]|^BC[0-9]|-PS")
  rm_markers = allMarkers[rm_markers]
  
  trimmed_marker_list = lapply(plot_marker_list, function(x){
    ind = which(x %in% rownames(datPat_all))
    tmp = make.names(x[ind])
    tmp = tmp[which(!tmp %in% rm_markers)]
    tmp = head(tmp, num_markers)
  }
  )
  pat_samps = qcMetrics %>%
    filter(contam_type == "Superficial.Layers_on") %>%
    #arrange(desc(Microglia)) %>%
    pull(sample_id)
  exp_mat = datPat_all[unlist(trimmed_marker_list), pat_samps]
  exp_mat = log2(exp_mat +1) 
  #exp_mat = cbind(annoPat_all[annoPat_all$transcriptomics_sample_id %in% pat_samps,], exp_mat) 
  
  order_samps = lapply(trimmed_marker_list, function(x){
    tmp = exp_mat[x,]
    tmp = tmp[,order(colMeans(tmp), decreasing = T)]
    colnames(tmp)
  }) 
  order_markers = lapply(trimmed_marker_list, function(x){
    tmp = exp_mat[x,]
    tmp = tmp[order(rowMeans(tmp), decreasing = F),]
    rownames(tmp)
  })
  
  pal = colorRampPalette(rev(brewer.pal(n = 21, name = "RdYlBu")))(20)
  
  m = datPat_all[trimmed_marker_list[[1]], pat_samps]
  m = log2(m +1) 
  df = melt(m) %>%
    mutate(Var1 = factor(Var1, levels = order_markers[[1]])) %>%
    mutate(Var2 = factor(Var2, levels = order_samps[[2]]))
  h1 = ggplot(df, aes(Var2, Var1, fill = value, color = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12)) +
    scale_color_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12), guide = "none") +
    guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
    labs(x = "", y = "Astrocyte\nmarkers", fill = "Log2\nExpression") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = unit(c(0,0.5,0,0.5), "cm"),
          legend.position = "top",
          legend.justification = "center",
          legend.text.align = 0,
          legend.title.align = 1,
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16)
    )
  
  m = datPat_all[trimmed_marker_list[[2]], pat_samps]
  m = log2(m +1) 
  df = melt(m) %>%
    mutate(Var1 = factor(Var1, levels = order_markers[[2]])) %>%
    mutate(Var2 = factor(Var2, levels = order_samps[[2]]))
  h2 = ggplot(df, aes(Var2, Var1, fill = value, color = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12)) +
    scale_color_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12), guide = "none") +
    guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
    labs(x = "", y = "Microglia\nmarkers", fill = "Log2\nExpression") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,0.5,0,0.5), "cm"),
          legend.position = "top",
          legend.justification = "center",
          legend.text.align = 0,
          legend.title.align = 1,
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16)
    )
  
  m = datPat_all[trimmed_marker_list[[3]], pat_samps]
  m = log2(m +1) 
  df = melt(m) %>%
    mutate(Var1 = factor(Var1, levels = order_markers[[3]])) %>%
    mutate(Var2 = factor(Var2, levels = order_samps[[2]]))
  h3 = ggplot(df, aes(Var2, Var1, fill = value, color = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12)) +
    scale_color_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12), guide = "none") +
    guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
    labs(x = "", y = "Oligodendrocyte\nmarkers", fill = "Log2\nExpression") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,0.5,0,0.5), "cm"),
          legend.position = "top",
          legend.justification = "center",
          legend.text.align = 0,
          legend.title.align = 1,
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16)
    )
  
  human_heatmap_p = ggarrange(h1, NULL, h2, NULL, NULL, ncol = 1, common.legend = T, legend = "none", align = "v", heights = c(1,-0.09,1,-0.09,1)) 
  human_heatmap_p = annotate_figure(human_heatmap_p, top = text_grob("Superficial Layer Neurons (Human)", hjust = 0.3, size = 16))
  human_heatmap_p 
}
