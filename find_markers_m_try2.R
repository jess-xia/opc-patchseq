suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(Seurat)
  library(edgeR)
  library(here)
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
})
'''

# read the data into R. *This step is slow*  
exons = read.csv(file = "/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_tasic_2018/mouse_VISp_2018-06-14_exon-matrix.csv", row.names = 1)
introns = read.csv(file = "/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_tasic_2018/mouse_VISp_2018-06-14_intron-matrix.csv", row.names = 1)
geneInfo = read.csv(file = "/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_tasic_2018/mouse_VISp_2018-06-14_genes-rows.csv", row.names = 1)
sampInfo = read.csv(file = "/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_tasic_2018/mouse_VISp_2018-06-14_samples-columns.csv", row.names = 1)

# convert the meta-data files into formats consistent with the rest of the analysis
sampInfo[is.na(sampInfo)] = 0
anno = auto_annotate(sampInfo) # Function in VENcelltypes
anno$sample_id = anno$sample_name

# convert the data into CPM(exons+introns) and format appropriately 
CPM = cpm(introns+exons)
rownames(CPM) = rownames(geneInfo)
colnames(CPM) = anno$sample_id

# format appropriately
data = as.data.frame(t(CPM))
data$sample_id = anno$sample_id

# write annotation file
write_feather(anno,"/external/rprshnas01/kcni/jxia/opc-ephys/anno.feather")
# write data file
write_feather(data,"/external/rprshnas01/kcni/jxia/opc-ephys/data.feather")
'''

Samp.dat = read_feather("/external/rprshnas01/kcni/jxia/opc-ephys/anno.feather") 
Expr.dat = feather("/external/rprshnas01/kcni/jxia/opc-ephys/data.feather") # FPKM

Samp.dat = Samp.dat[match(Expr.dat$sample_id,Samp.dat$sample_id),]

# Define a second annotation and data file with all clusters
ld = sort(unique(Samp.dat$cluster_label))
useClust2 = ld
for (val in c("ALM","Batch Grouping","Doublet","High Intron","Low Quality"))
  useClust2 = as.character(useClust2[!grepl(val,useClust2)]) #Keep in useClust2 if the values above are not in the cluster name
# THIS FUNCTION IS IN library(mfishtools)
kpSamp2 = subsampleCells(Samp.dat$subclass_label,100) #Subsets a categorical vector to include up to a maximum number of values for each category.
kpSamp2 = kpSamp2&is.element(Samp.dat$cluster_label,useClust2)



annoFACs2.m = Samp.dat[kpSamp2,]
datFACs2.m = as.matrix(Expr.dat[kpSamp2,names(Expr.dat)!="sample_id"])
rownames(datFACs2.m) = annoFACs2.m$sample_id
datFACs2.m = t(datFACs2.m)
annoFACs2.m$subclass_label = make.names(annoFACs2.m$subclass_label)
annoFACs2.m %<>% 
  relocate(subclass_label, .after = sample_name) %>%
  relocate(class_label, .after = subclass_label)




# A character vector of same length as names with each changed to a syntactically valid name
annoFACs2.m$subclass_label = make.names(annoFACs2.m$subclass_label)
# adjust subclass labels to match human MTG markers from Lee et al. (https://elifesciences.org/articles/65482#s2)
annoFACs2.m %<>%
  mutate(subclass_label_new = case_when(
    cluster_label %in% c("OPC Pdgfra Ccnb1", "OPC Pdgfra Grm5") ~ "OPC",
    TRUE ~ subclass_label #in all other cases, just set new column same as original
  )) %>%
  relocate(subclass_label_new, .after = subclass_label) #change column positions



# Define class labels
classBr.m = annoFACs2.m$subclass_label_new
classBr.m[annoFACs2.m$class_label!="Non-Neuronal"] = annoFACs2.m$class_label[annoFACs2.m$class_label!="Non-Neuronal"]
classBr.m = factor(classBr.m)
clustersF.m = factor(annoFACs2.m$subclass_label_new)
gc()


# find markers and format
markers = defineClassMarkers(datFACs2.m, clustersF.m, classBr.m, numMarkers = 50)
allMarkers = unique(unlist(markers))
markerTable = NULL
for (i in 1:length(markers))
  markerTable = cbind(markerTable, markers[[i]])
colnames(markerTable) = names(markers)
write.csv(markers,
          "/external/rprshnas01/kcni/jxia/opc-ephys/mouse_MTG_markers_calculated.csv",
          row.names = FALSE)
