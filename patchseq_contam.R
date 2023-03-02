library(RColorBrewer)
library(egg)


Samp.datp <- readRDS("~/opc-ephys/Samp_datp_m.rds")
Expr.datp <- readRDS("~/opc-ephys/Expr_datp_m.rds")

kpSampP = 1:dim(Samp.datp)[1] # use all cells

annoPat_all = Samp.datp[kpSampP,]
annoPat_all$dendcluster_color = annoPat_all$cluster_color
datPat_all = as.matrix(Expr.datp[kpSampP,names(Expr.datp)!="sample_id"])
rownames(datPat_all) = annoPat_all$transcriptomics_sample_id
datPat_all = t(datPat_all)


# Format the data for FACS and patch-seq
rownames(datPat_all)[rownames(datPat_all) == '03-Mar'] = "MARCH4" # fix error
tmp = datPat_all
rownames(tmp) = make.names(rownames(tmp))
pat_df = as.data.frame(t(tmp[allMarkers, annoPat_all$transcriptomics_sample_id])+1)
pat_df$sample_id = rownames(pat_df)
rownames(datFACs2.m)[rownames(datFACs2.m) == '03-Mar'] = "MARCH4" # fix error
tmp = datFACs2.m
rownames(tmp) = make.names(rownames(tmp))
facs_df = as.data.frame(t(tmp[allMarkers,])+1)
facs_df$sample_id = rownames(facs_df)
facs_df$major_type = as.character(classBr.m.m)
facs_df$contam_type = as.character(clustersF.m)





# define which subclass each patch-seq cell is assigned to, based on maximal marker expression
nm  = names(markers)
isOn = substr(nm,nchar(nm)-2,nchar(nm))=="_on"
useThese = nm[isOn&(!is.element(nm,paste0(nm,"_on")))]
useThese = setdiff(useThese,c("CR_on","Meis2_on")) # These types are rare and unlikely to be actually patched.
subclassDat = calcContamAllTypes(pat_df, markers[useThese])  # Identify subclass based on marker gene expression
subclass = colnames(subclassDat)[subclassDat %>% apply(1,which.max)]
subclass = gsub("_on","",subclass)
pat_df$contam_type = subclass

tmp2 = match(pat_df$contam_type, annoFACs2.m$subclass_label_new) 
pat_df$major_type  = as.character(classBr.m)[tmp2]
pat_df$contam_type = paste0(pat_df$contam_type,"_on")
# check
tmp = annoPat_all %>% 
  dplyr::rename(sample_id = transcriptomics_sample_id) %>%
  dplyr::select(sample_id, corresponding_AIT2.3.1_alias)
pat_df = merge(pat_df, tmp, by = "sample_id", all.y = FALSE)





# calculate the QC metrics 
range01 = function(x){(x-min(x))/(max(x)-min(x))}
# calculate 
qcMetrics = calculatePatchSeqQCMetrics2(pat_df, facs_df, markers)
# NMS score of 0.4 as a pass/fail call
qcMetrics$QC_pass = c(TRUE,FALSE)[(qcMetrics$marker_sum_norm<0.40)+1]
qcMetrics = dplyr::rename(qcMetrics, Microglia = Macrophage)
qcMetrics$Microglia = range01(qcMetrics$Microglia)
dir.create("./output", showWarnings = FALSE)
#write.csv(qcMetrics, file = "./output/human_qcMetrics.csv")
#write.csv(qcMetrics, file = "./output/mouse_qcMetrics.csv")




num_markers = 20
plot_cell_types = c('Sst_on','Astro','OPC')
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
  filter(contam_type == "Sst_on") %>%
  #arrange(desc(Microglia)) %>%
  pull(sample_id)
# Expression matrix with only rows which are genes in the marker list and columns in samples included
exp_mat = datPat_all[unlist(trimmed_marker_list), pat_samps]
exp_mat = log2(exp_mat +1) 
#exp_mat = cbind(annoPat_all[annoPat_all$transcriptomics_sample_id %in% pat_samps,], exp_mat) 

order_samps = lapply(trimmed_marker_list, function(x){
  tmp = exp_mat[x,] #Keep only rows containing marker genes
  tmp = tmp[,order(colMeans(tmp), decreasing = T)] #Reorder sample columns based on colMeans
  colnames(tmp)
}) 

View(tmp[,order(colMeans(tmp), decreasing = T)])
x = trimmed_marker_list[[1]]

try = list(c(1,2,3), c(1,2,3), c(1,2,3))
lapply(try, function(x){print})



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
m1 = ggplot(df, aes(Var2, Var1, fill = value, color = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12)) +
  scale_color_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12), guide = "none") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  labs(x = "", y = "SST Interneuron\nmarkers", fill = "Log2\nExpression") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.margin = unit(c(0,0.5,0,0.5), "cm"),
        legend.position = "bottom",
        legend.justification = "left",
        legend.title = element_text(vjust = 1.3),
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
m2 = ggplot(df, aes(Var2, Var1, fill = value, color = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12)) +
  scale_color_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12), guide = "none") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  labs(x = "", y = "Astro\nmarkers", fill = "Log2\nExpression") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.5,0,0.5), "cm"),
        legend.position = "bottom",
        legend.justification = "left",
        legend.title = element_text(vjust = 1.3),
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
m3 = ggplot(df, aes(Var2, Var1, fill = value, color = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12)) +
  scale_color_gradientn(colors = pal, breaks = seq(0,12,2), labels = seq(0,12,2), limits = c(0,12), guide = "none") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
  labs(x = "", y = "OPC\nmarkers", fill = "Log2\nExpression") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.5,0,0.5), "cm"),
        legend.position = "bottom",
        legend.justification = "left",
        legend.title = element_text(vjust = 1.3),
        legend.text.align = 0,
        legend.title.align = 1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)
  )

mouse_heatmap_p = ggarrange(m1, NULL, m2, NULL, m3, ncol = 1, common.legend = T, legend = "bottom", align = "v", heights = c(1,-0.09,1,-0.09,1)) 
mouse_heatmap_p = annotate_figure(mouse_heatmap_p, top = text_grob("SST Interneurons (Mouse)", hjust = 0.3, size = 16))
mouse_heatmap_p











