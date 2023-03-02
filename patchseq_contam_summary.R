library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

#write.csv(qcMetrics, file = "~/opc-ephys/output/mouse_qcMetrics.csv")
mouse_qc = read.csv(file = "./output/mouse_qcMetrics.csv") %>% mutate(Group = "Mouse")

joined <- left_join(mouse_qc, Samp_datp_m, by = c("sample_id" = "transcriptomics_sample_id"))

# By cluster

# Mean contamination by cluster
# GABAergic, Endothelial, Glutamatergic not included due to many -inf scores
data_long <- gather(joined, contam_glia, contam_score, c(Astro, Microglia:VLMC), factor_key=TRUE)
contam_sum <- data_long %>% 
  group_by(corresponding_AIT2.3.1_alias, contam_glia) %>%
  summarise(contam_score = mean(contam_score))
# Round contam scores to 2 decimal places for display on heatmap
contam_sum$contam_score <- round(contam_sum$contam_score ,digit=2)
# Remove clusters with no names
contam_sum <- contam_sum[!contam_sum$corresponding_AIT2.3.1_alias == '',]
# Add new column with subclass names taken from the first word of the cluster
contam_sum$subfromcluster <- sapply(contam_sum$corresponding_AIT2.3.1_alias, word)




# Order the glia based on average contamination
score_by_glia <- data_long %>% 
  group_by(contam_glia) %>%
  summarise(contam_score = mean(contam_score))
score_by_glia <- arrange(score_by_glia, desc(contam_score))
glia_order <- score_by_glia %>% 
  pull(contam_glia)


# Heatmap
# Plot average contamination score for each glia for each neural cell type
ggplot(contam_sum,                                
       aes(x=corresponding_AIT2.3.1_alias, y=factor(contam_glia, level=glia_order), fill = contam_score)) +
  geom_tile() +
  #geom_text(aes(label = contam_score), color = "white", angle = 90) +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab('contam_glia') +
  xlab('cluster')

# OPC only
contam_opc <- data_long %>%
  filter(contam_glia == "OPC") 
contam_opc <- contam_opc[!contam_opc$corresponding_AIT2.3.1_alias == '',] %>%
  subset(select = c(corresponding_AIT2.3.1_alias, contam_glia, contam_score))
contam_opc$subfromcluster <- sapply(contam_opc$corresponding_AIT2.3.1_alias, word)





# Boxplot for OPC (ordered)
ggplot(contam_opc, aes(x=reorder(corresponding_AIT2.3.1_alias, -contam_score, median), y=contam_score)) + 
  geom_boxplot(fill = "#4271AE", colour = "#1F3552", # Colors
               alpha = 0.9, outlier.colour = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab('opc_contam_score') +
  xlab('cluster')

# Boxplot for OPC (unordered)
ggplot(contam_opc, aes(x=corresponding_AIT2.3.1_alias, y=contam_score)) + 
  geom_boxplot(fill = "#4271AE", colour = "#1F3552", # Colors
               alpha = 0.9, outlier.colour = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab('opc_contam_score') +
  xlab('cluster')





# By subclass
contam_sum_sub <- data_long %>% 
  group_by(contam_type, contam_glia) %>%
  summarise(contam_score = mean(contam_score))
# Round contam scores to 2 decimal places for display on heatmap
contam_sum_sub$contam_score <- round(contam_sum_sub$contam_score ,digit=2)
contam_sum_sub <- contam_sum_sub[!contam_sum_sub$contam_type == '',]

# By subclass
contam_sum_sub_fromcluster <- contam_sum %>% 
  group_by(contam_glia, subfromcluster) %>%
  summarise(contam_score = mean(contam_score))
# Round contam scores to 2 decimal places for display on heatmap
contam_sum_sub_fromcluster$contam_score <- round(contam_sum_sub_fromcluster$contam_score ,digit=2)
contam_sum_sub_fromcluster <- contam_sum_sub_fromcluster[!contam_sum_sub_fromcluster$subfromcluster == '',]


# Heatmap
ggplot(contam_sum_sub,                                
       aes(x=contam_type, y=factor(contam_glia, level=glia_order), fill = contam_score)) +
  geom_tile() +
  geom_text(aes(label = contam_score), color = "white") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab('contam_glia') +
  xlab('subclass')

# Heatmap with subclass names taken from cluster names
ggplot(contam_sum_sub_fromcluster,                                
       aes(x=reorder(subfromcluster, -contam_score, median), y=factor(contam_glia, level=glia_order), fill = contam_score)) +
  geom_tile() +
  geom_text(aes(label = contam_score), color = "white") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab('contam_glia') +
  xlab('subclass_from_cluster_name')




# OPC only
contam_opc_sub <- data_long %>%
  filter(contam_glia == "OPC") 
contam_opc_sub <- contam_opc_sub[!contam_opc_sub$contam_type == '',] %>%
  subset(select = c(contam_type, contam_glia, contam_score)) 

# Boxplot for OPC (ordered)
ggplot(contam_opc_sub, aes(x=reorder(contam_type, -contam_score, median), y=contam_score)) + 
  geom_boxplot(fill = "#4271AE", colour = "#1F3552", # Colors
               alpha = 0.9, outlier.colour = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab('opc_contam_score') +
  xlab('subclass')


# Boxplot for OPC (ordered) - based on cluster substring
ggplot(contam_opc, aes(x=reorder(subfromcluster, -contam_score, median), y=contam_score)) + 
  geom_boxplot(fill = "#4271AE", colour = "#1F3552", # Colors
               alpha = 0.9, outlier.colour = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab('opc_contam_score') +
  xlab('subclass_from_cluster_name')





