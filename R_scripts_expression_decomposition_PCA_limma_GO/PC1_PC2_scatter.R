library(wesanderson)
library(ggrepel)
library(data.table)
library(dplyr)
library(svglite)
library(ggplot2)       # For plotting
library(ellipse)       # For ellipse calculations
library(grDevices) 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ellipse",version= '3.19')


colors=c('iPSC-derived neurons'='#469433',
'Fetal neurons'='#9c2725',
'Post-mortem neurons'='#034e91',
'NPC'='#e38519',
'iPSC'='#696969')

#colors=c('iPSC-derived neurons'='brown',
#         'Fetal neurons'='brown',
#         'Post-mortem neurons'='brown',
#         'NPC'='brown',
#         'iPSC'='brown')

shapes <- c('iPSC-derived neurons'=16,
            'Fetal neurons'=16,
            'Post-mortem neurons'=16,
            'NPC'=16,
            'iPSC'=16)

sizes <- c('iPSC-derived neurons'=3.5,
              'Fetal neurons'=3.5,
              'Post-mortem neurons'=3.5,
              'NPC'=3.5,
              'iPSC'=3.5)

shapes_sample <- c('Eres'=16,'Artrial'=17,'Ballarino'=17,'Rajarajan'=16,'Zaghi'=17,'Our data'=16,'Rahman'=16,'Hu'=15,'Tian'=17,'Heffel(adult)'=17,'Heffel(infant)'=17,'Heffel(2T)'=17,
                   'Heffel(3T)'=17,'Rahman'=16,'Wu'=16,'Li'=16,'Rajarajan'=16,'Ballerino'=17,'Zaghi'=17,'Rahman'=17,'Our_data'=16)
### 'o','*','*','o','*','o','o','s','*','*','*','*','*','o','o','o','o','*','*','*','o'
plot_pca_final_groups_human_withNames_CustomEllipse <- function(d1, plot_name, pc1_var, pc2_var, xaxis, yaxis) {
  d1$group <- factor(d1$group)
  num_groups <- length(unique(d1$group))
  
  if (num_groups > length(colors)) {
    stop("Not enough colors provided for the number of groups")
  }
  
  # Calculate ellipses manually for each group with at least 2 points
  ellipse_df <- d1 %>%
    group_by(group) %>%
    do({
      if (nrow(.) >= 2) {               
        cov_mat <- cov(.[, c(xaxis, yaxis)])
        center  <- colMeans(.[, c(xaxis, yaxis)])
        ellipse_points <- ellipse::ellipse(cov_mat, centre = center, level = 0.95)
        colnames(ellipse_points) <- c(xaxis, yaxis)
        data.frame(ellipse_points, group = .$group[1])
      } else {
        data.frame()
      }
    }) %>%
    ungroup()
  
  p1 <- ggplot(data = d1, aes_string(x = xaxis, y = yaxis, color = "group")) +
    geom_point(aes(shape = group, size = group)) +
    scale_shape_manual(values = shapes) +
    scale_size_manual(values = sizes) +
    xlab(paste0(xaxis, " ", pc1_var, "%")) +
    ylab(paste0(yaxis, " ", pc2_var, "%")) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12, family = 'Helvetica'),
          axis.title = element_text(size = 14, family = 'Helvetica')) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkgrey') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'darkgrey')
  
  if (nrow(ellipse_df) > 0) {
    p1 <- p1 + geom_polygon(data = ellipse_df, aes_string(x = xaxis, y = yaxis, fill = "group"), alpha = 0.1, show.legend = FALSE)
  }
  
  p_with_labels <- p1 + geom_text_repel(aes(label = sample_simple), force = 10,
                                        point.padding = 0.35, show.legend = TRUE)
  
  full_path <- sprintf("/media/anna/Samsung_T5/culture_vs_postmortem/%s", plot_name)
  print(full_path)
  ggsave(full_path, plot = p_with_labels, width = 8.1, height = 6)
  
  return(p_with_labels)
}

df_human <- fread('/media/anna/Samsung_T5/culture_vs_postmortem/PC1_PC2_data_for_autism.txt', sep='\t')
df_human
pca_human <- plot_pca_final_groups_human_withNames_CustomEllipse(df_human, "snp_autism_PCA_data_R.pdf", 10.4 , 7.8,  "PC1", "PC2")
pca_human

df_human <- fread('/media/anna/Samsung_T5/culture_vs_postmortem/PC1_PC2_data_for_fires.txt', sep='\t')
df_human
pca_human <- plot_pca_final_groups_human_withNames_CustomEllipse(df_human, "PC1_PC2_fires_restr_R.pdf", 22 , 16,  "PC1", "PC2")
pca_human



df_human <- fread('/media/anna/Samsung_T5/culture_vs_postmortem/PC1_PC2_data_for_BD_new_from_tab28_0.21_0.11_aggregated_check.txt', sep='\t')
df_human
pca_human <- plot_pca_final_groups_human_withNames_CustomEllipse(df_human, "snp_BD_tab_28_check_R.pdf",21,11,  "PC1", "PC2")
pca_human

df_human <- fread('/media/anna/Samsung_T5/culture_vs_postmortem/PC1_PC2_data_for_MDD_0.22_0.13_aggregated_check.txt', sep='\t')
df_human
pca_human <- plot_pca_final_groups_human_withNames_CustomEllipse(df_human, "snp_MDD_check_R.pdf",19.8,14.2,  "PC1", "PC2")
pca_human