## Load the necessary Packages

library(tidyverse)
library(microeco)
library(dplyr)
library(janitor)
library(agricolae)
library(ggplot2)
library(ggdendro)
library(magrittr)
library(ggtree)
library(BiocManager)

BiocManager::install('ggtree')

aa## Import and Read the Dataset (OTU and Sample data) 

otu_1 <- read.csv('Microbiome Analysis/otu_data.csv', 
                  row.names = 1, check.names = FALSE)

sample_1 <- read.csv('Microbiome Analysis/sample_data.csv',
                     row.names = 1, check.names = FALSE)

taxa_1 <- read.csv('Microbiome Analysis/taxa_data.csv',
                   row.names = 1, check.names = FALSE)

### View  and Inspect the imported datasets
## View the taxonomy datasets
taxa_1 |> 
  glimpse()

taxa_1 |> 
  view()

## View the metadata datasets
sample_1 |> 
  glimpse()

sample_1 |> 
  view()

## View the OTU/ASV datasets
otu_1 |> 
  glimpse()

otu_1 |> 
  view()

## Create a microtable Object
micro_table <- microtable$new(
  otu_1,
  sample_table = sample_1,
  tax_table = taxa_1,
  phylo_tree = NULL,
  rep_fasta = NULL,
  auto_tidy = TRUE
)

micro_table |> 
  str()

## To make the OTU and Sample Information consistent across all files
## in the dataset object
## Tidy dataset

taxa_1 |>
  tidy_taxonomy()

micro_table$tax_table |> 
  base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")

micro_table$filter_pollution(taxa = c('c__mitochondria',
                                      'c__chloroplast'))

micro_table$tidy_dataset()
print(micro_table)

## Normalize by relative abundance
micro_table$cal_abund() # Equivalent to total sum scaling (TSS)
class(micro_table$taxa_abund)

## check data structure
sample_1[1:3, ]
otu_1[1:10,]
taxa_1[1,]
print(micro_table)

#### Alpha_Diversity Indcices #######
## Calculate the Alpha-diversity Indices Values.
alpha_diversity <- trans_alpha$new(dataset = micro_table,
                                   group = 'Group')
print(alpha_diversity)

## Return t1$data_stat
alpha_diversity$data_stat[1:10,]

## Calculate the differences among groups using kruskal-walis Rank Sum Test
alpha_diversity$cal_diff(method = 'KW')
?trans_alpha

## Return alpha_diversity differences
alpha_diversity$res_diff |> 
  view()

## Convert to dataframe and save the Results
alpha_diversity_final <- as.data.frame((alpha_diversity$res_diff))

## write.csv(alpha_diversity_final, 'alpha_diversity.csv')

## Calculating differences among groups using ANOVA
alpha_diversity$cal_diff(method = 'anova')

## Return the result
alpha_diversity$res_diff

alpha_1 <- trans_alpha$new(dataset = micro_table, group = 'Group')
alpha_1$cal_diff(method = 'anova')
alpha_1$plot_alpha(measure = 'Chao1')
alpha_1$plot_alpha(measure = 'Chao1', order_x_mean = TRUE)


alpha_1$plot_alpha(pair_compare = TRUE, measure = 'Chao1',
                   shape = 'Group')

alpha_plot_1 <- alpha_1$plot_alpha(measure = 'Chao1', xtext_size = 15)
print(alpha_plot_1)

alpha_plot_1 + theme(axis.line = element_line(size = 0.5), panel.background = element_rect(fill = 'white'),
                     panel.grid.minor = element_line (colour = NA), axis.text.y = element_text(size = 15),
                     axis.title = element_text(size = 15), axis.text.x = element_text(size = 15),
                     legend.text = element_text(size = 13))

###### Beta Diversity #################
## Create a trans-beta object
## Measure parameter to invoke the distance matrix in dataset$beta_diversity 

micro_table$cal_betadiv(method = NULL, unifrac = F)
beta_diversity <- trans_beta$new(dataset = micro_table, group = 'Group',
                                 measure = 'bray')

## Use PCoA as an example, PCA or NMDS is also available
beta_diversity$cal_ordination(method = 'PCoA')

## Plot the PCoA result with confidence ellipse
pcoa_plot <- beta_diversity$plot_ordination(plot_color = 'Group', plot_shape = 'Group',
                                            plot_type = c('point', 'ellipse'))

print(pcoa_plot)

## Modify the plot with ggplot2
pcoa_plot + theme(axis.line = element_line(size = 0.5), panel.background = element_rect(fill = 'white'),
                  panel.grid.minor = element_line (colour = NA), axis.text.y = element_text(size = 15),
                  axis.title = element_text(size = 15), axis.text.x = element_text(size = 15),
                  legend.text = element_text(size = 13))

## Save diversity plot
ggsave('pcoa_plot.png',
       width = 30, height = 15, units = 'cm', limitsize = TRUE)

## Plotting PCA
beta_diversity$cal_ordination(method = 'PCA')

## Plot the PCA Results with confidence ellipse
pca_plot <- beta_diversity$plot_ordination(plot_color = 'Group', plot_shape = 'Group',
                                           plot_type = c('point', 'ellipse'))

print(pca_plot)

## Modify the plot with ggplot2
pca_plot + theme(axis.line = element_line(size = 0.5), panel.background = element_rect(fill = 'white'),
                 panel.grid.minor = element_line (colour = NA), axis.text.y = element_text(size = 15),
                 axis.title = element_text(size = 15), axis.text.x = element_text(size = 15),
                 legend.text = element_text(size = 13))

## Save diversity plot
ggsave('pca_plot.png',
       width = 30, height = 15, units = 'cm', limitsize = TRUE)

#### Comparing Group Distances ####
## Calculate and plot distances within groups
beta_diversity$cal_group_distance()

##Return the values of the group distance
beta_diversity$plot_group_distance(distance_pair_stat = TRUE)


## Clustering
## Use replace_name to set the label name, group parameter to set the color
beta_diversity$plot_clustering(group = 'Group', replace_name = 'Group')


## Save clustering plot
ggsave('bray_curtis_plot.png',
       width = 30, height = 15, units = 'cm', limitsize = TRUE)

### PERMANOVA ####
# manova for all groups when manova_all = TRUE
beta_diversity$cal_manova(manova_all = TRUE)
beta_diversity$res_manova

### Convert to a dataframe and save PERMANOVA for all
permanova_results <- as.data.frame(beta_diversity$res_manova)

write.csv(permanova_results, 'permanova_results.csv')

## manova for each paired groups
beta_diversity$cal_manova(manova_all = FALSE)
beta_diversity$res_manova

### Convert to a dataframe and save PERMANOVA for paired groups
permanova_results_paired <- as.data.frame(beta_diversity$res_manova)

write.csv(permanova_results_paired, 'permanova_paired.csv')

### Clade_label (level 5) represent phylum level in this analysis
## Requirement:ggtree package

## Calculate Relative Abundance
micro_table$cal_abund()

## Save relative abundance
gen_relative_abund <- as.data.frame(micro_table$taxa_abund$Genus)

write.csv(gen_relative_abund, 'Genus_RA.csv')

## Save Phylum Relative Abundance
phy_relative_abund <- as.data.frame(micro_table$taxa_abund$Phylum)

write.csv(phy_relative_abund, 'phylum_RA.csv')

## Merge samples as one community for each group
data_merge <- micro_table$merge_samples(use_group = 'Group')

t1 <- trans_venn$new(micro_table, ratio = NULL)

