## Load the necessary Packages

library(tidyverse)
library(microeco)
library(dplyr)
library(janitor)
library(agricolae)
library(ggplot2)

## Import and Read the Dataset (OTU and Sample data) 

otu_1 <- read.csv('Microbiome Analysis/otu_data.csv', 
                  row.names = 1, check.names = FALSE)

sample_1 <- read.csv('Microbiome Analysis/sample_data.csv',
                     row.names = 1, check.names = FALSE)

taxa_1 <- read.csv('Microbiome Analysis/taxa_data.csv',
                   row.names = 1, check.names = FALSE)

### View the imported datasets
taxa_1 |> 
  glimpse()

taxa_1 |> 
  view()

sample_1 |> 
  glimpse()

sample_1 |> 
  view()

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

micro_table$tidy_dataset()

## Normalize by relative abundance
micro_table$cal_abund() # Equivalent to total sum scaling (TSS)

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