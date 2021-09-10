library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)


### Set Working Directory + ID files
setwd("~/Dropbox/Banfield_Lab_Files/Manuscripts/2021_AMO_Archaea/Submissions/ISME/Submission_2/Revision_Work/amoA_pmoA_graftM/mapping_runs/")
file_list <- sub("./", "", list.dirs(recursive = F))

### Get rpL6 Abundance Data
rpL6_dat <- fread("~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Analysis/20_9_25_Cross_Biome_Mapping/Abundance_Analysis/Aggregated_Abundance_Table.txt")

### Metadata
metadata <- fread("~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Analysis/20_9_25_Cross_Biome_Mapping/Abundance_Analysis/Master_Metadata.txt")



####### AMO ANALYSIS ##########

### Build Master Tax List
Master_Tax <- vector(mode = "character")

for (i in 1:length(file_list)){
  
  foo <- fread(paste0(file_list[i],"/combined_count_table.txt"))
  Master_Tax <- c(Master_Tax, foo$ConsensusLineage)
  
}

Master_Tax <- data.table(ConsensusLineage = unique(Master_Tax))

## Get Individual Tax Columns
Master_Tax <- separate(data = Master_Tax, # split tax into individual columns
                      col = "ConsensusLineage",
                      into =  c("Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                      sep = "; ",
                      fill = "right",
                      remove = F)



### Build Master Abundance Table
Master_Abund <- Master_Tax

for (i in 1:length(file_list)){
  
  foo <- fread(paste0(file_list[i],"/combined_count_table.txt"))
  foo <- foo[,-1] # remove ID number column
  Master_Abund <- merge(Master_Abund, foo, by = "ConsensusLineage", all.x = T)
  
}


### Add Custom Tax Column
Master_Abund$Custom_Tax <- NA
Master_Abund$Custom_Tax <- ifelse(Master_Abund$Order == "o__Nitrososphaerales", "Nitrososphaerales", Master_Abund$Custom_Tax)
Master_Abund$Custom_Tax <- ifelse(Master_Abund$Order == "o__RBG-16-68-12", "Angelarchaeales", Master_Abund$Custom_Tax)
Master_Abund$Custom_Tax <- ifelse(Master_Abund$Kingdom == "k__Bacteria", "Bacteria", Master_Abund$Custom_Tax)
Master_Abund$Custom_Tax <- ifelse(is.na(Master_Abund$Phylum), "Unknown", Master_Abund$Custom_Tax)


### Aggregate data
Master_Abund_Agg <- Master_Abund[,-c(1:9)]
Master_Abund_Agg <- melt(Master_Abund_Agg)
Master_Abund_Agg[is.na(Master_Abund_Agg)] <- 0
Master_Abund_Agg <- data.table(Master_Abund_Agg %>%
                                 group_by(Custom_Tax, variable) %>%
                                 summarise(counts = sum(value)))

### Merge in metadata
Master_Abund_Agg <- merge(Master_Abund_Agg, metadata, by.x = "variable", by.y = "Map_Sample", all.x = T)


### Create Normalized Data Column
Master_Abund_Agg$Norm_counts <- Master_Abund_Agg$counts / Master_Abund_Agg$Reads



### All raw cts plot
ggplot(Master_Abund_Agg, aes(x = Sample_Name, y = counts, fill = Custom_Tax)) +
  geom_bar(stat = "identity") +
  #scale_y_sqrt() +
  xlab(NULL) +
  ylab("Total Amo/Pmo Reads") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1),
        panel.grid.major.x = element_blank()) 

### All norm cts plot
ggplot(Master_Abund_Agg, aes(x = Sample_Name, y = Norm_counts, fill = Custom_Tax)) +
  geom_bar(stat = "identity") +
  #scale_y_sqrt() +
  xlab(NULL) +
  ylab("Normalized Amo/Pmo Reads") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1),
        panel.grid.major.x = element_blank()) 


### Filtered for only archaea plot normalized counts
Master_Abund_Agg_onlyArch <- subset(Master_Abund_Agg, Custom_Tax %in% c("Angelarchaeales", "Nitrososphaerales"))

ggplot(Master_Abund_Agg_onlyArch, aes(x = Sample_Name, y = Norm_counts, fill = Custom_Tax)) +
  geom_bar(stat = "identity") +
  #scale_y_sqrt() +
  xlab(NULL) +
  ylab("Normalized AmoA/PmoA Reads") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2e-05)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none") 

#ggsave("../amoA_pmoA_graftM_counts_plot.pdf", width = 12, height = 4)



####### rpL6 ANALYSIS ##########

### Reformat taxonomy in rpL6_dat
colnames(rpL6_dat)[1] <- "ConsensusLineage"

rpL6_dat <- separate(data = rpL6_dat, # split tax into individual columns
                     col = "ConsensusLineage",
                     into =  c("Kingdom", "Phylum", "Class", "Order"), # note taxonomy of this file has already been pruned a little
                     sep = ";",
                     fill = "right",
                     remove = F)

### Keep Order Lineage, melt, and aggregate data
rpL6_dat_Agg <- rpL6_dat[,-c(1:4)]
rpL6_dat_Agg <- melt(rpL6_dat_Agg)
#Master_Abund_Agg[is.na(Master_Abund_Agg)] <- 0
rpL6_dat_Agg <- data.table(rpL6_dat_Agg %>%
                                 group_by(Order, variable) %>%
                                 summarise(counts = sum(value)))

### Merge in metadata
rpL6_dat_Agg <- merge(rpL6_dat_Agg, metadata, by.x = "variable", by.y = "Map_Sample", all.x = T)

### Filter to Angel + Nitrososphaerales
rpL6_dat_Agg_sub <- subset(rpL6_dat_Agg, Order %in% c("o__Nitrososphaerales", "o__RBG-16-68-12"))

### Create Normalized Data Column
rpL6_dat_Agg_sub$Norm_counts <- rpL6_dat_Agg_sub$counts / rpL6_dat_Agg_sub$Reads

### Rename Orders
rpL6_dat_Agg_sub$Order <- ifelse(rpL6_dat_Agg_sub$Order == "o__RBG-16-68-12", "Angelarchaeales",rpL6_dat_Agg_sub$Order)
rpL6_dat_Agg_sub$Order <- ifelse(rpL6_dat_Agg_sub$Order == "o__Nitrososphaerales", "Nitrososphaerales",rpL6_dat_Agg_sub$Order)

#### Plotting

### All raw cts plot
ggplot(rpL6_dat_Agg_sub, aes(x = Sample_Name, y = counts, fill = Order)) +
  geom_bar(stat = "identity") +
  #scale_y_sqrt() +
  xlab(NULL) +
  ylab("Total rpL6 Reads") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1),
        panel.grid.major.x = element_blank()) 

### All norm cts plot
ggplot(rpL6_dat_Agg_sub, aes(x = Sample_Name, y = Norm_counts, fill = Order)) +
  geom_bar(stat = "identity") +
  #scale_y_sqrt() +
  xlab(NULL) +
  ylab("Normalized rpL6 Reads") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2e-05)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none") 

#ggsave("../rpL6_graftM_counts_plot.pdf", width = 12, height = 4)



##### Perform Correlations Between Data

### Angel Cor
Master_Angel_pmo <- subset(Master_Abund_Agg_onlyArch, Custom_Tax == "Angelarchaeales")
Master_Angel_rpL6 <- subset(rpL6_dat_Agg_sub, Order == "Angelarchaeales")
Angel_cor_dat <- data.table(Sample = Master_Angel_pmo$variable,
                            Site = Master_Angel_pmo$Site,
                            Norm_pmo_cts = Master_Angel_pmo$Norm_counts,
                            Norm_rpL6_cts = Master_Angel_rpL6$Norm_counts)

cor.test(Angel_cor_dat$Norm_pmo_cts, Angel_cor_dat$Norm_rpL6_cts, method = "spearman")

ggplot(Angel_cor_dat, aes(x = Norm_pmo_cts, y = Norm_rpL6_cts, fill = Site)) +
  geom_point(shape = 21, size = 3, alpha = 0.8) +
  scale_x_sqrt() +
  scale_y_sqrt() +
  scale_fill_brewer(palette = "Set3")


### Nitroso Cor
Master_Nitroso_pmo <- subset(Master_Abund_Agg_onlyArch, Custom_Tax == "Nitrososphaerales")
Master_Nitroso_rpL6 <- subset(rpL6_dat_Agg_sub, Order == "Nitrososphaerales")
Nitroso_cor_dat <- data.table(Sample = Master_Nitroso_pmo$variable,
                            Site = Master_Nitroso_pmo$Site,
                            Norm_pmo_cts = Master_Nitroso_pmo$Norm_counts,
                            Norm_rpL6_cts = Master_Nitroso_rpL6$Norm_counts)

cor.test(Nitroso_cor_dat$Norm_pmo_cts, Nitroso_cor_dat$Norm_rpL6_cts, method = "spearman")

ggplot(Nitroso_cor_dat, aes(x = Norm_pmo_cts, y = Norm_rpL6_cts, fill = Site)) +
  geom_point(shape = 21, size = 3, alpha = 0.8) +
  #scale_x_sqrt() +
  #scale_y_sqrt() +
  scale_fill_brewer(palette = "Set3")







######## Analyze Data From rpL6 Correlation Network 

#### Load Data

## Pairwise Correlations Angel
rpL6_prop_Angel <- fread("~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Analysis/20_9_25_Cross_Biome_Mapping/Abundance_Analysis/L6_Proportionality_filt_Angel.txt")

## Pairwise Correlations Nitroso
rpL6_prop_Nitroso <- fread("~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Analysis/20_9_25_Cross_Biome_Mapping/Abundance_Analysis/L6_Proportionality_filt_Nitroso.txt")
rpL6_prop_Nitroso <- subset(rpL6_prop_Nitroso, Partner_Phylum != "Unknown")

## Node Statistics
rpL6_prop_node_stats <- read.csv("../L6_graftM_network/node_list_wStats.csv")

## Custom colors for Phyla
phy_colors <- c(Acidobacteriota = "#8DC63F", Actinobacteriota = "#1B75BC", Armatimonadota = "#BE1E2D", Bacteroidota = "#00A79D",
                Chloroflexota = "#F16447", "CSP1-3" = "#ED1C24", Desulfobacterota = "#FFFFB3", Dormibacterota = "#006838", Gemmatimonadota = "#EC008C",
                Methylomirabilota = "#92278F", Myxococcota = "#774E24", Nitrospirota = "#32C4E8", Patescibacteria = "#E3C0D3",
                Planctomycetota = "#90AAE1", Proteobacteria = "#C2B59B", Thermoplasmatota = "#F9ED32", Thermoproteota = "#FBB040",
                Verrucomicrobiota = "#39B54A")


### Plot Associations - Angel
ggplot(rpL6_prop_Angel, aes(x = reorder(Partner, -Propr), y = Propr, fill = Partner_Phylum)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  xlab(NULL) +
  ylab("Rho") +
  scale_fill_manual(values = phy_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
                                   panel.grid.major.x = element_blank()) 

#ggsave("../L6_graftM_network/Angel_Associations.pdf", width = 12, height = 6)

### Plot Associations - Nitroso
ggplot(rpL6_prop_Nitroso, aes(x = reorder(Partner, -Propr), y = Propr, fill = Partner_Phylum)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  xlab(NULL) +
  ylab("Rho") +
  scale_fill_manual(values = phy_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()) 

#ggsave("../L6_graftM_network/Nitroso_Associations.pdf", width = 12, height = 6)


#### Plot Network Node Properties 

## Closeness centrality
ggplot(rpL6_prop_node_stats, aes(x = as.factor(modularity_class), y = closnesscentrality)) +
  geom_boxplot() +
  geom_jitter() 

## Bridging centrality
ggplot(rpL6_prop_node_stats, aes(x = as.factor(modularity_class), y = bridgingcentrality)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_bw() +
  xlab("Sub-Network Module") +
  ylab("Bridging Centrality")

#ggsave("../L6_graftM_network/Bridging_Centrality.pdf", width = 3, height = 3)

# First check for global significnat difference 
kruskal.test(bridgingcentrality ~ modularity_class, data = rpL6_prop_node_stats) ## Is Significant
# Next check for pairwise significant differences 
bar <- pairwise.wilcox.test(rpL6_prop_node_stats$bridgingcentrality, rpL6_prop_node_stats$modularity_class,
                            p.adjust.method = "fdr")
bar


## Betweeness centrality 
ggplot(rpL6_prop_node_stats, aes(x = as.factor(modularity_class), y = betweenesscentrality)) +
  geom_boxplot() +
  geom_jitter() 

# First check for global significnat difference 
kruskal.test(betweenesscentrality ~ modularity_class, data = rpL6_prop_node_stats) ## Is Significant
# Next check for pairwise significant differences 
bar <- pairwise.wilcox.test(rpL6_prop_node_stats$betweenesscentrality, rpL6_prop_node_stats$modularity_class,
                            p.adjust.method = "fdr")
bar
