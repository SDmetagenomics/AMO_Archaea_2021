library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(indicspecies)
library(superheat)



######## IMPORT AND SET UP ALL THE DATA ##########

'%notin%' <- Negate('%in%')

### Set WD
setwd("~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Analysis/20_11_18_Prot_Clust_Fam_Enrichment_Analysis/")


### Import Gene Annotation Tables
master_anno_filt <- fread("../20_9_15_General_Metabolic_Analysis/Annotations/MASTER/master_annotations_filt_FINAL.txt")
best_genomes <- unique(master_anno_filt$genome)


### Import Detailed Genome Tax and filter for best genomes
genome_tax <- fread("../../Genome_Data/master_genome_tax_table.txt")
genome_tax <- subset(genome_tax, genome %in% best_genomes)


### Generate count of genomes at order taxonomic level and merge in
tax_order_counts <- genome_tax %>% count(tax_Order)
colnames(tax_order_counts) <- c("tax_Order", "tax_Order_counts")
genome_tax <- merge(genome_tax, tax_order_counts, by = "tax_Order")


### Filter out genomes with bad tax + those with < 3 genomes in order
genome_tax <- subset(genome_tax, tax_Domain != "d__Unk")
genome_tax <- subset(genome_tax, tax_Domain == "d__Archaea")
genome_tax <- subset(genome_tax, tax_Order_counts >= 3)


### Generate KEGG Lookup Table
kegg_anno <- data.table(KO = master_anno_filt$KO, KO_anno = master_anno_filt$KO_Anno)
kegg_anno <- kegg_anno[order(kegg_anno$KO),]
kegg_anno <- kegg_anno[!duplicated(kegg_anno$KO),]


### Get Organism AMO Group Variables
AMO_groupings <- fread("AMO_Groupings.txt") 



#### Generate KEGG Count Table by Genome (Only Good Genomes) + Remove outgroup genome and order level taxa < 3 genomes
#### NOTE: We are going to analyze all KOs in this case even if they did not pass threshold but passed e-val threshold

### Generate small table for reformatting + remove NA values 
small_table <- master_anno_filt[,c("ORF", "genome", "KO")]
small_table <- small_table[!is.na(small_table$KO),]
sum(is.na(small_table$fam))


### Filter to just genomes that will be analyzed
small_table_filt <- subset(small_table, genome %in% genome_tax$genome)


### Identify and filter out any KOs that occur less than N (10) times 

## Count KOs 
KO_Counts <- plyr::count(small_table_filt$KO)
KO_Counts_filt <- subset(KO_Counts, freq >= 10)

## Filter KOs from small table using above count filtering
small_table_filt <- subset(small_table_filt, KO %in% KO_Counts_filt$x)

  
  
### Translate small table to long format for all KEGG KOs
KO_hit_table <- dcast(small_table_filt, genome ~ KO,
                       value.var = "KO", fun.aggregate = length)




#### Perform ISA Analaysis 

### Merge genome_tax with subfam_hit_table_ISA   
ISA_input_KO <- merge(genome_tax, KO_hit_table, by = "genome")

### Prepare ISA data
ISA_counts_KO <- ISA_input_KO[,10:ncol(ISA_input_KO)]
ISA_groups_KO <- ISA_input_KO$tax_Order

### Run ISA
ISA_output_KO <- multipatt(x = ISA_counts_KO,
                        cluster = ISA_groups_KO,
                        func = "IndVal.g",
                        restcomb = c(1:25,262), #262 is just the grouping of amo genomes together
                        max.order = 2,
                        control = how(nperm=9999),
                        print.perm = TRUE)




### Generate output summaries and plots

## Generate Tables + FDR correction
ISA_output_KO_tab <- data.table(KO = rownames(ISA_output_KO$sign), ISA_output_KO$sign)
ISA_output_KO_tab$FDR <- p.adjust(ISA_output_KO_tab$p.value, method = "fdr")


## Merge Annotations into ISA output table
# Subfam summaries + remove columns
ISA_output_KO_tab <- merge(ISA_output_KO_tab, kegg_anno, by = "KO")
ISA_output_KO_tab <- ISA_output_KO_tab[,-c(2:26)]


## Filter ISA output to correct FDR and stat 
ISA_output_KO_tab_filt <- subset(ISA_output_KO_tab, FDR <= 0.05 & stat >= 0.4)

### Generate Final Output
ISA_output_KO_tab_final <- subset(ISA_output_KO_tab_filt, index %in% c(14,17,26))

### Write Final Output
#fwrite(ISA_output_KO_tab_final, "~/Desktop/ISA_output_KEGG.txt", sep = "\t")






### Draw Heatmap of Significant KOs

## Set up data
foo <- which(colnames(ISA_input_KO) %in% c("genome", ISA_output_KO_tab_final$KO))
bar <- ISA_input_KO[,..foo]
nex <- as.matrix(bar[,2:ncol(bar)],rownames = bar$genome)
t_nex <- t(nex)

## create a strength variable for side plot if needed 
# strength <- data.frame(KO = colnames(bar[,-1]))
# strength <- merge(strength, ISA_output_KO_tab_final[,c("KO","stat")], by = "KO")

## Create high level taxonomy grouping
ISA_input_KO_2 <- merge(ISA_input_KO, AMO_groupings, by = "genome")

# ## KOs on bottom / Orgs on side
# superheat(sqrt(nex),
#           pretty.order.cols = T,
#           pretty.order.rows = T,
#           heat.na.col = "white",
#           heat.pal = c("white", "steelblue4"),
#           heat.lim = c(0,5),
#           #heat.pal.values = c(0, 0.5, 1))
#           membership.cols = ISA_output_KO_tab_final$index,
#           membership.rows = ISA_input_KO_2$AMO_Group_Phy_2,
#           #col.dendrogram = T,
#           dist.method = "euc",
#           #bottom.label.col = c("","green","blue"),
#           linkage.method = "ward.D2",
#           scale = F,
#           left.label.size = 0.1,
#           left.label.text.size = 3,
#           force.bottom.label = TRUE,
#           bottom.label = "variable",
#           bottom.label.text.size = 3,
#           bottom.label.text.angle = 90,
#           bottom.label.size = 0.2, 
#           #yt = strength$stat,
#           yt.plot.type = "bar")

#png(filename = "KOfams/Significant_KOs_Heatmap.png", width = 10000, height = 6000, res = 600)
superheat(sqrt(t_nex),
          pretty.order.cols = T,
          pretty.order.rows = T,
          heat.na.col = "white",
          heat.pal = c("white", "steelblue4"),
          heat.lim = c(0,5),
          #heat.pal.values = c(0, 0.5, 1))
          membership.rows = ISA_output_KO_tab_final$index,
          membership.cols = ISA_input_KO_2$AMO_Group_Phy_2,
          dist.method = "euc",
          linkage.method = "ward.D2",
          scale = F,
          bottom.label.size = 0.1,
          bottom.label.text.size = 3,
          force.left.label = TRUE,
          left.label = "variable",
          left.label.text.size = 2,
          left.label.size = 0.1,
          left.label.col = "white",
          grid.hline.col = "white",
          grid.vline.col = "white")

#dev.off()




### Test for Statistical Enrichment of KEGG Functional Categories in each Group

### Subset KOs in each cluster into seperate dfs 
index_14 <- subset(ISA_output_KO_tab_final, index == 14)
index_17 <- subset(ISA_output_KO_tab_final, index == 17)
index_26 <- subset(ISA_output_KO_tab_final, index == 26)

### Get Files That Have KO to FunCat Correspondance 
file_list <- list.files("KOfams/KEGG_FunCat_Custom/FunCat_Groups/Final_Sets/")

## Create Output Dataframes
all_nr_annos <- data.table()
funcat_counts <- data.table()

for (i in 1:length(file_list)){
  
  tmp_annos <- fread(paste0("KOfams/KEGG_FunCat_Custom/FunCat_Groups/Final_Sets/",file_list[i]))
  
  anno_name <- as.character(tmp_annos[1,2])
  
  tmp_annos <- tmp_annos[!duplicated(tmp_annos$KO)]
  
  all_nr_annos <- rbind(all_nr_annos, tmp_annos)
  
  # Calculate Counts for each set 
  tot_cnt <- sum(tmp_annos$KO %in% ISA_output_KO_tab$KO)
  
  gp14_cnt <- sum(tmp_annos$KO %in% index_14$KO)
  
  gp17_cnt <- sum(tmp_annos$KO %in% index_17$KO)
  
  gp26_cnt <- sum(tmp_annos$KO %in% index_26$KO)
  
  # Create and append to dataframe
  
  output_i <- data.table(Group = anno_name,
                         Overall_Count = tot_cnt,
                         Overall_Total = 3412,
                         #Overall_Frac = Overall_Count / Overall_Total,
                         gp14_Count = gp14_cnt,
                         gp14_Total = 75,
                         #gp14_Frac = gp14_Count / gp14_Total,
                         gp17_Count = gp17_cnt,
                         gp17_Total = 78,
                         #gp17_Frac = gp17_Count / gp17_Total,
                         gp26_Count = gp26_cnt,
                         gp26_Total = 48)
  #gp26_Frac = gp26_Count / gp26_Total)
  
  funcat_counts <- rbind(funcat_counts, output_i)
  
  
}

### Create a Non-Redundant list of each KO with a single column containing all its functional categories
master_nr_annos <- aggregate(all_nr_annos$Group, list(all_nr_annos$KO), 
                             FUN = function(X) paste(unique(X), collapse="|"))
colnames(master_nr_annos) <- c("KO", "FunCat")


### Merge Functional Categories with ISA Output
ISA_output_KO_tab_final_2 <- merge(ISA_output_KO_tab_final, master_nr_annos, by = "KO", all.x = T)



### Write Important Outputs ***** THE FIRST OUTPUT HERE IS USED TO MANUALLY GENERATE A LONG FORMAT TABLE FOR THE PIE PLOTS
#fwrite(funcat_counts, "KOfams/FunCat_counts.txt", sep = "\t")
#fwrite(ISA_output_KO_tab_final_2, "KOfams/ISA_output_KEGG.txt", sep = "\t")





### NOT INCLUDED IN MANUSCRIPT, FUNCTIONAL CATEGORIES ARE NOT DESCRIPTIVE ENOUGH

### Perform Statistical Testing on Each Functional Category for each group 

# # make funcat_counts a data.frame not data.table
# funcat_counts <- as.data.frame(funcat_counts)
# 
# # create output df
# funcat_sig <- data.table()
# 
# # Run loop with stats 
# for (i in 1:nrow(funcat_counts)){
#   
#   func_name <- as.character(funcat_counts$Group[i])
#   
#   gp_14_fish <- fisher.test(matrix(c(funcat_counts[i,4],
#                                      funcat_counts[i,2] - funcat_counts[i,4],
#                                      funcat_counts[i,5] - funcat_counts[i,4],
#                                      funcat_counts[i,3] - funcat_counts[i,2] - (funcat_counts[i,5] - funcat_counts[i,4])), 2, 2))
#   
#   gp_17_fish <- fisher.test(matrix(c(funcat_counts[i,6],
#                                      funcat_counts[i,2] - funcat_counts[i,6],
#                                      funcat_counts[i,7] - funcat_counts[i,6],
#                                      funcat_counts[i,3] - funcat_counts[i,2] - (funcat_counts[i,7] - funcat_counts[i,6])), 2, 2))
#   
#   gp_26_fish <- fisher.test(matrix(c(funcat_counts[i,8],
#                                      funcat_counts[i,2] - funcat_counts[i,8],
#                                      funcat_counts[i,9] - funcat_counts[i,8],
#                                      funcat_counts[i,3] - funcat_counts[i,2] - (funcat_counts[i,9] - funcat_counts[i,8])), 2, 2))
#   
#   tmp_out <- data.table(Group = func_name,
#                         gp14_p = gp_14_fish$p.value,
#                         gp14_or = gp_14_fish$estimate,
#                         gp17_p = gp_17_fish$p.value,
#                         gp17_or = gp_17_fish$estimate,
#                         gp26_p = gp_26_fish$p.value,
#                         gp26_or = gp_26_fish$estimate)
#   
#   
#   funcat_sig <- rbind(funcat_sig, tmp_out)
#   
# }








##### DO PLOTTING FOR FIG 4C 



# ### Dumbell Plots to show difference of each FunCat in each group from the background frequency
# funcat_frac_long <- fread("~/Desktop/TMP_Fig4/Enrichment_and_FunCat/KO_funcat_frac_long.txt")
# 
# ggplot(funcat_frac_long) +
#   geom_segment(aes(x=FunCat, xend=FunCat, y=Group_Frac, yend=Overall_Frac), color="grey") +
#   geom_point(aes(x=FunCat, y=Group_Frac), size = 3,fill = "steelblue4", shape = 21) +
#   geom_point(aes(x=FunCat, y=Overall_Frac), size = 3,  color = "grey") +
#   facet_wrap(.~Group, ncol = 1) +
#   theme_bw() +
#   xlab(NULL) +
#   ylab("Fraction of KOs in Category") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))



### PIE Plots from counts of FunCat observations 
### Long format table was created manually from FunCat_counts.txt
long_format_counts <- fread("KOfams/FunCat_counts_long.txt")



##**** NOTE: Saved as 6x8 landscape PDF

gp14_plot <- subset(long_format_counts, Group == "Gp 14")

ggplot(gp14_plot, aes(x="", y = Count, fill = FunCat)) +
  geom_bar(stat="identity", width=1, color = "black") +
  coord_polar("y", start=0) +
  #geom_text(aes(label = paste0(gp14_Count, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  #facet_wrap(.~Group, ncol = 1) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_brewer(palette = "Paired")

gp17_plot <- subset(long_format_counts, Group == "Gp 17")

ggplot(gp17_plot, aes(x="", y = Count, fill = FunCat)) +
  geom_bar(stat="identity", width=1, color = "black") +
  coord_polar("y", start=0) +
  #geom_text(aes(label = paste0(gp14_Count, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  #facet_wrap(.~Group, ncol = 1) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_brewer(palette = "Paired")



gp26_plot <- subset(long_format_counts, Group == "Gp 26")

ggplot(gp26_plot, aes(x="", y = Count, fill = FunCat)) +
  geom_bar(stat="identity", width=1, color = "black") +
  coord_polar("y", start=0) +
  #geom_text(aes(label = paste0(gp14_Count, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  #facet_wrap(.~Group, ncol = 1) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_brewer(palette = "Paired")


# ggplot(long_format_counts, aes(x="", y = Count, fill = FunCat)) +
#   geom_bar(stat="identity", width=1, color = "black") +
#   coord_polar("y", start=0) +
#   #geom_text(aes(label = paste0(gp14_Count, "%")), position = position_stack(vjust=0.5)) +
#   labs(x = NULL, y = NULL, fill = NULL) +
#   theme_classic() +
#   facet_wrap(.~Group, ncol = 1, scales = "free") +
#   theme(axis.line = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank()) +
#   scale_fill_brewer(palette = "Paired")
