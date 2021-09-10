library(ggplot2)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(multcompView)
library(rcompanion)




### Function for using multicompview see line 308
tri.to.squ <- function(x){
  rn<-row.names(x)
  cn<-colnames(x)
  an<-unique(c(cn,rn))
  myval<-x[!is.na(x)]
  mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}


######## IMPORT AND SET UP ALL THE DATA ##########

'%notin%' <- Negate('%in%')

### Set WD
setwd("~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Analysis/20_9_14_Cu_Protein_Analysis/")

### Import Gene Annotation Tables
#master_anno <- fread("../20_9_15_General_Metabolic_Analysis/Annotations/MASTER/master_annotations_all.txt")
master_anno_filt <- fread("../20_9_15_General_Metabolic_Analysis/Annotations/MASTER/master_annotations_filt_fixed.txt")
best_genomes <- unique(master_anno_filt$genome)

### Import Detailed Genome Tax and filter for best genomes
genome_tax <- fread("../../Genome_Data/master_genome_tax_table.txt")
genome_tax <- subset(genome_tax, genome %in% best_genomes)

## Generate count of genomes at order taxonomic level and merge in
tax_order_counts <- genome_tax %>% count(tax_Order)
colnames(tax_order_counts) <- c("tax_Order", "tax_Order_counts")
genome_tax <- merge(genome_tax, tax_order_counts, by = "tax_Order")

## Filter out genomes with bad tax 
genome_tax <- subset(genome_tax, tax_Domain != "d__Unk")
genome_tax <- subset(genome_tax, tax_Domain == "d__Archaea")

### Import Subfam Summary Tables
#subfam_summary <- fread("../20_9_15_General_Metabolic_Analysis/Annotations/MASTER/subfam_summary_all.txt")
subfam_summary_filt <- fread("../20_9_15_General_Metabolic_Analysis/Annotations/MASTER/subfam_summary_filt.txt")

### Import Subfam Hit Tables
fam_hit_table_filt <- fread("../20_9_15_General_Metabolic_Analysis/Annotations/MASTER/fam_hit_table_filt.txt")
subfam_hit_table_filt <- fread("../20_9_15_General_Metabolic_Analysis/Annotations/MASTER/subfam_hit_table_filt.txt")

### Import Has AMO Table

### Generate Annotation Tables 

## KEGG
#kegg_anno <- data.table(KO = master_anno_filt$KO, KO_anno = master_anno_filt$KO_Anno)
#kegg_anno <- kegg_anno[order(kegg_anno$KO),]
#kegg_anno <- kegg_anno[!duplicated(kegg_anno$KO),]

## arCOG
#arcog_anno <- data.table(arCOG = master_anno_filt$arCOG, arCOG_FC = master_anno_filt$arCOG_FC, arCOG_anno = master_anno_filt$arCOG_Anno)
#arcog_anno <- arcog_anno[order(arcog_anno$arCOG),]
#arcog_anno <- arcog_anno[!duplicated(arcog_anno$arCOG),]
  









###### Assessment of fam00018 as a proxy for BCPs - Looks only at proteins with domain level annotations ######

### 01 ### We must determine if there is a protein fam that is representative of all BCPs in our analysis 
### 01 ### First we find all proteins that have Cu-binding like domains from pfam and assess how this is distirbuted at the family level

## Name domains from copper_bind clan to be searched
search_domains <- c("Copper-bind", "COX2", "COX_ARM", "Cu-oxidase", "Cu-oxidase_2", "Cu-oxidase_3", "Cu_bind_like", "Cupredoxin_1", "CzcE", "DP-EP", "Ephrin", "hGDE_N", "PAD_N", "PixA", "SoxE")

## Run search for loop - loop is used to search for each pattern individually using %like%
search_out <- data.table()

for (i in 1:length(search_domains)){
  tmp <- master_anno_filt[master_anno_filt$PFAMs %like% search_domains[i],]
  search_out <- rbind(search_out, tmp)
}

## Remove Duplicated ORFs + Create Tabular Output
search_out <- search_out[!duplicated(search_out$ORF),]
BCP_prot_total <- nrow(search_out)

## Generate counts of hit genes in each fam
BCP_count_by_fam <- plyr::count(search_out$fam)  ### There are 100 proteins with BCP domains not in a fam (NA)
colnames(BCP_count_by_fam) <- c("fam","BCP_count")

## Gather all genes from fams where at least one hit to above domains was registered
fam_grab <- na.omit(unique(search_out$fam))
all_BCP_fam_prot <- master_anno_filt[master_anno_filt$fam %in% fam_grab,]

## Generate Total number of proteins in each fam with a BCP hit and merge into BCP_count_by_fam
all_BCP_fam_counts <- plyr::count(all_BCP_fam_prot$fam)
colnames(all_BCP_fam_counts) <- c("fam","total_count")
BCP_count_by_fam <- merge(BCP_count_by_fam, all_BCP_fam_counts, by = "fam", all.x = T)
BCP_count_by_fam$NotBCP_count <- BCP_count_by_fam$total_count - BCP_count_by_fam$BCP_count
BCP_count_by_fam$BCP_frac <- BCP_count_by_fam$BCP_count/BCP_count_by_fam$total_count
BCP_count_by_fam$fam00018 <- ifelse(BCP_count_by_fam$fam == "fam00018", TRUE, FALSE)
BCP_count_by_fam <- as.data.table(BCP_count_by_fam)

## Plot number of BCP in each fam - PLOT 1
ggplot(na.omit(BCP_count_by_fam), aes(x = reorder(fam, -BCP_count), y = log10(BCP_count + 1))) +
  geom_bar(stat = "identity", fill = "steelblue4") +
  #scale_y_log10() +
  xlab(NULL) +
  ylab("Number of BCP Domain Containing Proteins (Log10 Scale)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Plot Fraction of BCP in each fam - PLOT 2
#ggplot(BCP_count_by_fam, aes(x = reorder(fam, -BCP_frac), y = BCP_frac)) +
#  geom_bar(stat = "identity", fill = "steelblue4") +
#  xlab(NULL) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Generate plot showing colored fraction of total protein count that has BCP domains
BCP_col_frac_plot <- melt(BCP_count_by_fam, id.vars = c("fam"), measure.vars = c("BCP_count", "NotBCP_count"))
BCP_col_frac_plot <- na.omit(BCP_col_frac_plot)

## Plot Fraction of BCP in each fam - PLOT 3
ggplot(BCP_col_frac_plot, aes(x = reorder(fam, -value), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab(NULL) +
  ylab("Number of Proteins") +
  scale_fill_manual(values = c("firebrick3", "steelblue4")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## For Proteins in Fams Plot BCP in fam00018 vs Not in fam00018 - PLOT 4
ggplot(na.omit(BCP_count_by_fam), aes(x = fam00018, y = BCP_count, fill = fam00018)) +
  geom_bar(stat = "identity") +
  xlab("In fam00018") +
  ylab("Number of BCP Containing Proteins") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("steelblue4", "firebrick3"))
  
aggregate(BCP_count ~ fam00018, BCP_count_by_fam, sum) ### fam00018 Fraction is 90.2%


### 01 ### Explanation:
### 01 ### Using the above code we have looked for all the PFAM domains that are associated with BCPs
### 01 ### we then figured out and plotted the number of BCP proteins that were found in each protein family where
### 01 ### at least one was identified (PLOT 1). We also plotted the fraction of each protein family that is made up of BCP domain
### 01 ### containing protens (PLOT 2).  At this point we can see that fam00018 has many BCP domain containing proteins and a very 
### 01 ### large fraction of the proteins in this family are those which contain BCP domains. Other families either have very few
### 01 ### BCP domain containing proteins relative to their total size or are rather small families. Next we plot the total 
### 01 ### number of proteins in each family and color based on how many are BCP domain containing protens (PLOT 3). Here it is very
### 01 ### obvious that fam00018 contains many BCP domain contating proteins and other fams that contain them do not. Finally we calculate
### 01 ### for all BCP containing proteins that are in fams how many are in fam00018 vs other fams (PLOT 4). We find that > 90 % of 
### 01 ### BCP containing proteins assigned to fams are in fam00018. This makes fam00018 a good proxy for this analysis.




### 02 ### We have now established that fam00018 is the one we want to analyze because it contains all the BCPs
### 02 ### Second we show the domain distribution within fam00018

## Get all proteins in fam00018
fam00018_genes_filt <- subset(master_anno_filt, fam == "fam00018")
fam00018_genes_filt <- merge(fam00018_genes_filt, genome_tax, by = "genome", all.x = T)

## Plot Number of fam00018 proteins with each identified domain architecture - PLOT 1
fam00018_genes_domfreq <- plyr::count(fam00018_genes_filt$PFAMs)
colnames(fam00018_genes_domfreq) <- c("Pfam", "Prot_Count")
fam00018_genes_domfreq[is.na(fam00018_genes_domfreq)] <- "No_Annotation"
# all - PLOT 1a
ggplot(fam00018_genes_domfreq, aes(x = reorder(Pfam, -Prot_Count), y = Prot_Count)) +
  geom_bar(stat = "identity", fill = "steelblue4") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# architectures in > 5 proteins - PLOT 1b
ggplot(fam00018_genes_domfreq[fam00018_genes_domfreq$Prot_Count >= 5,], aes(x = reorder(Pfam, -Prot_Count), y = Prot_Count)) +
  geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
  xlab("Protein Domain Architecture") +
  ylab("Number of fam00018 Proteins") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

## Plot Domain Distribution for each subfam - PLOT 2
fam00018_subfam_domfreq <- data.table(fam00018_genes_filt %>% 
                                        group_by(subfam, PFAMs) %>%
                                        summarise(prot_count = n()))
fam00018_subfam_domfreq[is.na(fam00018_subfam_domfreq)] <- "No_Annotation"
fam00018_subfam_domfreq <- merge(fam00018_subfam_domfreq, subfam_summary_filt, by = "subfam")
# all - PLOT 2a
ggplot(fam00018_subfam_domfreq, aes(x = subfam, y = prot_count, fill = PFAMs)) +
  geom_bar(stat = "identity", position = "fill") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# subfams with > 5 proteins PLOT 2b
ggplot(fam00018_subfam_domfreq[fam00018_subfam_domfreq$count >= 5,], aes(x = subfam, y = prot_count, fill = PFAMs)) +
  geom_bar(stat = "identity", position = "fill") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Plot Domain Count for each subfam - PLOT 3 (Has NA removed and should be zero if no domains are found)
subfam_grab <- unique(fam00018_genes_filt$subfam)
fam00018_subfam_domfreq2 <- subset(subfam_summary_filt, subfam %in% subfam_grab)
# all - PLOT 3a
ggplot(fam00018_subfam_domfreq2, aes(x = reorder(subfam, -PFAM_num), y = PFAM_num)) +
  geom_bar(stat = "identity") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# subfams with > 5 proteins - PLOT 3b
ggplot(fam00018_subfam_domfreq2[fam00018_subfam_domfreq2$count >= 5,], aes(x = reorder(subfam, -PFAM_num), y = PFAM_num)) +
  geom_bar(stat = "identity") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



### 02 ### Explanation:
### 02 ### Using the above code we first get all proteins that are in fam00018. Next we count all the unique domain
### 02 ### architectures that are found for each protein and plot them (PLOT 1). Here we see that the vast majority of 
### 02 ### proteins only contain the Copper-bind domain, many have no annotation, and the remainder are either COX2 or  
### 02 ### other types of single / multi domain Cu oxidase or cupredoxin. This allows us to see the composition of fam00018
### 02 ### as one associated with basically 3 things (1 domain Cu protiens, 2-3 domain Cu proteins, and COX2 proteins). Next 
### 02 ### we want to investigate the purity of domain annotations in the subfams, so we count the number of unique domain 
### 02 ### architectures in each subfam and plot them as a part of a whole (PLOT 2). We also plot the number of unique domains found
### 02 ### in each subfam after excluding No_Annotation (PLOT 3). Using PLOT 2 and PLOT 3 it is possible to see that most subfams contain
### 02 ### either one or 2 PFAM architectures, and often those architectures are related.




### 03 ### We have now established that fam00018 contains many subfams with BCP domains and that those subfams are reasonably pure in domain architecture
### 03 ### Third we want to look at the distribution of fam00018 proteins across taxonomic lineages in our study

## Get all proteins in fam00018
fam00018_genes_filt <- subset(master_anno_filt, fam == "fam00018")
fam00018_genes_filt <- merge(fam00018_genes_filt, genome_tax, by = "genome", all.x = T)

## Count number of fam00018 proteins in each genome that has them, merge into full taxonomy of genomes, add zeros for those that have none
fam00018_genes_genomecount <- plyr::count(fam00018_genes_filt$genome)
colnames(fam00018_genes_genomecount) <- c("genome","fam00018_count")
fam00018_genes_genomecount <- merge(genome_tax, fam00018_genes_genomecount, by = "genome", all.x = T)
fam00018_genes_genomecount[is.na(fam00018_genes_genomecount)] <- 0

## Counts for all Order level taxa together (groups w/ >= 3 genomes) - PLOT 1
ggplot(fam00018_genes_genomecount[fam00018_genes_genomecount$tax_Order_counts >= 3,],
       aes(x = reorder(tax_Order, -fam00018_count), y = fam00018_count, fill = tax_Phylum)) +
  geom_jitter(width = .1, size = 1, alpha = 0.8, height = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  xlab(NULL) +
  ylab("Fam00018 Proteins (Counts)") +
  scale_fill_manual(values = c("steelblue4", "firebrick3")) +
  #theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

## Counts for Thermoplasmatota  - PLOT 2
ggplot(fam00018_genes_genomecount[fam00018_genes_genomecount$tax_Order_counts >= 3 & fam00018_genes_genomecount$tax_Phylum == "p__Thermoplasmatota",], 
       aes(x = reorder(tax_Order, -fam00018_count), y = fam00018_count)) +
  geom_jitter(width = .1, size = 2, alpha = 0.6, height = 0) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "steelblue4") +
  xlab(NULL) +
  ylab("Fam00018 Proteins (Counts)") +
  ylim(c(0, 15)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#ggsave("~/Desktop/TMP_Fig3/Cu_Prot_Statistics/Fig3A_v1.pdf", width = 6, height = 4)

## Counts for Thermoproteota  - PLOT 3
ggplot(fam00018_genes_genomecount[fam00018_genes_genomecount$tax_Order_counts >= 3 & fam00018_genes_genomecount$tax_Phylum == "p__Thermoproteota",], 
       aes(x = reorder(tax_Order, -fam00018_count), y = fam00018_count)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "firebrick3") +
  geom_jitter(width = .1, size = 1.5, alpha = 0.6, height = 0) +
  xlab(NULL) +
  ylab("Fam00018 Proteins (Counts)") +
  ylim(c(0,35)) +
  #scale_fill_manual(values = c("steelblue4")) +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 




##########$$$ REVIEW UPDATE 
##########$$$


## Counts for Thermoproteota Order  - PLOT 4
Nitroso_compare <- subset(fam00018_genes_genomecount, tax_Order == "o__Nitrososphaerales")
#fwrite(Nitroso_compare, "~/Dropbox/Banfield_Lab_Files/Manuscripts/2021_AMO_Archaea/Submissions/ISME/Submission_2/Revision_Work/Fam00018_work/Nitrososphaerales_Fam00018_summary.txt", sep = "\t")
Nitroso_hasAMO <- fread("~/Dropbox/Banfield_Lab_Files/Manuscripts/2021_AMO_Archaea/Submissions/ISME/Submission_2/Revision_Work/Fam00018_work/Nitrososphaerales_Fam00018_summary_HasAMO.txt")


ggplot(Nitroso_compare,aes(x = reorder(tax_Family, -fam00018_count), y = fam00018_count)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "firebrick3") +
  geom_jitter(width = .1, size = 1.5, alpha = 0.6, height = 0) +
  xlab(NULL) +
  ylab("Fam00018 Proteins (Counts)") +
  ylim(c(0,35)) +
  #scale_fill_manual(values = c("steelblue4")) +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("~/Dropbox/Banfield_Lab_Files/Manuscripts/2021_AMO_Archaea/Submissions/ISME/Submission_2/Revision_Work/Fam00018_work/Nitroso_fam00018_counts.pdf", width = 6, height = 4)

ggplot(Nitroso_hasAMO,aes(x = reorder(Has_AMO_Fam, -fam00018_count), y = fam00018_count, fill = tax_Family )) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "grey") +
  geom_jitter(shape = 21, width = .1, size = 2, alpha = 0.6, height = 0) +
  xlab(NULL) +
  ylab("Fam00018 Proteins (Counts)") +
  ylim(c(0,35)) +
  scale_fill_brewer(palette = "Set1") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("~/Dropbox/Banfield_Lab_Files/Manuscripts/2021_AMO_Archaea/Submissions/ISME/Submission_2/Revision_Work/Fam00018_work/Nitroso_AMO_Compare.pdf", width = 6, height = 4)

wilcox.test(fam00018_count ~ Has_AMO_Fam, data = Nitroso_hasAMO)



## Statistical testing for differences between fam00018 genome counts by order
# Remove all orders with less than 3 genomic representatives
stat_input_1 <- subset(fam00018_genes_genomecount, tax_Order_counts >= 3)

## Do testing across all orders - STAT 1
# First check for global signficnat differences 
kruskal.test(fam00018_count ~ tax_Order, data = stat_input_1) ## Is Significant
# Next check for pairwise significant differences 
pairwise.wilcox.test(stat_input_1$fam00018_count, stat_input_1$tax_Order,
                     p.adjust.method = "fdr")

## Do testing across THERMOPLASMATOTA - STAT 2
stat_input_2 <- subset(stat_input_1, tax_Phylum == "p__Thermoplasmatota")
# First check for global signficnat differences 
kruskal.test(fam00018_count ~ tax_Order, data = stat_input_2) ## Is Significant
# Next check for pairwise significant differences 
foo <- pairwise.wilcox.test(stat_input_2$fam00018_count, stat_input_2$tax_Order,
                     p.adjust.method = "fdr")
foo
multcompLetters(tri.to.squ(foo$p.value))

## Do testing across THERMOPROTEOTA - STAT 3
stat_input_3 <- subset(stat_input_1, tax_Phylum == "p__Thermoproteota")
# First check for global signficnat differences 
kruskal.test(fam00018_count ~ tax_Order, data = stat_input_3) ## Is Significant
# Next check for pairwise significant differences 
foo <- pairwise.wilcox.test(stat_input_3$fam00018_count, stat_input_3$tax_Order,
                     p.adjust.method = "fdr")
foo
multcompLetters(tri.to.squ(foo$p.value))



### 03 ### Explanation:
### 03 ### 
### 03 ### 
### 03 ### 
### 03 ### 
### 03 ### 
### 03 ### 
### 03 ###
### 03 ### 
### 03 ###




### 04 ### We have now established that fam00018 proteins are significantly enriched in both Angelarcheales and Nitrosospheria
### 04 ### Fourth we will look at the taxonomic distributions of types of copper proteins across order level lineages

## Get all proteins in fam00018
fam00018_genes_filt <- subset(master_anno_filt, fam == "fam00018")
fam00018_genes_filt <- merge(fam00018_genes_filt, genome_tax, by = "genome", all.x = T)

## Get Subfam Summaries for those within fam00018
subfam_grab <- unique(fam00018_genes_filt$subfam)
fam00018_subfam_summaries <- subset(subfam_summary_filt, subfam %in% subfam_grab)

## Identify Subfams with >= 5 proteins and pull those proteins
subfam_n5_grab <- subset(fam00018_subfam_summaries, count >= 5)$subfam
fam00018_genes_n5subfam <- subset(fam00018_genes_filt, subfam %in% subfam_n5_grab)

## Load in manual annotation of subfams with > 5 proteins (HHsearch vs. PFAM was used to get BCP domain number)
fam00018_subfam_manualanno <- fread("All_fam00018_n5_subfams/fam00018_subfam_manualAnno.txt")

## Merge manual annotations to proteins from subfams with >= 5 proteins 
fam00018_genes_n5subfam <- merge(fam00018_genes_n5subfam, fam00018_subfam_manualanno, by = "subfam")


## What is the distribution of domain architecture annotations (Detailed) across order level taxa - PLOT 1
fam00018_genes_AnnoCount <- data.table(fam00018_genes_n5subfam %>%
                                        group_by(genome, tax_Order, Anno) %>%
                                        summarise(count = n()))

fam00018_genes_AnnoCount <- subset(fam00018_genes_AnnoCount, Anno != "Unknown")

ggplot(fam00018_genes_AnnoCount, aes(x = tax_Order, y = count, fill = Anno)) +
  geom_jitter(alpha = 0.5, shape = 21, height = 0.2, width = 0.2, color = "black", size = 2) +
  facet_wrap(.~Anno) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


## What is the distribution of domain architecture annotations (Simple) across order level taxa - PLOT 2
fam00018_genes_SimpleAnnoCount <- data.table(fam00018_genes_n5subfam %>%
                                         group_by(genome, tax_Order, Simple_Anno, Func_Homlog) %>%
                                         summarise(count = n()))

fam00018_genes_SimpleAnnoCount <- subset(fam00018_genes_SimpleAnnoCount, Simple_Anno != "Unknown")

ggplot(fam00018_genes_SimpleAnnoCount, aes(x = tax_Order, y = count, fill = Func_Homlog)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5, shape = 21, height = 0.2, width = 0.2, color = "black", size = 2) +
  facet_wrap(.~Simple_Anno) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## What are the protien sizes and their functions by domain architecture (Simple) and across order level taxa - PLOT 3a
plot3 <- subset(fam00018_genes_n5subfam, Simple_Anno != "Unknown")

ggplot(plot3, aes(x = tax_Order, y = length, fill = Func_Homlog)) +
  geom_jitter(alpha = 0.5, shape = 21, height = 0.2, width = 0.2, color = "black", size = 2) +
  facet_wrap(.~Simple_Anno) +
  scale_y_log10() +
  xlab(NULL) +
  ylab("Protein Length (aa)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Frequency of each Simple Annotation across order level groups - PLOT 3b + PLOT 4
fam00018_genes_SimpleAnnoFreq <- data.table()

# loop through each annotation so we can get zero counts for genomes that do not have them
for (i in unique(fam00018_genes_SimpleAnnoCount$Simple_Anno)){
  
  tmp <- subset(fam00018_genes_SimpleAnnoCount, Simple_Anno == i)
  tmp <- tmp[,-2]
  tmp <- merge(genome_tax, tmp, by = "genome", all.x = T)
  tmp$Simple_Anno <- ifelse(is.na(tmp$Simple_Anno),i,tmp$Simple_Anno)
  tmp$count <- ifelse(is.na(tmp$count), 0 ,tmp$count)
  
  fam00018_genes_SimpleAnnoFreq <- rbind(tmp, fam00018_genes_SimpleAnnoFreq)
}

# add column for presence absence and remove order taxa with < 3 genomes
fam00018_genes_SimpleAnnoFreq$present <- ifelse(fam00018_genes_SimpleAnnoFreq$count > 0, 1, 0)
fam00018_genes_SimpleAnnoFreq <- subset(fam00018_genes_SimpleAnnoFreq, tax_Order_counts >= 3)

# Summarise counts with mean and SD and plot - (PLOT 3b)
plot_means <- groupwiseMean(count ~ Simple_Anno + tax_Order, data = fam00018_genes_SimpleAnnoFreq)

ggplot(plot_means, aes(x = reorder(tax_Order, -Mean), y = Mean, fill = Simple_Anno)) +
  geom_bar(stat = "identity") +
  facet_wrap(.~Simple_Anno) +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Summarise counts with mean and SD and plot using stacked bar - (PLOT 3c)
ggplot(plot_means, aes(x = Simple_Anno, y = Mean, fill = tax_Order)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab(NULL) +
  #facet_wrap()
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





##########$$$ REVIEW UPDATE 
##########$$$

# Summarize subtypes in Nitrosopsherales order
fam00018_genes_SimpleAnnoFreq_nitroso <- subset(fam00018_genes_SimpleAnnoFreq, tax_Order == "o__Nitrososphaerales")

ggplot(fam00018_genes_SimpleAnnoFreq_nitroso, aes(x = tax_Family, y = count, fill = Simple_Anno)) +
  #geom_boxplot() +
  #geom_jitter(alpha = 0.5, shape = 21, height = 0.0, width = 0.2, color = "black", size = 2) +
  geom_bar(position = "fill", stat = "identity") +
  #facet_wrap(.~Simple_Anno) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




# Summarise frequency of genomes that have a protein with that annotation - (PLOT 4)
plot_penetrance <- groupwiseMean(present ~ Simple_Anno + tax_Order + tax_Phylum, data = fam00018_genes_SimpleAnnoFreq)

ggplot(plot_penetrance, aes(x = reorder(tax_Order, -Mean), y = Mean, fill = Simple_Anno)) +
  geom_bar(stat = "identity") +
  facet_wrap(.~Simple_Anno) +
  xlab(NULL) +
  ylab("Fraction of Genomes with BCP Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Summarise frequency of genomes that have a protein with that annotation; Thermoprotea Only - (PLOT 4a)
plot_penetrance_thermoprot <- subset(plot_penetrance, tax_Phylum == "p__Thermoproteota" & Simple_Anno %in% c("1D_Cu_Small", "1D_Cu_Med", "2D_Cu"))
#plot_penetrance_thermoprot <- subset(plot_penetrance, tax_Phylum == "p__Thermoproteota")

## Using barplot with taxa as X axis 
ggplot(plot_penetrance_thermoprot, aes(x = reorder(tax_Order, -Mean), y = Mean, fill = Simple_Anno)) +
  geom_bar(stat = "identity") +
  facet_wrap(.~Simple_Anno, nrow = 3) +
  xlab(NULL) +
  ylab("Fraction of Genomes with BCP Type") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

#ggsave("~/Desktop/TMP_Fig3/Cu_Prot_Statistics/Fig3C_v1.pdf", width = 4, height = 6)

## Using dotplot with type as x axis
# ggplot(plot_penetrance_thermoprot, aes(x = Simple_Anno, y = Mean, fill = tax_Order)) +
#   geom_point(shape = 22, size = 4) +
#   xlab(NULL) +
#   ylab("Fraction of Genomes with BCP Type") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Summarise frequency of genomes that have a protein with that annotation; Thermoprotea Only - (PLOT 4b)
plot_penetrance_thermoplas <- subset(plot_penetrance, tax_Phylum == "p__Thermoplasmatota" & Simple_Anno %in% c("1D_Cu_Small", "1D_Cu_Med", "2D_Cu"))
plot_penetrance_thermoplas <- subset(plot_penetrance, tax_Phylum == "p__Thermoplasmatota")

ggplot(plot_penetrance_thermoplas, aes(x = reorder(tax_Order, -Mean), y = Mean, fill = Simple_Anno)) +
  geom_bar(stat = "identity") +
  facet_wrap(.~Simple_Anno,nrow = 3) +
  xlab(NULL) +
  ylab("Fraction of Genomes with BCP Type") +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

#ggsave("~/Desktop/TMP_Fig3/Cu_Prot_Statistics/Fig3B_v1.pdf", width = 4, height = 6)





# Try with all
# plot_penetrance_all <- subset(plot_penetrance, Simple_Anno %in% c("1D_Cu_Small", "1D_Cu_Med", "2D_Cu"))
# 
# ggplot(plot_penetrance_all, aes(x = reorder(tax_Order, -Mean), y = Mean, fill = Simple_Anno)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(.~Simple_Anno + tax_Phylum, nrow = 3, scales = "free") +
#   xlab(NULL) +
#   ylab("Fraction of Genomes with BCP Type") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         legend.position = "none")






## Plot counts per genome of only 2D and 3D Cu proteins in boxplot format to excentuate findings - PLOT 5
Simple_Anno_Select <- c("2D_Cu", "3D_Cu", "1D_Cu_Med")
plot_counts_2D3D <- subset(fam00018_genes_SimpleAnnoFreq, Simple_Anno %in% Simple_Anno_Select)

ggplot(plot_counts_2D3D, aes(x = reorder(tax_Order, -count), y = count, fill = tax_Phylum)) +
  geom_jitter(width = .1, size = 1, alpha = 0.8, height = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  xlab(NULL) +
  ylab("Protein Counts per Genome") +
  scale_fill_manual(values = c("steelblue4", "firebrick3")) +
  facet_wrap(.~Simple_Anno, ncol = 1) +
  #theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


## Do statistical testing to look for enrichment of proteins 
SimpleAnno_Stats <- list()

for (i in unique(fam00018_genes_SimpleAnnoFreq$Simple_Anno)){
  
  tmp <- subset(fam00018_genes_SimpleAnnoFreq, Simple_Anno == i) 
  tmp_krusk <- kruskal.test(count ~ tax_Order, data = tmp)
  
  if(tmp_krusk$p.value <= 0.05){
    SimpleAnno_Stats[[i]] <-pairwise.wilcox.test(tmp$count, tmp$tax_Order, p.adjust.method = "fdr")
  }
  
  else{SimpleAnno_Stats[[i]] <- "KW-Test Not Significant"}
}




### 04 ### Explanation:
### 04 ### 
### 04 ### 
### 04 ### 
### 04 ### 
### 04 ### 
### 04 ### 
### 04 ###
### 04 ### 
### 04 ###




### 05 ### We have now established that 1D_Cu_Med + 2D_Cu + 3D_Cu are statistically more prevelant in the AMO containing
### 05 ### species specificallty. Now we will do minimal analysis of their subfams and create a metadata file that will
### 05 ## allow us to pull the sequences associated with them.

## Get all genes + metadata for those with Simple_Anno == 1D_Cu_Med | 2D_Cu | 3D_Cu
Simple_Anno_Select <- c("2D_Cu", "3D_Cu", "1D_Cu_Med")
fam00018_genes_n5subfam_targeted <- subset(fam00018_genes_n5subfam, Simple_Anno %in% Simple_Anno_Select)

## Write out data for each class of BCP Individually
#fwrite(subset(fam00018_genes_n5subfam_targeted, Simple_Anno == "1D_Cu_Med"), "All_fam00018_Class_Analysis_and_Trees/1D_Cu_Med/fam00018_genes_n5subfam_1DCuMed.txt", quote = F, sep = "\t", na = NA)
#fwrite(subset(fam00018_genes_n5subfam_targeted, Simple_Anno == "2D_Cu"), "All_fam00018_Class_Analysis_and_Trees/2D_Cu/fam00018_genes_n5subfam_2DCu.txt", quote = F, sep = "\t", na = NA)
#fwrite(subset(fam00018_genes_n5subfam_targeted, Simple_Anno == "3D_Cu"), "All_fam00018_Class_Analysis_and_Trees/3D_Cu/fam00018_genes_n5subfam_3DCu.txt", quote = F, sep = "\t", na = NA)
#fwrite(subset(fam00018_genes_n5subfam_targeted, Simple_Anno == "2D_Cu" | Simple_Anno == "3D_Cu"),"All_fam00018_Class_Analysis_and_Trees/2D_3D_Cu/fam00018_genes_n5subfam_2D_3DCu.txt", quote = F, sep = "\t", na = NA)

## Get all subfam summaries for those associated with above genes 
subfam_grab <- unique(fam00018_genes_n5subfam_targeted$subfam)
fam00018_genes_n5subfam_targeted_summaries <- subset(fam00018_subfam_summaries, subfam %in% subfam_grab)

## Plot taxonomic composition of each subfam - PLOT 1
plot_BCP_target_subfamTax <- data.table(fam00018_genes_n5subfam_targeted %>%
                                          group_by(subfam, Simple_Anno, tax_Order) %>%
                                          summarise(count = n()))

ggplot(plot_BCP_target_subfamTax, aes(x = reorder(subfam, -count), y = count, fill = tax_Order)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(.~Simple_Anno, ncol = 1, scales = "free") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 





 

##### NOTE TO SELF GCA_004377185.1 HAS NO fam00018 proteins! 



### Get number of unique angelarcheales genomes with subfam
bar <- data.table(subfam = fam00018_subfam_domfreq2$subfam, Ang_Genomes = 0)


for (i in 1:nrow(bar)){
  tmp_sfm <- subset(master_anno_filt, subfam == bar$subfam[i] & order_tax == "RBG-16-68-12")
  
  if(nrow(tmp_sfm) == 0){
    next
  } 
  
  bar[i,2] <- length(unique(tmp_sfm$genome))
  
  
}

fwrite(bar, "~/Desktop/Ang_genomes_by_subfam.txt", sep ="\t")
