library(data.table)
library(dplyr)
library(ggplot2)
library(gggenes)
#library(superheat)

######## IMPORT AND SET UP ALL THE DATA ##########

'%notin%' <- Negate('%in%')

### Set WD
setwd("~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Analysis/20_9_14_Cu_Protein_Analysis/")

### Location of protein annotations
prot_seqs <- "~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Genome_Data/All_Proteins_fixed.faa"

### Import Gene Annotation Tables
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


### Get only annotations from Angelarcheales
AA_genomes <- subset(genome_tax, tax_Order == "o__RBG-16-68-12")
AA_anno_filt <- subset(master_anno_filt, genome %in% AA_genomes$genome)









####### FIRST WE GET THE AMO ABC CONTIG




####### Full amoABC Contig containing any of the members from all orgs

### Known Contigs containing Amo A B or C 
amoABC_contig <- c("RifSed_csp2_13ft_2_scaffold_3917_curated_", "14_0903_16_30cm_scaffold_19786_curated_", "14_0929_05_40cm_scaffold_25826_", "15con2h1.8.AGTTCC_scaffold_83_curated_",
                   "15con2h1.8.AGTTCC_scaffold_21785_curated_", "PLM3-1_170_b2_sep16_scaffold_385_curated_", "PLM3-1_200_b2_sep16_scaffold_7115_curated_", "RifSed_csp1_16ft_3_scaffold_52442_curated_",
                   "Sage1_R_100_16_scaffold_15587_curated_", "Sage1_R_20_16_scaffold_1211_curated_", "Sage2_L_100_16_scaffold_4332_curated_", "Sage2_R_115_16_scaffold_14192_curated_",
                   "14_0929_12_30cm_scaffold_18749_", "14_0929_05_40cm_scaffold_4006_curated_", "GCA_005878525.1_VBMB01000052.1_", "14_0927_05_20cm_scaffold_35140_",
                   "GCA_005878985.1_VBMA01000129.1_", "14_0927_12_40cm_scaffold_29139_", "14_0903_13_40cm_scaffold_10481_curated_", "15_D_Rain_20_5_09082015_scaffold_8_curated_")

## Run Grep
#for (i in 1:length(amoABC_contig)){
#  system(paste0("ggrep '",amoABC_contig[i],"' ",prot_seqs," | sed 's/>//' | awk '{print $1,$3,$5,$7}' >> Ang_fam00018_subfam_loci/loci_info_amoABC/amoABC_contigs.txt"))
#}

## Pull in All found Loci, rename columns, and add directionality for each gene
amoABC_loci <- fread("Ang_fam00018_subfam_loci/amoABC_contigs.txt")
colnames(amoABC_loci) <- c("ORF", "start", "end", "direction")
amoABC_loci$strand <- ifelse(amoABC_loci$direction == 1, "forward", "reverse")


### Attach metadata to Genes in Loci (This will remove species genome redundancy)
amoABC_loci <- merge(amoABC_loci, AA_anno_filt, by = "ORF")
amoABC_loci <- amoABC_loci[order(amoABC_loci$Order),]

### Fix Subfam/Fam NA annotations
amoABC_loci$subfam <- ifelse(is.na(amoABC_loci$subfam) == T,"Unk",amoABC_loci$subfam)
amoABC_loci$fam <- ifelse(is.na(amoABC_loci$fam) == T,"Unk",amoABC_loci$fam)

### Export for Manual Curation
fwrite(amoABC_loci,  "Ang_fam00018_subfam_loci/loci_info_amoABC/amoABC_loci.txt", sep = "\t")


#### AT THIS POINT YOU OPEN THE ABOVE FILE AND DO MANUAL ANNOTATION 
### 1) we remove all genes outside of 10 gene distance away from the amoCAXB operon
### 2) we select an orientation for amoCAXB and make the coordinates of contigs where this orientation is reversed negative so the plot can be drawn properly



## Plot Data
amoABC_input <- fread("Ang_fam00018_subfam_loci/loci_info_amoABC/amoABC_loci_plot_n10.txt")


# ggplot(amoABC_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno)) +
#   geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
#   facet_wrap(~ genome, scales = "free", ncol = 1) +
#   #scale_fill_manual(values = cols) +
#   theme_genes()


## Plot centered on amoB
dummies <- make_alignment_dummies(
  amoABC_input,
  aes(xmin = start, xmax = end, y = genome, id = Man_Anno),
  on = "amoB"
)

ggplot(amoABC_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno, label = Man_Anno)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_gene_label() +
  geom_blank(data = dummies) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  theme(legend.position = "none")

#ggsave("Ang_fam00018_subfam_loci/loci_info_amoABC/amoABC_loci_plot_n10.pdf", height = 10, width = 40)









####### NOW WE WILL LOOK FOR BCPs FLAKING CoxA 




### Get Loci for coxA - We will evaluate which BCPs are flaking coxA (not coxAC)

## First evaluate best way to find coxA 
foo <- subset(AA_anno_filt, KO == "K02274") # This is best 
foo <- subset(AA_anno_filt, fam == "fam00219")

## Pull annotation and contigs
coxA_anno <- subset(AA_anno_filt, KO == "K02274")
coxA_contig <- sub("_[0-9]*$","_",x = coxA_anno$ORF,perl = T)

## Run grep to pull all contig genes
for (i in 1:length(coxA_contig)){
  system(paste0("ggrep '",coxA_contig[i],"' ",prot_seqs," | sed 's/>//' | awk '{print $1,$3,$5,$7}' >> Ang_fam00018_subfam_loci/loci_info_coxA/coxA_contigs.txt"))
}

## Pull in All found Loci
coxA_loci <- fread("Ang_fam00018_subfam_loci/loci_info_coxA/coxA_contigs.txt")
colnames(coxA_loci) <- c("ORF", "start", "end", "direction")
coxA_loci$strand <- ifelse(coxA_loci$direction == 1, "forward", "reverse")

## Attach metadata/annotations to Genes in Loci
coxA_loci <- merge(coxA_loci, AA_anno_filt, by = "ORF")
coxA_loci <- coxA_loci[order(coxA_loci$Order),]

## Fix Subfam/Fam NA annotations
coxA_loci$subfam <- ifelse(is.na(coxA_loci$subfam) == T,"Unk",coxA_loci$subfam)
coxA_loci$fam <- ifelse(is.na(coxA_loci$fam) == T,"Unk",coxA_loci$fam)

## Export for Manual Curation
fwrite(coxA_loci,  "Ang_fam00018_subfam_loci/loci_info_coxA/coxA_loci.txt", sep = "\t")


#### AT THIS POINT YOU OPEN THE ABOVE FILE AND DO MANUAL ANNOTATION 
### 1) we remove all genes outside of 10 gene distance away from the amoCAXB operon
### 2) we select an orientation for amoCAXB and make the coordinates of contigs where this orientation is reversed negative so the plot can be drawn properly
### 3) we also resolve genomes that have more than one... we actually just give any genomes with 2 contigs different names (ie. genome a and b)


#### IMPORT MANUALLY CURRATED GENE LIST

## Gene list with 10 genes on either side of coxA locus 
coxA_input <- fread("Ang_fam00018_subfam_loci/loci_info_coxA/coxA_loci_plot_n10.txt")



## Plot centered on coxA
dummies <- make_alignment_dummies(
  coxA_input,
  aes(xmin = start, xmax = end, y = genome, id = Man_Anno),
  on = "CoxA"
)

ggplot(coxA_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno, label = Man_Anno)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_gene_label() +
  geom_blank(data = dummies) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  theme(legend.position = "none")

ggsave("Ang_fam00018_subfam_loci/loci_info_coxA/coxA_loci_plot_n10.pdf", height = 20, width = 40)








####### NOW WE WILL LOOK FOR BCPs FLAKING COMPLEX III 




### Get Loci for complex III - To identify complex III we will use the larger cytochrome B protein and test two different arCOGs
foo <- subset(AA_anno_filt, arCOG == "arCOG01721") # we use the shorter but better hitting cyt b subunit (this is the C-terminal portion)
foo <- subset(AA_anno_filt, arCOG == "arCOG04594")


### Get Loci for cplxIII
cplxIII_anno <- subset(AA_anno_filt, arCOG == "arCOG01721")
cplxIII_contig <- sub("_[0-9]*$","_",x = cplxIII_anno$ORF,perl = T)

## Run Grep
# for (i in 1:length(cplxIII_contig)){
#   system(paste0("ggrep '",cplxIII_contig[i],"' ",prot_seqs," | sed 's/>//' | awk '{print $1,$3,$5,$7}' >> Ang_fam00018_subfam_loci/loci_info_cplxIII/cplxIII_contigs.txt"))
# }

## Pull in All found Loci
cplxIII_loci <- fread("Ang_fam00018_subfam_loci/loci_info_cplxIII/cplxIII_contigs.txt")
colnames(cplxIII_loci) <- c("ORF", "start", "end", "direction")
cplxIII_loci$strand <- ifelse(cplxIII_loci$direction == 1, "forward", "reverse")


### Attach metadata to Genes in Loci
cplxIII_loci <- merge(cplxIII_loci, AA_anno_filt, by = "ORF")
cplxIII_loci <- cplxIII_loci[order(cplxIII_loci$Order),]

### Fix Subfam/Fam NA annotations
cplxIII_loci$subfam <- ifelse(is.na(cplxIII_loci$subfam) == T,"Unk",cplxIII_loci$subfam)
cplxIII_loci$fam <- ifelse(is.na(cplxIII_loci$fam) == T,"Unk",cplxIII_loci$fam)

### Export for Manual Curation
fwrite(cplxIII_loci,  "Ang_fam00018_subfam_loci/loci_info_cplxIII/cplxIII_loci.txt", sep = "\t")

#### AT THIS POINT YOU OPEN THE ABOVE FILE AND DO MANUAL ANNOTATION 
### 1) we remove all genes outside of 22 gene distance away from the cytB_N-term protein (need this weird length to capture full complex I)
### 2) we select an orientation for amoCAXB and make the coordinates of contigs where this orientation is reversed negative so the plot can be drawn properly
### 3) we also resolve genomes that have more than one... we actually just give any genomes with 2 contigs different names (ie. genome a and b)


#### IMPORT MANUALLY CURRATED GENE LIST

## Gene list with 10 genes on either side of cplxIII locus 
cplxIII_input <- fread("Ang_fam00018_subfam_loci/loci_info_cplxIII/cplxIII_loci_plot_n22.txt")



## Plot centered on cytB_N
dummies <- make_alignment_dummies(
  cplxIII_input,
  aes(xmin = start, xmax = end, y = genome, id = Man_Anno),
  on = "cytB_N"
)

ggplot(cplxIII_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno, label = Man_Anno)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_gene_label() +
  geom_blank(data = dummies) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  theme(legend.position = "none")

#ggsave("Ang_fam00018_subfam_loci/loci_info_cplxIII/cplxIII_loci_plot_n22.pdf", height = 20, width = 40)






####### NOW WE WILL PLOT 3 REPRESENTATIVE OPERONS TOGETHER FOR FIG 4

## Locus List
combined_input <- fread("Ang_fam00018_subfam_loci/loci_combiend_forfig/Combined_Operon_Plot.txt")

ggplot(combined_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno, label = Man_Anno, forward = direction)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_gene_label() +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() 


ggsave("Ang_fam00018_subfam_loci/loci_combiend_forfig/Loci_Figure4A_Raw.pdf", height = 7, width = 12)


























### Get Loci for subfam16590
subfam16590_anno <- subset(AA_anno_filt, subfam == "subfam16590")
subfam16590_contig <- sub("_[0-9]*$","_",x = subfam16590_anno$ORF,perl = T)

## Run Grep
#for (i in 1:length(subfam16590_contig)){
#  system(paste0("ggrep '",subfam16590_contig[i],"' ",prot_seqs," | sed 's/>//' | awk '{print $1,$3,$5,$7}' >> Ang_fam00018_subfam_loci/subfam16590_contigs.txt"))
#}

## Pull in All found Loci
subfam16590_loci <- fread("Ang_fam00018_subfam_loci/subfam16590_contigs.txt")
colnames(subfam16590_loci) <- c("ORF", "start", "end", "direction")
subfam16590_loci$strand <- ifelse(subfam16590_loci$direction == 1, "forward", "reverse")


### Attach metadata to Genes in Loci
subfam16590_loci <- merge(subfam16590_loci, AA_anno_filt, by = "ORF")
subfam16590_loci <- subfam16590_loci[order(subfam16590_loci$Order),]

### Fix Subfam/Fam NA annotations
subfam16590_loci$subfam <- ifelse(is.na(subfam16590_loci$subfam) == T,"Unk",subfam16590_loci$subfam)
subfam16590_loci$fam <- ifelse(is.na(subfam16590_loci$fam) == T,"Unk",subfam16590_loci$fam)

### Export for Manual Curation
#fwrite(subfam16590_loci,  "~/Desktop/subfam16590_loci.txt", sep = "\t")




## Test Plot
subfam16509_input <- fread("~/Desktop/TMP_Fig4/subfam16590_type1_locus.txt")


ggplot(subfam16509_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_manual(values = cols) +
  theme_genes()




dummies <- make_alignment_dummies(
  subfam16509_input,
  aes(xmin = start, xmax = end, y = genome, id = Man_Anno),
  on = "Pyc"
)

ggplot(subfam16509_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_blank(data = dummies) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes()

### Plot only contig from PLM3-1_170_b2_sep16_Maxbin2_081_curated

## Subset
subfam16509_input_subset <- subset(subfam16509_input, genome == "PLM3-1_170_b2_sep16_Maxbin2_081_curated")

## Plot
ggplot(subfam16509_input_subset, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno, forward = direction)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_manual(values = cols) +
  theme_genes()


ggsave("~/Desktop/test_genes2.pdf", width = 10, height = 6)





### Combined annotated plot of PLM3-1_170_b2_sep16_Maxbin2_081_curated Contigs
PLM_contigs <- fread("~/Desktop/TMP_AMO_Fig4/Combined_Operon_Plot.txt")

## Plot
ggplot(PLM_contigs, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno, forward = direction)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_manual(values = cols) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

ggsave("~/Desktop/combined_operons.pdf", width = 12, height = 6)




















### Get Loci for coxAC 
coxAC_anno <- subset(AA_anno_filt, fam == "fam00219" & length >= 600)
coxAC_contig <- sub("_[0-9]*$","_",x = coxAC_anno$ORF,perl = T)

## Run Grep
for (i in 1:length(coxAC_contig)){
  system(paste0("ggrep '",coxAC_contig[i],"' ",prot_seqs," | sed 's/>//' | awk '{print $1,$3,$5,$7}' >> Ang_fam00018_subfam_loci/coxAC_contigs.txt"))
}

## Pull in All found Loci
coxAC_loci <- fread("Ang_fam00018_subfam_loci/coxAC_contigs.txt")
colnames(coxAC_loci) <- c("ORF", "start", "end", "direction")
coxAC_loci$strand <- ifelse(coxAC_loci$direction == 1, "forward", "reverse")


### Attach metadata to Genes in Loci
coxAC_loci <- merge(coxAC_loci, AA_anno_filt, by = "ORF")
coxAC_loci <- coxAC_loci[order(coxAC_loci$Order),]

### Fix Subfam/Fam NA annotations
coxAC_loci$subfam <- ifelse(is.na(coxAC_loci$subfam) == T,"Unk",coxAC_loci$subfam)
coxAC_loci$fam <- ifelse(is.na(coxAC_loci$fam) == T,"Unk",coxAC_loci$fam)

### Export for Manual Curation
fwrite(coxAC_loci,  "~/Desktop/coxAC_loci.txt", sep = "\t")




## Test Plot
coxAC_input <- fread("~/Desktop/TMP_AMO_Fig4/coxAC_loci_plot.txt")


ggplot(coxAC_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_manual(values = cols) +
  theme_genes()


## Plot centered on coxAC
dummies <- make_alignment_dummies(
  coxAC_input,
  aes(xmin = start, xmax = end, y = genome, id = Man_Anno),
  on = "CoxAC"
)

ggplot(coxAC_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno, label = Man_Anno)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_gene_label() +
  geom_blank(data = dummies) +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes()








# ### Get Loci for amoC
# fam00948_anno <- subset(AA_anno_filt, fam == "fam00948")
# fam00948_contig <- sub("_[0-9]*$","_",x = fam00948_anno$ORF,perl = T)
# 
# ## Run Grep
# for (i in 1:length(fam00948_contig)){
#   system(paste0("ggrep '",fam00948_contig[i],"' ",prot_seqs," | sed 's/>//' | awk '{print $1,$3,$5,$7}' >> Ang_fam00018_subfam_loci/fam00948_contigs.txt"))
# }
# 
# ## Pull in All found Loci
# fam00948_loci <- fread("Ang_fam00018_subfam_loci/fam00948_contigs.txt")
# colnames(fam00948_loci) <- c("ORF", "start", "end", "direction")
# fam00948_loci$strand <- ifelse(fam00948_loci$direction == 1, "forward", "reverse")
# 
# 
# ### Attach metadata to Genes in Loci
# fam00948_loci <- merge(fam00948_loci, AA_anno_filt, by = "ORF")
# fam00948_loci <- fam00948_loci[order(fam00948_loci$Order),]
# 
# ### Fix Subfam/Fam NA annotations
# fam00948_loci$subfam <- ifelse(is.na(fam00948_loci$subfam) == T,"Unk",fam00948_loci$subfam)
# fam00948_loci$fam <- ifelse(is.na(fam00948_loci$fam) == T,"Unk",fam00948_loci$fam)
# 
# ### Export for Manual Curation
# fwrite(fam00948_loci,  "~/Desktop/fam00948_loci.txt", sep = "\t")
# 
# 
# 
# 
# ## Test Plot
# fam00948_input <- fread("~/Desktop/fam00948_type1_locus.txt")
# 
# 
# ggplot(fam00948_input, aes(xmin = start, xmax = end, y = genome, fill = Man_Anno)) +
#   geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
#   facet_wrap(~ genome, scales = "free", ncol = 1) +
#   #scale_fill_manual(values = cols) +
#   theme_genes()





















#### THIS NEEDS TO MOVE TO Cu PROTEIN ANALYSIS 

cheese <- subset(AA_anno_filt, fam =="fam00018")
nex <- plyr::count(cheese, vars = c("genome", "subfam"))
jim <- dcast(nex, subfam ~ genome)
jim[is.na(jim)] <- 0

jim_mat <- as.matrix(jim[,2:ncol(jim)], rownames = jim$subfam)


superheat(jim_mat,
          pretty.order.cols = T,
          pretty.order.rows = T,
          col.dendrogram = T,
          dist.method = "euc",
          linkage.method = "ward.D2",
          scale = T,
          left.label.size = 0.1,
          left.label.text.size = 3,heat.lim = c(0,10),
          #force.bottom.label = TRUE,
          bottom.label.text.size = 4,
          bottom.label.text.angle = 90,
          bottom.label.size = 0.2)

superheat(t(jim_mat),
          pretty.order.cols = T,
          pretty.order.rows = T,
          row.dendrogram = T,
          col.dendrogram = T,
          dist.method = "euc",
          linkage.method = "ward.D2",
          scale = F,
          left.label.size = 0.3,
          left.label.text.size = 3,heat.lim = c(0,10),
          #force.bottom.label = TRUE,
          bottom.label.text.size = 4,
          bottom.label.text.angle = 90,
          bottom.label.size = 0.2)




### Get Loci for subfam16590
subfam17112_anno <- subset(AA_anno_filt, subfam == "subfam17112")
subfam17112_contig <- sub("_[0-9]*$","_",x = subfam17112_anno$ORF,perl = T)

## Run Grep
#for (i in 1:length(subfam17112_contig)){
#  system(paste0("ggrep '",subfam17112_contig[i],"' ",prot_seqs," | sed 's/>//' | awk '{print $1,$3,$5,$7}' >> Ang_fam00018_subfam_loci/subfam17112_contigs.txt"))
#}

## Pull in All found Loci
subfam17112_loci <- fread("Ang_fam00018_subfam_loci/subfam17112_contigs.txt")
colnames(subfam17112_loci) <- c("ORF", "start", "end", "direction")
subfam17112_loci$strand <- ifelse(subfam17112_loci$direction == 1, "forward", "reverse")


### Attach metadata to Genes in Loci
subfam17112_loci <- merge(subfam17112_loci, AA_anno_filt, by = "ORF")
subfam17112_loci <- subfam17112_loci[order(subfam17112_loci$Order),]

### Fix Subfam/Fam NA annotations
subfam17112_loci$subfam <- ifelse(is.na(subfam17112_loci$subfam) == T,"Unk",subfam17112_loci$subfam)
subfam17112_loci$fam <- ifelse(is.na(subfam17112_loci$fam) == T,"Unk",subfam17112_loci$fam)

### Export for Manual Curation
#fwrite(subfam16590_loci,  "~/Desktop/subfam16590_loci.txt", sep = "\t")