library(data.table)
library(ggplot2)
library(propr)


########### Import Needed Data

### Set Working Directory
setwd("~/Dropbox/Banfield_Lab_Files/Projects/AMO_Archaea/Analysis/20_9_25_Cross_Biome_Mapping/Abundance_Analysis/")


### Non-Redundant Taxonomy List
tax_template <- fread("Tax_Template.txt", sep = "\t")


### Abundance Data
Abundance_List <- list()
Abundance_List[[1]] <- fread("Mapping_Data_Raw/Master_Table_Angelo.txt")
Abundance_List[[2]] <- fread("Mapping_Data_Raw/Master_Table_East_River.txt")
Abundance_List[[3]] <- fread("Mapping_Data_Raw/Master_Table_Ocean.txt")
Abundance_List[[4]] <- fread("Mapping_Data_Raw/Master_Table_Rifle.txt")
Abundance_List[[5]] <- fread("Mapping_Data_Raw/Master_Table_Riv_Sage.txt")


### Metadata
metadata <- fread("Master_Metadata.txt")



############ Aggregate All Data and Create Structures for Subsequent Analysis

### Aggregate All Count Data
Master_Abund <- tax_template

for (i in 1:length(Abundance_List)){
  
  Master_Abund <- merge(Master_Abund, Abundance_List[[i]], by = "ConsensusLineage", all.x = T)
  
}

### Replace NA values with Zeros
Master_Abund[is.na(Master_Abund)] <- 0


### Aggregate Counts by Taxonomy at Order Level
Agg_Abund <- aggregate(. ~ Domain + Phylum + Class + Order, Master_Abund[,-c(1,6,7)], sum)


### Create Compound Naming and then create transposed matrix with proper row and column names
Agg_Abund$Compound_Name <-  paste(Agg_Abund$Domain,Agg_Abund$Phylum,Agg_Abund$Class,Agg_Abund$Order, sep = ";")
rownames(Agg_Abund) <- Agg_Abund$Compound_Name

Agg_Abund_mat <- t(Agg_Abund[,-c(1:4,ncol(Agg_Abund))])

# Save abundance table used for analysis 
#fwrite(as.data.frame(t(Agg_Abund_mat)), "~/Desktop/TMP_Fig2/Abundance/Aggregated_Abundance_Table.txt", sep = "\t", row.names = T)


############# Perform Plotting of Read Fractions Associated with Angelarcheales and Nitrosospheria

### Build Data Frame of Values
Angel_Frac <- Agg_Abund[which(Agg_Abund$Order == "o__RBG-16-68-12"),5:(ncol(Agg_Abund)-1)]/colSums(Agg_Abund[,5:(ncol(Agg_Abund)-1)])
Angel_Name <- rownames(Angel_Frac)
mean(as.numeric(Angel_Frac[1,]))
sd(as.numeric(Angel_Frac[1,]))

Nitroso_Frac <- Agg_Abund[which(Agg_Abund$Order == "o__Nitrososphaerales"),5:(ncol(Agg_Abund)-1)]/colSums(Agg_Abund[,5:(ncol(Agg_Abund)-1)])
Nitroso_Name <- rownames(Nitroso_Frac)
mean(as.numeric(Nitroso_Frac[1,]))
sd(as.numeric(Nitroso_Frac[1,]))

L6_Fractions <- rbind(Angel_Frac, Nitroso_Frac)
L6_Fractions$Group <- c("Angelarcheales", "Nitrososphaerales")


L6_Fractions <- melt(L6_Fractions)
colnames(L6_Fractions) <- c("Group", "Map_Sample", "Fraction")

### Merge in Metadata
L6_Fractions <- merge(L6_Fractions, metadata, by = "Map_Sample")


### Plot raw abundances
ggplot(L6_Fractions, aes(x = Sample_Name, y = Fraction, fill = Group)) +
  geom_bar(stat = "identity") +
  #scale_y_sqrt() +
  xlab(NULL) +
  ylab("Fraction of L6 Reads") +
  scale_fill_manual(values = c("firebrick3", "steelblue")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1),
        panel.grid.major.x = element_blank()) 

#ggsave("~/Desktop/TMP_Fig2/Abundance/L6_Abundance_Plot.pdf", width = 14, height = 4)


############# Perform Correlation analysis of Angel + Nitroso
## Total samples == 185

### Subset to taxa that show at least 5 counts in 33% of samples
keep <- apply(Agg_Abund_mat, 2, function(x) sum(x >= 5) >= 61) ## At least 5 counts in 33% of samples (61 samples)
keep_mat <-  Agg_Abund_mat[,keep]


# test_mat <- as.matrix(data.frame(A = c(0,0,1,0,0),
#                                  B = c(0,0,0,0,0),
#                                  C = c(0,0,0,0,0)))
# test_rho <- propr(test_mat, metric = "rho")
# test_rho@results


### Plot this distribution of what fraction of summed reads were lost for each sample after filtering 
mean(rowSums(keep_mat) / rowSums(Agg_Abund_mat))
sd(rowSums(keep_mat) / rowSums(Agg_Abund_mat))
boxplot(rowSums(keep_mat) / rowSums(Agg_Abund_mat))



### Calculate proportionality metric rho
### propr expects samples in rows and "genes" or features in columns
set.seed(123)
L6_rho <- propr(keep_mat, metric = "rho", p = 1000)

### Identify max and min rho to 
max_rho <- max(abs(L6_rho@results$propr))
min_rho <- min(abs(L6_rho@results$propr))

### Identify FDR cutoff < 0.001
updateCutoffs(L6_rho, cutoff = seq(from = 0.05, to = 0.40,by = 0.01), ncores = 12)
## cutoff is rho > 0.25 for FDR <0.001


### Output results table - Unfiltered
L6_rho_res <- L6_rho@results
hist(L6_rho_res$propr, breaks = 50)
L6_rho_res$Partner <- colnames(keep_mat)[L6_rho_res$Partner]
L6_rho_res$Pair <- colnames(keep_mat)[L6_rho_res$Pair]
fwrite(L6_rho_res, "~/Desktop/TMP_Fig2/Abundance/L6_Proportionality.txt", sep = "\t")

### Output results table - filtered
L6_rho_res_filt <- subset(L6_rho_res, propr >= 0.25)
fwrite(L6_rho_res_filt, "~/Desktop/TMP_Fig2/Abundance/L6_Proportionality_filt.txt", sep = "\t")

### Do some individual plotting if desired 
numb <- 1755 #1755 - This is the comparison between angelarcheales and Nitrosospherales # row number in L6_rho_res table 

## Plotting of raw values
qplot(x = keep_mat[,L6_rho@results$Partner[numb]],
      y = keep_mat[,L6_rho@results$Pair[numb]], log = "xy") +
  geom_smooth()

## Plotting of clr normalized values - Nitrosopherales vs Angelarcheales 
nice_plot <- data.frame(x = L6_rho@logratio[,L6_rho@results$Partner[numb]],
                        y = L6_rho@logratio[,L6_rho@results$Pair[numb]])

ggplot(nice_plot, aes (x = x , y = y)) +
  geom_point(shape = 21, fill = "firebrick4", alpha = 0.8) +
  geom_smooth(method = "lm", linetype = 2) +
  theme_bw() +
  xlab("Normalized Ca. Angelarcheales Abundance") +
  ylab("Normalized Nitrosospherales Abundance")


cor.test(L6_rho@logratio[,L6_rho@results$Partner[numb]],L6_rho@logratio[,L6_rho@results$Pair[numb]], method = "spear")





### Additional Plots for Supplementary Material

## Plotting of clr normalized values - Angelarcheales vs Nitrospirales
numb <- 1756 #This is the comparison between Angelarcheales and Nitrospirales # row number in L6_rho_res table 

nice_plot <- data.frame(x = L6_rho@logratio[,L6_rho@results$Partner[numb]],
                        y = L6_rho@logratio[,L6_rho@results$Pair[numb]])

ggplot(nice_plot, aes (x = x , y = y)) +
  geom_point(shape = 21, fill = "steelblue", alpha = 0.8) +
  geom_smooth(method = "lm", linetype = 2) +
  theme_bw() +
  xlab("Normalized Ca. Angelarcheales Abundance") +
  ylab("Normalized Nitrospirales Abundance")


## Plotting of clr normalized values - Nitrosospherales vs Nitrospirales
numb <- 990 #This is the comparison between Nitrosospherales vs Nitrospirales # row number in L6_rho_res table 

nice_plot <- data.frame(x = L6_rho@logratio[,L6_rho@results$Partner[numb]],
                        y = L6_rho@logratio[,L6_rho@results$Pair[numb]])

ggplot(nice_plot, aes (x = x , y = y)) +
  geom_point(shape = 21, fill = "orange", alpha = 0.8) +
  geom_smooth(method = "lm", linetype = 2) +
  theme_bw() +
  xlab("Normalized Nitrospirales Abundance") +
  ylab("Normalized Nitrosospherales Abundance") + 
  coord_flip()


