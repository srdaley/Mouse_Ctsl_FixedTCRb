# Morisita Horn analyses on entire repertoires 
# of M2 CD4SP thymocytes
# are outputs of this script.

# In RStudio, install and/or load "tidyverse" and "abdiv". 
if (!require("tidyverse"))
  install.packages("tidyverse")
library(tidyverse)

if (!require("abdiv"))
  install.packages("abdiv")
library(abdiv)

# Set working directory to the folder that contains data
data_folder <- "" # **Instruction: Make "data_folder" the path to the folder containing the file called "fig3.clonotype.table.txt"
data_folder <- "/Users/daleys6/Documents/TCR_seq/ctsl/250103data/" 
# The preceding line is the path to the folder on the author's computer
setwd(data_folder)
df <- read.delim(file = "ctsl.tcra.clonotype.table.txt", header = T)

# For record keeping, set doa (date of analysis) and toa (type of analysis).
doa <- "250105" 
toa <- ".ctsl.mh"
destination_folder <- paste0("/Users/daleys6/Documents/TCR_seq/",
                             doa, toa,"/")
paste0(doa,toa)
#create new folder with the result 
setwd(destination_folder)

# Set "Genotypes of Interest (GOI)", "Cell Type of Interest (CTOI)" and "Tissue of Interest (TOI)"   
GOI <- c("WT", "CtslDeltaTEC")
CTOI <- "M2" 
TOI <- "Thymus"
###################################################################
# use min_count to vary the minimum count threshold
min_count <- 1
df2 <- filter(df, 
              Genotype %in% GOI,
              Cell_type == CTOI,
              Tissue == TOI,
              Count >= min_count) 

colnames(df2)
# Collapse to clonotype (unique combination of V, J and CDR3aa sequence).
df3 <- df2 %>%
  group_by(V, CDR3.amino.acid.sequence, J, cw, Genotype) %>%
  summarise(count = sum(Count)) %>%
  rename(Sample = cw)

# Set order of factor(Sample)
# By determining the Genotype of each Sample
df3 %>% group_by(Genotype, Sample) %>% summarise(clones = n())
# Samples 45-47 are WT and Samples 48-50 are CtslDeltaTEC
unique(df3$Sample)
SOI <- c("30903b-045", "30903b-046",  "30903b-047", 
         "30903b-048", "30903b-049", "30903b-050" )
df3$Sample <- factor(df3$Sample, levels = SOI)
df3 <- arrange(df3, Sample)
# Make a table with one row per clonotype and one column per sample
df4 <- df3 %>%
  pivot_wider(id_cols = c(V, CDR3.amino.acid.sequence, J),
              names_from = Sample,
              values_from = count,
              values_fill = 0)

# Convert data columns from integer vectors to numeric vectors.
# This enables larger numbers to be calculated. 
df4 <- df4 %>%
  mutate(across(where(is.integer), as.numeric))

x <- as.matrix(df4[,4:ncol(df4)])

# Use the horn_morisita command in abdiv 
res <- matrix(nrow = length(SOI), ncol = length(SOI))
for (a in 1:length(SOI)) {
  for(b in 1:length(SOI)){
    res[a,b] <- 1 - horn_morisita(x[,a], x[,b])
  }
}
colnames(res) <- SOI
rownames(res) <- SOI
# Write results to a Source Data file.
write.table(res, file = paste(doa, toa, "abdiv_mh_weighted.txt", sep = "."), 
            col.names = NA, sep = "\t", quote = F)
