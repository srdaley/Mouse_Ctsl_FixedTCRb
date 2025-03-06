# Morisita Horn analyses comparing TCR repertoires expressed by 
# CD5-hi and CD5-lo subsets of M2 CD4SP thymocytes from Ctsl+/+ mice
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
doa <- "250202" 
toa <- ".ctsl.mh.cd5_hi_lo"
destination_folder <- paste0("/Users/daleys6/Documents/TCR_seq/",
                             doa, toa,"/")
paste0(doa,toa)
#create new folder with the result 
setwd(destination_folder)

# Set "Cell Type of Interest (CTOI)" 
CTOI <- c("M2 CD5lo", "M2 CD5hi")

###################################################################
# use min_count to vary the minimum count threshold
min_count <- 1
df2 <- filter(df, 
              Cell_type %in% CTOI,
              Count >= min_count) 

# Collapse to clonotype (unique combination of V, J and CDR3aa sequence).
df3 <- df2 %>%
  group_by(V, CDR3.amino.acid.sequence, J, cw, Cell_type) %>%
  summarise(count = sum(Count)) %>%
  rename(Sample = cw)

# Set order of factor(Sample)
# By determining the CD5 status of each Sample
df3 %>% group_by(Cell_type, Sample) %>% summarise(clones = n())
# Odd Samples are CD5lo and Even Samples are CD5hi
unique(df3$Sample)
SOI <- c("31411-001", "31411-003", "31411-005", "31411-007",
         "31411-002", "31411-004", "31411-006", "31411-008")
df3$Sample <- factor(df3$Sample, levels = SOI)
df3 <- df3 %>% arrange(Sample)
# Make a dataframe with one row per clonotype
# and one column per Sample 
df4 <- df3 %>%
  pivot_wider(id_cols = c(V, CDR3.amino.acid.sequence, J),
              names_from = Sample,
              values_from = count,
              values_fill = 0)
# Convert data columns from integer vectors to numeric vectors
df4 <- df4 %>%
  mutate(across(where(is.integer), as.numeric))

x <- as.matrix(df4[,4:ncol(df4)])

results <- matrix(nrow = length(SOI), ncol = length(SOI))
for (a in 1:length(SOI)) {
  for(b in 1:length(SOI)){
    results[a,b] <- 1 - horn_morisita(x[,a], x[,b])
  }
}
colnames(results) <- SOI
rownames(results) <- SOI

# Write results to a Source Data file.
write.table(results, file = paste0(doa, toa, ".txt"), col.names = NA, sep = "\t", quote = F)

