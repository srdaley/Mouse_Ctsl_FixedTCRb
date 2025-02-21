# Note: The senior author, Ludger Klein, prepared Figures 5c, d, e independently. 
# This code can be used to reproduce those results, 
# with the caveat that the percentage of reads produced here should differ slightly from those shown in the paper.
# This is because the read counts were normalised to "counts per million" (cpm) within each sample in Ludger Klein's dataset 
# whereas the dataset used in this code used the raw UMI counts (without normalisation).

# This code was used to determine the frequency of reads 
# corresponding to "recurrent" 'Ctsl-dependent' versus 'Ctsl-independent' TCRs 
# in WT CD4SP thymocyte samples from Fixed-β Ctsl+/+ mice (n = 3).  

# In RStudio, install and/or load "tidyverse", "stringr", and "scales" 
if (!require("tidyverse"))
  install.packages("tidyverse")
library(tidyverse)

if (!require("stringr"))
  install.packages("stringr")
library(stringr)

if (!require("scales"))
  install.packages("scales")
library(scales)

# Set working directory to the folder that contains data
data_folder <- "" # **Instruction: Make "data_folder" the path to the folder containing the file called "ctsl.tcra.clonotype.table.txt"
setwd(data_folder)
metadata <- read.delim(file = "ctsl.metadata.txt", header = T)
df <- read.delim(file = "ctsl.tcra.clonotype.table.txt", header = T)

# Define the "recurrent" TCRs
# Set "Genotypes of Interest (GOI)", "Cell Type of Interest (CTOI)" and "Tissue of Interest (TOI)"   
GOI <- c("WT", "CtslDeltaTEC")
CTOI <- "M2" 
TOI <- "Thymus"
###################################################################
dr2 <- filter(df, 
              Genotype %in% GOI,
              Cell_type == CTOI,
              Tissue == TOI) 

# Collapse to clonotype (unique combination of V, J and CDR3aa sequence).
dr3 <- dr2 %>%
  group_by(V, CDR3.amino.acid.sequence, J, cw) %>%
  summarise(count = sum(Count)) %>%
  rename(Sample = cw)

# Set order of factor(Sample)
SOI <- c("30903b-045", "30903b-046",  "30903b-047", 
         "30903b-048", "30903b-049", "30903b-050" )
dr3$Sample <- factor(dr3$Sample, levels = SOI)
dr3 <- arrange(dr3, Sample)
# Set factors and levels
WTSOI <- SOI[1:3] # First 3 samples are WT
KOSOI <- SOI[4:6] # Last 3 samples are CtslDeltaTEC

# Set "Categories of Interest (COI)"
COI <- c("wtpri", "shared", "kopri")
##########################################################################################################
# Make a list of TCR clonotypes detected exclusively in 3 or more samples.
# Display 1 TCR clonotype per row and 1 column per sample 
dr4 <- dr3 %>%
  pivot_wider(names_from = Sample,
              values_from = count,
              values_fill = 0)

# Convert counts to 0 (absent) or 1 (present)
dr5 <- dr4 
for (cn in 4:ncol(dr5)){
  dr5[,cn] <- as.numeric(as.logical(dr5[,cn] > 0))
}

# Filter TCR clonotypes detected in 3 or more samples
dr5$rs <- rowSums(dr5[,SOI])
dr6 <- dr5 %>% filter(rs>2) 
# We call those TCRs "recurrent TCRs".
# Categorise recurrent TCR clonotypes into "wt_private", "shared", or "ko_private" 
dr6$wtrs <- rowSums(dr6[,WTSOI])
dr6$kors <- rowSums(dr6[,KOSOI])
dr6$cat <- factor("shared", levels = COI)
dr6$cat[dr6$wtrs == 0] <- "kopri" 
dr6$cat[dr6$kors == 0] <- "wtpri" 

# Get columns necessary to identify and categorise "recurrent" TCRs
dr7 <- dr6 %>% select(V, CDR3.amino.acid.sequence, J, cat) %>% 
  filter(cat == "wtpri" | cat == "shared") # Exclude the "newcomer" TCRs

# Merge dr7 with df to get a complete list of detections of Ctsl-dependent and Ctsl-independent TCRs 
dr8 <- merge(dr7, df)

# Collapse to TCR + sample combinations
dr9 <- dr8 %>%
  group_by(V, CDR3.amino.acid.sequence, J, cw, cat) %>%
  summarise(count = sum(Count))

# Attach Cell_type column
all_rec <- metadata %>%
  distinct(cw, Cell_type) %>%
  merge(dr9) 
# all_rec contains all detections of each "recurrent" TCR
# across the entire dataset of 22 samples
# and includes a column indicating whether the TCR is 'Ctsl-dependent' (wtpri) or 'Ctsl-independent' (shared)
# Next focus on TCRs expressed by in pooled spleen and lymph node CD4+ T cells from Fixed-β Ctsl+/+ mice (n = 4).
# Set up character vectors to obtain data from samples of interest

unique(metadata$Cell_type)
CTOI <- "M2" 
met1 <- metadata %>%
  filter(Cell_type == CTOI,
         Genotype == "WT") %>%
  select(Cell_type, cw)
CWOI <- met1$cw

# Filter for CTOI
sd <- all_rec %>%
  filter(cw %in% CWOI) 

# Collapse to one clonotype per row and
# filter for clonotypes that were present in at least 3 samples
sd1 <- sd %>%
  group_by( Cell_type, cat, V, CDR3.amino.acid.sequence, J) %>%
  summarise(count = sum(count))  


# Summarise the percentage of clonotypes and reads
# that were 'Ctsl-dependent' (wtpri) or 'Ctsl-independent' (shared)
sd2 <- sd1 %>%
  group_by(cat) %>%
  summarise(
    clones = n(),
    reads = sum(count)
  ) %>%
  mutate(
    pc_clones = 100 * clones / nrow(sd1),
      pc_reads = 100 * reads / (sum(sd1$count))
  )
write.table(sd2, file = "Fig5e.thymic.txt",
            sep = "\t", quote = F, row.names = F)
