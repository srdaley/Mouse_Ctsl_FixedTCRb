# Note: The senior author, Ludger Klein, prepared Figures 5c, d, e independently. 
# This code can be used to reproduce those results,
# with the caveat that the percentage of reads produced here should differ slightly from those shown in the paper.
# This is because the read counts were normalised to "counts per million" (cpm) within each sample in Ludger Klein's dataset 
# whereas the dataset used in this code used the raw UMI counts (without normalisation).

# This code was used to identify "recurrent" TCRs,
# that are also among the top ten most frequent ‘public’ clonotypes 
# (defined by presence in all four replicates of
# expanded LLO-Tet+ cells from Fixed-β Ctsl+/+ mice 7d after systemic immunization with LLO190-201),
# and plot the frequency of 'Ctsl-dependent' and 'Ctsl-independent' TCRs in a pie chart.

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
# Next focus on TCRs expressed by expanded LLO-Tet+ cells 
# from Fixed-β Ctsl+/+ mice (n = 4) 7d after systemic immunization with LLO190-201 
# Find the top ten most frequent ‘public’ clonotypes (defined by presence in all four replicates). 
# Set up character vectors to obtain data from samples of interest

unique(metadata$Cell_type)
CTOI <- "LLO-tetr+ CD4+ Tconv"
met1 <- metadata %>%
  filter(Cell_type == CTOI) %>%
  select(Cell_type, cw)
CWOI <- met1$cw

# Filter for CTOI
sd <- all_rec %>%
  filter(cw %in% CWOI) 
# Make a dataframe containing counts 
# with one row per TCR clonotype
# and one column per sample
sd1 <- sd %>%
  pivot_wider(id_cols = c("V", "CDR3.amino.acid.sequence", "J", "cat"),
              names_from = cw,
              values_from = count,
              values_fill = 0)

# To filter TCRs present in all 4  samples,
# convert counts to 0 (absent) or 1 (present)
for (cn in 5:ncol(sd1)){
  sd1[,cn] <- as.numeric(as.logical(sd1[,cn] > 0))
}

sd1$rs <- rowSums(sd1[,CWOI])

# Deal with CD5lo first
sd2 <- sd1 %>% 
  filter(rs == 4) 

# Delete the present/absent data and reattach the count data
sd3 <- sd2 %>%
  select(1:4) %>%
  merge(sd)

# Collapse to one clonotype per row
sd4 <- sd3 %>%
  group_by(Cell_type, cat, V, CDR3.amino.acid.sequence, J) %>%
  summarise(count = sum(count)) %>%
  arrange(desc(count))

# Get the top 10
sd5 <- head(sd4, 10) 
sd5$pc_reads = 100 * sd5$count / (sum(sd5$count))
write.table(sd5, file = "Fig5d.Top10.txt",
            sep = "\t", quote = F, row.names = F)
# Summarise the percentage of clonotypes and reads
# that were 'Ctsl-dependent' (wtpri) or 'Ctsl-independent' (shared)
sd6 <- sd5 %>%
  group_by(cat) %>%
  summarise(
    clones = n(),
    reads = sum(count)
  ) %>%
  mutate(
    pc_reads = 100 * reads / (sum(sd5$count))
  )
write.table(sd6, file = "Fig5d.Ctsl_dependency.txt",
            sep = "\t", quote = F, row.names = F)
