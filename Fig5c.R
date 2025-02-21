# Note: The senior author, Ludger Klein, prepared Figures 5c, d, e independently. 
# This code can be used to reproduce those results.

# This code was used to identify "recurrent" TCRs,
# that are also "natural CD5lo" or "natural CD5hi" TCRs,
# defined as being present in 3 or more samples of the appropriate "natural CD5" subset (hi or low),
# and 0 samples of the opposite "natural CD5" subset.
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
metadata <- read.delim(file = "ctsl.metadata.txt", header = T)
df <- read.delim(file = "ctsl.tcra.clonotype.table.txt", header = T)

# Define the "recurrent" TCRs.
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
# How many "natural CD5lo" and "natural CD5hi" TCRs are present in all_rec?
# Set up character vectors to obtain data from samples of interest

unique(metadata$Cell_type)
CTOI <- c("M2 CD5lo","M2 CD5hi" )
met1 <- metadata %>%
  filter(Cell_type %in% CTOI) %>%
  select(Cell_type, cw)
CWOI <- met1$cw
CWLO <- CWOI[met1$Cell_type == CTOI[1]]
CWHI <- CWOI[met1$Cell_type == CTOI[2]]

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

# To filter "natural CD5lo/hi" TCRs 
# defined as being present in 3 or more samples of the appropriate "natural CD5" subset (hi or low),
# and 0 samples of the opposite "natural CD5" subset,
# convert counts to 0 (absent) or 1 (present)
for (cn in 5:ncol(sd1)){
  sd1[,cn] <- as.numeric(as.logical(sd1[,cn] > 0))
}

sd1$lors <- rowSums(sd1[,CWLO])
sd1$hirs <- rowSums(sd1[,CWHI])
# Deal with CD5lo first
sdlo <- sd1 %>% 
  filter(lors > 2, hirs == 0) 

# Delete the present/absent data and reattach the count data
lo1 <- sdlo %>%
  select(1:4) %>%
  merge(sd)

# Collapse to one clonotype per row
lo2 <- lo1 %>%
  group_by(Cell_type, cat, V, CDR3.amino.acid.sequence, J) %>%
  summarise(count = sum(count),
            incidences = n())
unique(lo2$incidences) # sanity check!

# Summarise the percentage of clonotypes and reads
# that were 'Ctsl-dependent' (wtpri) or 'Ctsl-independent' (shared)
CD5LO <- lo2 %>%
  group_by(cat) %>%
  summarise(
    clones = n(),
    reads = sum(count)
  ) %>%
  mutate(
    pc_clones = 100 * clones / nrow(lo2),
    pc_reads = 100 * reads / (sum(lo2$count))
  )
write.table(CD5LO, file = "Fig5c.cd5lo.txt",
            sep = "\t", quote = F, row.names = F)
# Repeat for CD5hi
sdhi <- sd1 %>% 
  filter(hirs > 2, lors == 0)

# Delete the present/absent data and reattach the count data
hi1 <- sdhi %>%
  select(1:4) %>%
  merge(sd)

# Collapse to one clonotype per row
hi2 <- hi1 %>%
  group_by(Cell_type, cat, V, CDR3.amino.acid.sequence, J) %>%
  summarise(count = sum(count),
            incidences = n())
unique(hi2$incidences) # sanity check!

# Summarise the percentage of clonotypes and reads
# that were 'Ctsl-dependent' (wtpri) or 'Ctsl-independent' (shared)
CD5HI <- hi2 %>%
  group_by(cat) %>%
  summarise(
    clones = n(),
    reads = sum(count)
  ) %>%
  mutate(
    pc_clones = 100 * clones / nrow(hi2),
    pc_reads = 100 * reads / (sum(hi2$count))
  )
write.table(CD5HI, file = "Fig5c.cd5hi.txt",
            sep = "\t", quote = F, row.names = F)
