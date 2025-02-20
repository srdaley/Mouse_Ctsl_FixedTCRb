# This code was used to plot Trav and Traj gene segment usage
# in sub-repertoires of interest from the thymus
# as a function of their chromosomal distribution. 

# Sub-repertoires of interest from the thymus are:
# TCRs found in 3 WT samples but 0 out of 3 KO samples ("Ctsl-dependent") and
# TCRs found in 3 KO samples but 0 out of 3 WT samples ("Newcomer").

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
data_folder <- "" # **Instruction: Make "data_folder" the path to the folder containing the file called "fig3.clonotype.table.txt"
data_folder <- "/Users/daleys6/Documents/TCR_seq/ctsl/250103data/" 
# The preceding line is the path to the folder on the author's computer
setwd(data_folder)
df <- read.delim(file = "ctsl.tcra.clonotype.table.txt", header = T)

doa <- "250201" 
toa <- ".ctsl.cumsum.subrep"
destination_folder <- paste0("/Users/daleys6/Documents/TCR_seq/",
                             doa, toa,"/")
paste0(doa,toa)
# Create new folder with the result
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

# Collapse to clonotype (unique combination of V, J and CDR3aa sequence).
df3 <- df2 %>%
  group_by(V, CDR3.amino.acid.sequence, J, cw) %>%
  summarise(count = sum(Count)) %>%
  rename(Sample = cw)

# Set order of factor(Sample)
SOI <- c("30903b-045", "30903b-046",  "30903b-047", 
         "30903b-048", "30903b-049", "30903b-050" )
df3$Sample <- factor(df3$Sample, levels = SOI)
df3 <- arrange(df3, Sample)
# Set factors and levels
WTSOI <- SOI[1:3] # First 3 samples are WT
KOSOI <- SOI[4:6] # Last 3 samples are CtslDeltaTEC

# Set "Categories of Interest (COI)"
COI <- c("wtpri", "shared", "kopri")
##########################################################################################################
# Make a list of TCR clonotypes detected in 3 or more samples.
# TCRs detected only in Fixedβ Ctsl-delta-TEC samples will be red,
# TCRs detected only in Fixedβ littermate control samples will be blue, and
# TCRs detected in both groups will be gray. 
# Exclude TCR clonotypes that were detected in one sample only.

# Display 1 TCR clonotype per row and 1 column per sample 
dc2 <- df3 %>%
  pivot_wider(names_from = Sample,
              values_from = count,
              values_fill = 0)

# Convert counts to 0 (absent) or 1 (present)
dc3 <- dc2 
for (cn in 4:ncol(dc3)){
  dc3[,cn] <- as.numeric(as.logical(dc3[,cn] > 0))
}

# Filter TCR clonotypes detected in 2 or more mice
dc3$rs <- rowSums(dc3[,SOI])
dc4 <- dc3 %>% filter(rs>2)
# dc4 contains a list of TCR clonotypes that were detected in 3 or more samples.
# We call those TCRs "recurrent TCRs".
# Categorise TCR clonotypes into "wt_private", "shared", or "ko_private" 
dc4$wtrs <- rowSums(dc4[,WTSOI])
dc4$kors <- rowSums(dc4[,KOSOI])
dc4$cat <- factor("shared", levels = COI)
dc4$cat[dc4$wtrs == 0] <- "kopri" 
dc4$cat[dc4$kors == 0] <- "wtpri" 

# Get columns necessary to reattach the count data in df3
dc5 <- dc4 %>% select(V, CDR3.amino.acid.sequence, J, cat)

# dc5 contains the "recurrent" TCR clonotypes detected in 3 or more samples.
# Merge dc5 with df3, which contains the count data.
d3 <- merge(dc5, df3)

# Identify unique combinations of "cat" and "Sample" in d3 
d3$agg <- paste(d3$cat, d3$Sample, sep = "_")
unique(d3$agg)
# Here, we are not interested in the "shared" TCRs,
# so reduce the dataframe to those rows with "agg" strings that contain either "wtpri" or "kopri"
AOI <- c("wtpri_30903b-045", 
         "wtpri_30903b-046",  
         "wtpri_30903b-047", 
         "kopri_30903b-048",   
         "kopri_30903b-049",
         "kopri_30903b-050")
d4 <- d3 %>% filter(agg %in% AOI)
d4$agg <- factor(d4$agg, levels = AOI)

# As our file containing the Traj and Traj chromosomal locations
# does not have allele information,
# remove the allele information from the clonotype table.
d4$j1 <- substr(d4$J, 1, str_locate(d4$J, "[*]") - 1)
d4$j2 <- paste0("TRA", str_replace(d4$j1, "-",""))
j_names <- unique(d4$j2)

# Check that each row has 1 TRAJ call only 
max(str_count(j_names, "J")) # This command should return a result of 1

# Check TRAJ names match the TRAJ names in our mouse TRAV/TRAJ chromosomal location file
# If necessary, download the file called "imgt_mousevjlocus_mixcr.txt" from https://github.com/srdaley/CathepsinL_TCRseq  
locn <- read.delim(file = "imgt_mousevjlocus_mixcr.txt", header = TRUE, skip = 1)
#locn <- read.delim(file = "../r_codes/imgt_mousevjlocus_mixcr.txt", header = TRUE, skip = 1)
qj <- tibble(j_names,match(j_names, locn$IMGT_gene_name))
# qj should have a number for each value of "j_names"
# Rename "j2" column to "j"
d4 <- d4 %>% rename(j = j2)

# For Trav analysis, EXCLUDE rows with > 1 Trav call
d4$v_calls <- str_count(d4$V, ",") + 1
d5 <- d4 %>% filter(v_calls == 1)

# Add "TRA", remove first hyphen, and remove allele info from "V" column
d5$v1 <- str_replace(d5$V, "-", "")
d5$v2 <- paste0("TRA", substr(d5$v1, 1, str_locate(d5$v1, "[*]") - 1))
d5$v3 <- str_replace_all(d5$v2, "/", "-")
v_names <- unique(d5$v3)

# check TRAV names match our mouse TRAV/TRAJ chromosomal location file
qv <- tibble(v_names,match(v_names, locn$IMGT_gene_name))
# qv should have a number for each value of "j_names"
# Rename "v3" column to "v"
d5 <- d5 %>% rename(v = v3)
nrow(d4)
nrow(d5)
#############################################################################
# Use d5 with TRAV genes in the "v" column for TRAV analysis. 
# Use d4 with TRAJ genes in the "j" column for TRAJ analysis. 
#############################################################################
# For each unique combination of sample and sub-repertoire, 
# determine the number of unique TCR cDNA molecules or UMIs (count)
# that use each TRAV segment
# and make a dataframe with TRAV segments in rows 
# and samples in columns
df1_v <- d5 %>%
  group_by(agg, v) %>%
  summarise(Count = sum(count)) %>%
  pivot_wider(id_cols = v, names_from = agg, 
              values_from = Count, values_fill = 0 ) 

# Convert counts to per-sample frequency
for (i in 2:ncol(df1_v)) {
  df1_v[,i] <- df1_v[,i] / sum(df1_v[,i])}

# Make IMGT_gene_name a character vector
locn$IMGT_gene_name <- as.character(locn$IMGT_gene_name)

# Attach chromosomal location info for each TRAV segment
# and arrange rows from proximal (top) to distal (bottom)
v_posn <- merge(df1_v, locn, by.x = "v", by.y = "IMGT_gene_name") %>%
  arrange(desc(IMGT_gene_order))

# Write a source data file
v_o <- v_posn[, c(1,8,2:7)]
colnames(v_o)[1] <- "V element"
write.table(v_o, file = paste0(doa, toa, ".v.raw.txt"),
            row.names = F, quote = F, sep = "\t")

# Determine the cumulative frequency of UMIs that use each TRAV segment
for (cn in 2:(ncol(v_posn)-1)){
  v_posn[,cn] <- cumsum(v_posn[,cn])
}

# Make a long dataframe for plotting
vp <- v_posn %>%
  pivot_longer(names_to = "sample_subrep", values_to = "freq", 2:(ncol(v_posn)-1)) 
vp$subrep <- substr(vp$sample_subrep, 1, str_locate(vp$sample_subrep, "_")-1)


# Summarise by sub-repertoire and determine the mean and range
vp1 <- vp %>% 
  group_by(subrep,  IMGT_gene_order, v) %>%
  summarise(
    mean_freq = mean(freq),
    max_freq = max(freq),
    min_freq = min(freq)
  )

# For brevity of plot, remove "TRAV"
vp1$v_lab <- gsub( "TRAV","",vp1$v)

# Make x-axis run from proximal (left) to distal (right)  
xof <- vp1 %>%
  arrange(desc(IMGT_gene_order)) 

# Get x-axis labels
xo <- unique(xof$v_lab) 

# Set order of x-axis labels 
vp1$v_lab <- factor(vp1$v_lab, levels = xo)

unique(vp1$subrep)
pal <- c("wtpri" = "royalblue", 
         "kopri" = "firebrick1")

vp1$subrep <- factor(vp1$subrep, levels = COI)
# Plot it
pdf(file = paste (doa, toa,  "all.v.pdf", sep = "."), 
    height = 3.5, width = 9, useDingbats = F)
g <- ggplot(vp1, aes(x=v_lab, 
                     y=mean_freq, color = subrep, group = subrep)) + 
  geom_line () +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x =element_text(size = 8, angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20, hjust = .5, face = "italic"),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor.y  = element_blank(),
        legend.title = element_blank())+
  scale_x_discrete(labels=xo) +
  labs(y="Cumulative proportion of\nTCR cDNAs",
       title = "TRAV") +
  scale_colour_manual(values = pal) +
  geom_errorbar(aes(ymin = min_freq, ymax = max_freq))
g
dev.off()
# The plot above has labels.
pdf(file = paste (doa, toa,  "all.v.clean.pdf", sep = "."), 
    height = 3.5, width = 7, useDingbats = F)
g <- ggplot(vp1, aes(x=v_lab, 
                     y=mean_freq, color = subrep, group = subrep)) + 
  geom_line () +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x =element_text(size = 8, angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor.y  = element_blank(),
        legend.position = "none")+
  scale_x_discrete(labels=xo) +
  scale_colour_manual(values = pal) +
  geom_errorbar(aes(ymin = min_freq, ymax = max_freq))
g
dev.off()
# The plot above has no y-axis labels or legend.
##################################################################################
# Repeat the analysis for TRAJ segment usage.
# For each combination of sample and sub-repertoire, 
# determine the number of TCR cDNA molecules or UMIs (count)
# that use each TRAJ segment
# and make a dataframe with TRAJ segments in rows 
# and samples in columns
df1_j <- d4 %>%
  group_by(agg, j) %>%
  summarise(count = sum(count)) %>%
  pivot_wider(id_cols = j, names_from = agg,
              values_from = count, values_fill = 0 ) 

# Convert counts to per-sample frequency
for (i in 2:ncol(df1_j)) {
  df1_j[,i] <- df1_j[,i] / sum(df1_j[,i])}

# Attach chromosomal location info for each TRAJ segment
# and arrange rows from proximal (top) to distal (bottom)
j_posn <- merge(df1_j, locn, by.x = "j", by.y = "IMGT_gene_name") %>%
  arrange(IMGT_gene_order)

# Write a source data file
j_o <- j_posn[, c(1,8,2:7)]
colnames(j_o)[1] <- "J element"
write.table(j_o, file = paste0(doa, toa, ".j.raw.txt"),
            row.names = F, quote = F, sep = "\t")

# Determine the cumulative frequency of UMIs that use each TRAJ segment
for (cn in 2:(ncol(j_posn)-1)){
  j_posn[,cn] <- cumsum(j_posn[,cn])
}

# Make a long dataframe for plotting
jp <- j_posn %>%
  pivot_longer(names_to = "sample_subrep", values_to = "freq", 
               2:(ncol(j_posn)-1)) 
jp$subrep <- factor(substr(jp$sample_subrep, 1, str_locate(jp$sample_subrep, "_") -1), 
                    levels = COI)

# Summarise by sub-repertoire and determine the mean and range
jp1 <- jp %>% 
  group_by(subrep,  IMGT_gene_order, j) %>%
  summarise(
    mean_freq = mean(freq),
    max_freq = max(freq),
    min_freq = min(freq)
  )

# For brevity of plot, remove "TRAJ"
jp1$j_lab <- gsub( "TRAJ","",jp1$j)

# Make x-axis run from proximal (left) to distal (right)  
jof <- jp1 %>%
  arrange(IMGT_gene_order) 

# Get x-axis labels
jo <- unique(jof$j_lab) 

# Set order of x-axis labels 
jp1$j_lab <- factor(jp1$j_lab, levels = jo)
jp1$subrep <- factor(jp1$subrep, levels = COI)

# Plot it
pdf(file = paste (doa, toa, "all.j.pdf", sep = "."), 
    height = 3.5, width = 9, useDingbats = F)
g <- ggplot(jp1, aes(x=j_lab, 
                     y=mean_freq, color = subrep, group = subrep)) + 
  geom_line () +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x =element_text(size = 8, angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20, hjust = .5, face = "italic"),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor.y  = element_blank(),
        legend.title = element_blank())+
  scale_x_discrete(labels=jo) +
  labs(y="Cumulative proportion of\nTCR cDNAs",
       title = "TRAJ") +
  scale_colour_manual(values = pal) +
  geom_errorbar(aes(ymin = min_freq, ymax = max_freq))
g
dev.off()
# The above plot has labels.
pdf(file = paste (doa, toa, "all.j.clean.pdf", sep = "."), 
    height = 3.5, width = 7, useDingbats = F)
g <- ggplot(jp1, aes(x=j_lab, 
                     y=mean_freq, color = subrep, group = subrep)) + 
  geom_line () +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x =element_text(size = 8, angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=20, hjust = .5, face = "italic"),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor.y  = element_blank(),
        legend.position = "none")+
  scale_x_discrete(labels=jo) +
  scale_colour_manual(values = pal) +
  geom_errorbar(aes(ymin = min_freq, ymax = max_freq))
g
dev.off()
# The above plot has no y-axis labels or legend.
##################################################################################
# Get AUC for each combination of sample and sub-repertoire
va <- v_posn[, names(v_posn) %in% vp$sample_subrep]
# Get the sum of cumulative frequencies for each sample
varo <- tibble(names(va), colSums(va))
colnames(varo) <- c("agg", "sigv")

# To calculate AUC, divide the sum of cumulative frequencies 
# by the number of gene segments detected in the dataset
varo$v_auc <- varo$sigv/nrow(va)


##################################################################################
# Repeat for TRAJ 
ja <- j_posn[, names(j_posn) %in% jp$sample_subrep]
jaro <- tibble(names(ja), colSums(ja))
colnames(jaro) <- c("agg", "sigj")
jaro$j_auc <- jaro$sigj/nrow(ja)

#####################################################################################
# Merge TRAV and TRAJ results
out <- merge(varo,jaro, by = c("agg")) 
# Attach Subrepertoire information
SROI <- c("Ctsl-dependent", "Ctsl-independent", "Newcomer")
out$cat <- factor(substr(out$agg, 1, 5), levels = COI)
out$Subrepertoire <- SROI[match(out$cat, COI)]
out$Sample <- factor(substr(out$agg, 7, nchar(out$agg)), levels = SOI)
out <- out %>% arrange(Sample)
colnames(out)
out1 <- out[, c ("Subrepertoire", "Sample", "v_auc", "j_auc")]
out2 <- out1 %>%
  rename(AUCforTRAV = v_auc,
         AUCforTRAJ = j_auc)
# Make a source data file containing the area under the curve results.
write.table(out2, file = paste0(doa, toa, "auc.stats.txt"), sep = "\t", quote = F, row.names = F)
