# This code was used to analyse the V-J junction.
# in sub-repertoires of interest from the thymus.
# Sub-repertoires of interest from the thymus are:
# TCRs found in 3 WT samples and 0 KO samples ("Ctsl-dependent") and
# TCRs found in 3 KO samples and 0 WT samples ("Newcomer").

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
df <- read.delim(file = "ctsl.tcra.clonotype.table.txt", header = T) %>%
  rename(Sample = cw)

doa <- "250201" 
toa <- ".ctsl.vjj.subrep"
destination_folder <- paste0("/Users/daleys6/Documents/TCR_seq/",
                             doa, toa,"/")
paste0(doa,toa)
# Create new folder with the result
setwd(destination_folder)

# Plan:
# Get a list of unique CDR3nt sequences that encode "recurrent TCRs" that are "private". 
# This means TCR clonotypes detected in 2 or more WT samples only (Ctsl-dependent) OR
# detected in 2 KO samples only ("Newcomer TCRs"). 
# Identify all unique TCR clones, i.e. unique combinations of V, J, and CDR3nucleotide
# Analyse those CDR3nt sequences using the IMGT JunctionAnalysis tool in batches of 5000  (lines x - y of this code).
# Versions: IMGT/JunctionAnalysis program version: 2.3.1 (2 December 2024)  
# IMGT/JunctionAnalysis reference directory release: 202449-1 (2 December 2024)
# Options on the https://www.imgt.org/IMGT_jcta/analysis page:
# Species: Mus musculus (strain C57BL6/J)
# Locus: TRA
# List of all eligible D-GENE: No
# Colored IMGT AA classes and histogram: No
# Ouput (sic) order: Same order as input
# 5' and 3' ends of the JUNCTION: May start and/or end with any codon
# Number of accepted mutations: 1  in 3'V-REGION and 1 in 5'J-REGION
# Delimitation of 3'V-REGION, D-REGION and 5'J-REGION: Stop trimming with the first encountered identical nucleotide
# Otherwise, the default options were used.

# wherein "Vda" shows the number of nucleotides deleted from (-) or P nucleotides added to (+) the TRAV gene, 
# "Nnt" shows the number of N nucleotides present, and
# Jda shows the number of nucleotides deleted from (-) or P nucleotides added to (+) the TRAJ gene 
##########################################################################################################

# Make a list of TCR clonotypes detected exclusively in 2 or more Fixedβ Ctsl-delta-TEC samples (red)
# or 2 or more Fixedβ littermate control (blue) samples. 
# Exclude TCR clonotypes that were detected in one sample only.
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
  group_by(V, CDR3.amino.acid.sequence, J, Sample) %>%
  summarise(count = sum(Count)) 

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
# Make a list of TCR clonotypes detected exclusively in 3 or more samples.
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

# Filter TCR clonotypes detected in 3 or more mice
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

# Get columns necessary to reattach the count data 
dc5 <- dc4 %>% select(V, CDR3.amino.acid.sequence, J, cat)

# dc5 contains the "recurrent" TCR clonotypes detected in 3 or more samples.
# Merge dc5 with df2, which contains the count data for each clone.
d3 <- merge(dc5, df2)

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

# Exclude TCR clonotypes with > 1 TRAV or > 1 TRAJ call
d4$v_calls <- str_count(d4$V, ",") + 1
d4$j_calls <- str_count(d4$J, ",") + 1
d5 <- d4 %>% filter(v_calls == 1, j_calls == 1)
nrow(d4) - nrow(d5)


# Get desired columns from d5 
d6 <- d5[ ,c("Sample", "agg", "V", "J", "CDR3.nucleotide.sequence", "Count")]
# d6 contains a row for each instance of unique TCR clones
# that encode a TCR clonotype detected in 3 or more samples.

# Now get unique combinations of "V", "J", and "CDR3.nucleotide.sequence"
d7 <- d6 %>%
  group_by(V, J, CDR3.nucleotide.sequence) %>%
  summarise(clones = n()) %>%
  arrange(desc(clones))

# Duplicate clones occupy two rows in the clonotype table,
# but they have identical nucleotide sequences. Why?

# Duplicate clones arise in the original dataset
# due to the presence and absence of the sequence alignment with TRAC
# resulting in 3 outcomes for the column called C: c("C*01", "", "C*01,C*01").

# The column called C is absent from the clonotype table that we have used as the starting point for all our analyses. 

# This has no effect on analyses based on the CDR3 amino acid sequences used to distinguish TCR clonotypes.

# For analyses based on CDR3 nucleotide sequences, such as the V-J junction analyses, 
# duplicate clones should be collapsed into a single clone,
# with a count that equals the sum of the counts of each member of the duplicate clone.

# Format the TCR clonotypes in d7 for input to IMGT Junction Analysis tool.
# Make a column for the sequence identifier
d7$seq_no <- paste0(">seq", rownames(d7))

# Add characters to match input format required by the IMGT JunctionAnalysis tool.
d7$TRAV <- paste0("TRAV", substr(d7$V,3, nchar(d7$V)))
d7$TRAJ <- paste0("TRAJ", substr(d7$J,3, nchar(d7$J)))
d7$seq_id <- paste(d7$seq_no, d7$TRAV, d7$TRAJ, sep = ", ")

# Convert to FASTA format
X <- d7[,c("seq_id", "CDR3.nucleotide.sequence")]
D <- do.call(rbind, lapply(seq(nrow(X)), function(i) t(X[i, ])))

# The IMGT JunctionAnalysis tool handles a maximum of 5000 sequences.
# So we need to write 6 files containing up to 5000 sequences each.
write.table(head(D,1e4), file = paste(doa, toa, "00001_05000.txt", sep = "."), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(D[10001:20000,], file = paste(doa, toa, "05001_10000.txt", sep = "."), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(D[20001:30000,], file = paste(doa, toa, "10001_15000.txt", sep = "."), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(D[30001:nrow(D),], file = paste(doa, toa, "15001_15128.txt", sep = "."), row.names = FALSE, col.names = FALSE, quote = FALSE)



# Read the results of IMGT JunctionAnalysis back in.
res1 <- read.delim(file = "imgt_ja_00001_05000.txt", header = F)
res2 <- read.delim(file = "imgt_ja_05001_10000.txt", header = F)
res3 <- read.delim(file = "imgt_ja_10001_15000.txt", header = F)
res4 <- read.delim(file = "imgt_ja_15001_15128.txt", header = F)

my.list <- list(res1, res2, res3, res4)
res <- bind_rows(my.list)
colnames(res) <- c("null", "seq", "Vname", "V3", "P1", "N", "P2", "J5", "Jname", "Vmut", "Jmut", "Ngc", "decrypt")
n_distinct(res$seq)
# IMGT Junction Analysis tool did not return a result for 34 sequences.

# As "res" does not contain the original CDR3nt sequence,
# we can use "seq" to attach the junctional statistics to each CDR3nt sequence in d7
d7$seq <- substr(d7$seq_no, 2, nchar(d7$seq_no))
d8 <- merge(res, d7, by = "seq")

# Keep a record of these results.
# The datafile is called "250107..ctsl.vjj.subrep.imgt.ja.fullresults.txt"
# and it is available for download from https://github.com/srdaley/CathepsinL_TCRseq.
 write.table(d8, file = paste(doa, toa, "imgt.ja.fullresults.txt", sep = "."), row.names = F, quote = F, sep = "\t")

# Convert the output from IMGT Junction Analysis into a format that shows
# the estimated number of nucleotides added or deleted at the V-J junction
# Select only the necessary columns 
colnames(d8)
d9 <- d8[, c("V", "J", "CDR3.nucleotide.sequence", "decrypt")]

# We will use the "decrypt" column to get the results of interest.
# Remove first and last parentheses
d9$dec0 <- substr(d9$decrypt, 2, nchar(d9$decrypt)-1)

# Replace all parentheses with underscores
d9$dec1 <- str_replace_all(d9$dec0, "[)]", "_")
d9$dec2 <- str_replace_all(d9$dec1, "[{]", "_")
d9$dec3 <- str_replace_all(d9$dec2, "[}]", "_")
d9$dec4 <- str_replace_all(d9$dec3, "[()]", "_")

junc_stats <- d9 %>% 
  separate(dec4, into = c("Vgl", "Vda", "Nnt", "Jda", "Jgl"), sep = "_") %>%
  # Select only the necessary columns
  select("V", "J", "CDR3.nucleotide.sequence" ,"Vda", "Nnt", "Jda") 
# wherein "Vda" shows the number of nucleotides deleted from (-) or P nucleotides added to (+) the TRAV gene, 
# "Nnt" shows the number of N nucleotides present, and
# Jda shows the number of nucleotides deleted from (-) or P nucleotides added to (+) the TRAJ gene 
###################################################################################################################
# junc_stats contains the junctional statistics for each CDR3nt sequence along with V and J info ##################
###################################################################################################################
# Merge the junctional analysis results (junc_stats) with d6
df <- merge(d6, junc_stats, by = c("V", "J", "CDR3.nucleotide.sequence"))
df$Tda <- abs(as.integer(df$Vda)) + abs(as.integer(df$Nnt)) + abs(as.integer(df$Jda))
df$subrep <- substr(df$agg, 1, 5)
# For each combination of sample and sub-repertoire, 
# determine the number of TCR cDNA molecules or UMIs (count)
# in which the V-J junction was modified by 
# deletion or addition of a defined number of nucleotides. 
df1 <- df %>%
  group_by(Sample, subrep, Tda) %>%
  summarise(count = sum(Count)) %>%
  pivot_wider(id_cols = Tda, names_from = c("Sample", "subrep"),
              values_from = count, values_fill = 0 ) 

# Convert counts to per-sample frequency
for (i in 2:ncol(df1)) {
  df1[,i] <- df1[,i] / sum(df1[,i])}

# Arrange rows from lowest to highest number of nucleotides deleted or added
df1 <- arrange(df1, Tda)

# Pool values of total deletions or additions (Tda) > 15 together
df2 <- tail(df1, nrow(df1) - 16)
for (cn in 2:(ncol(df2))){
  df2[,cn] <- cumsum(df2[,cn])
}
df2z <- tail(df2,1)
df3 <- df1 %>%
  head(16) %>%
  bind_rows(df2z)
df3$Tda[nrow(df3)] <- ">15"
xo <- as.character(c(seq(0,15), ">15"))
df3$Tda <- factor(df3$Tda, levels = xo)

# Write a table for source data
vj_out <- df3
colnames(vj_out)[1] <- "Total deletions or additions"
write.table(vj_out, file = paste0(doa, toa,".source.txt"),
            sep = "\t", quote = F, row.names = F)

# Determine the cumulative frequency of UMIs 
# in which the V-J junction was modified by 
# deletion or addition of a defined number of nucleotides. 
for (cn in 2:(ncol(df3))){
  df3[,cn] <- cumsum(df3[,cn])
}

# Make a long dataframe for plotting
df4 <- df3 %>%
  pivot_longer(names_to = "sample_subrep", values_to = "freq", 
               2:(ncol(df3))) 
df4$subrep <- factor(substr(df4$sample_subrep,  nchar(df4$sample_subrep) - 4, nchar(df4$sample_subrep)), 
                    levels = COI)

# Summarise by sub-repertoire and determine the mean and range
df5 <- df4 %>% 
  group_by(subrep,  Tda) %>%
  summarise(
    mean_freq = mean(freq),
    max_freq = max(freq),
    min_freq = min(freq)
  )

df5$Tda <- factor(df5$Tda, levels = xo)
pal <- c("wtpri" = "royalblue1", "kopri" = "firebrick1")
# Plot it
pdf(file = paste0(doa, toa, ".pdf"), 
    height = 3.5, width = 9, useDingbats = F)
g <- ggplot(df5, aes(x=Tda, 
                     y=mean_freq, color = subrep, group = subrep)) + 
  geom_line () +
  theme_bw() +
  theme(
        axis.text.x =element_text(size = 8, angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20, hjust = .5, face = "italic"),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor.y  = element_blank(),
        legend.title = element_blank())+
  scale_x_discrete(labels=xo) +
  labs(y="Cumulative proportion of\nTCR cDNAs",
       x = "# nucleotides deleted or added at V-J junction") +
  scale_colour_manual(values = pal) +
  geom_errorbar(aes(ymin = min_freq, ymax = max_freq))
g
dev.off()
# The output pdf file above has labels. 
pdf(file = paste0(doa, toa, ".clean.pdf"), 
    height = 3.5, width = 7, useDingbats = F)
g <- ggplot(df5, aes(x=Tda, 
                     y=mean_freq, color = subrep, group = subrep)) + 
  geom_line () +
  theme_bw() +
  theme(axis.title = element_blank(),
    axis.text.x =element_text(size = 16, angle = 90, vjust = .5, hjust = 1),
    axis.text.y = element_blank(),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor.y  = element_blank(),
    legend.position = "none")+
  scale_x_discrete(labels=xo) +
  scale_colour_manual(values = pal) +
  geom_errorbar(aes(ymin = min_freq, ymax = max_freq))
g
dev.off()
# The output pdf file above has no y-axis labels or legend.
##################################################################################
# Get AUC for each combination of sample and sub-repertoire
df6 <- df3[, unique(df4$sample_subrep)]
# Get the sum of cumulative frequencies for each sample
df7 <- tibble(names(df6), colSums(df6))
colnames(df7) <- c("sample_subrep", "cumsum")

# To calculate AUC, divide the sum of cumulative frequencies 
# by the number of denominations
df7$auc <- df7$cumsum/nrow(df6)

# Prepare  for output
SROI <- c("Ctsl-dependent", "Ctsl-independent", "Newcomer")
df7$subrep <- factor(substr(df7$sample_subrep, nchar(df7$sample_subrep) - 4, nchar(df7$sample_subrep)), 
                     levels = COI)
df7$Subrepertoire <- SROI[match(df7$subrep, COI)]
df7$Sample <- substr(df7$sample_subrep, 1, str_locate(df7$sample_subrep, "_") -1)
df8 <- df7 %>%
select(Subrepertoire, Sample, auc) %>%
  rename(AreaUnderCurve = auc)
# Write a Source Data file showing the area under the curve results.
write.table(df8, file = paste0(doa, toa, "auc.stats.txt"), sep = "\t", quote = F, row.names = F)
