# This code was used to identify "recurrent" TCRs,
# defined as being present in 3 or more samples,
# sorted from the thymus as "M2" CD4SP cells,
# and plot the frequency of each TCR clonotype in WT vs CtslDeltaTEC (KO) samples.

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
toa <- ".ctsl.subrep.xy"
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
# Make a list of TCR clonotypes detected exclusively in 3 or more samples.
# TCRs detected only in 3 Fixedβ Ctsl-delta-TEC samples will be red,
# TCRs detected only in 3 Fixedβ littermate control samples will be blue, and
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

# Filter TCR clonotypes detected in 3 or more samples
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

# Get columns necessary to reattach the count data in dc2
dc5 <- dc4 %>% select(V, CDR3.amino.acid.sequence, J, cat)

# dc5 contains the "recurrent" TCR clonotypes detected in 3 or more samples.
# Merge dc5 with dc2, which contains the count data.
d8 <- merge(dc5, dc2)

# FOR SCATTER PLOTS WE NEED:
# TCR frequency, defined as the frequency of UMIs encoding each TCR clonotype
# as a percentage of all UMIs encoding recurrent TCRs from samples with a given genotype. 

# Aggregate counts within genotype 
d8$wtrs <- rowSums(d8[,WTSOI])
d8$kors <- rowSums(d8[,KOSOI])

# Convert counts to frequency
d9 <- d8[,c(1:4,ncol(d8)-1,ncol(d8))]

for (cn in (ncol(d9)-1):ncol(d9)){
  d9[,cn] <- d9[,cn] / sum(d9[,cn])
}

# Add names of Subrepertoires
# Designate "Subrepertoires of Interest (SROI)"
SROI <- c("Ctsl-dependent", "Ctsl-independent", "Newcomer")
d9$Subrepertoire <- SROI[match(d9$cat, COI)]
sourcedata <- d9 %>% 
  select(Subrepertoire, V,  CDR3.amino.acid.sequence, J, wtrs, kors) %>%
  rename(TCRfreqWT = wtrs,
         TCRfreqCtslDeltaTEC = kors) %>%
  arrange(Subrepertoire, desc(TCRfreqWT), desc(TCRfreqCtslDeltaTEC))
write.table(sourcedata, file = paste0(doa, toa, ".source.txt"), sep = "\t", quote = F, row.names = F)

# Determine the number of TCR clonotypes in each colour-coded category ("cat")
ds1 <- d9 %>% 
  group_by(cat) %>% 
  summarise(clones = n()) 

# Prepare the "d9" dataframe for plotting  
d9$cat <- factor(d9$cat, levels = COI)

# Determine the minimum frequency observed
mf <- 1/max(colSums(d8[,c((ncol(d8)-1):ncol(d8))]))

# Give "private" TCRs a "frequency" 
# one order of magnitude lower than the minimum frequency observed.

for (rn in 1:nrow(d9)){
  if (d9$wtrs[rn] == 0){
    d9$wtrs[rn] <- (.1 * mf) 
  } else (d9$wtrs[rn] == d9$wtrs[rn])
}

for (rn in 1:nrow(d9)){
  if (d9$kors[rn] == 0){
    d9$kors[rn] <- (.1 * mf) 
  } else (d9$kors[rn] == d9$kors[rn])
}

# Log transform 
d9$wt_lt <- log10(d9$wtrs)
d9$ko_lt <- log10(d9$kors)

# For TCRs observed in one condition only, 
# vary the frequency in the "non-observed" condition.   

d9$shi <- runif(nrow(d9), -.25, .25)

for (rn in 1:nrow(d9)){
  if (d9$wt_lt[rn] == log10((.1 * mf))){
    d9$wt_lt[rn] <- d9$wt_lt[rn] + d9$shi[rn]
  } else (d9$wt_lt[rn] == d9$wt_lt[rn])
}


for (rn in 1:nrow(d9)){
  if (d9$ko_lt[rn] == log10((.1 * mf))){
    d9$ko_lt[rn] <- d9$ko_lt[rn] + d9$shi[rn]
  } else (d9$ko_lt[rn] == d9$ko_lt[rn])
}

pal <- c("royalblue1", "gray", "firebrick1" )

pdf(file = paste(doa, toa, min_count, "xy.scat.pdf", sep = "."),
    height = 4, width = 4, useDingbats = F)
g <- ggplot(d9, aes(x = ko_lt, y = wt_lt, color = cat)) +
  geom_point(shape = ".") +
  scale_x_continuous(breaks = c(log10(.1*mf),seq(-6,-1)),
                     labels = c("n.d.",1,10,100,1000,10000,100000),
                     limits = c(log10(.01*mf), .5 + max(d9$ko_lt))) +
  scale_y_continuous(breaks = c(log10(.1*mf),seq(-6,-1)),
                     labels = c("n.d.",1,10,100,1000,10000,100000),
                     limits = c(log10(.01*mf), .5 + max(d9$wt_lt))) +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = .2),
        plot.title = element_text(hjust = .5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Frequency in Ctsl^TEC (x10^-6)",
       y = "Frequency in WT (x10^-6)",
       title = paste(nrow(d9), "recurrent TCRs")) +
  geom_hline(yintercept = log10(mf/2), linetype = "dashed", linewidth =.1 ) +
  geom_vline(xintercept =  log10(mf/2), linetype = "dashed", linewidth =.1 ) +
  annotate("text", x = log10(.02*mf), y = .5 + max(d9$wt_lt), color = "royalblue1", label = paste(ds1$clones[ds1$cat == COI[1]])) +
  annotate("text", x = -2.5, y = .5 + max(d9$wt_lt), color = "black", label = paste(ds1$clones[ds1$cat == COI[2]], "'shared' TCRs")) +
  annotate("text", x = .2 + max(d9$ko_lt), y = log10(.05*mf), color = "firebrick1", label = paste(ds1$clones[ds1$cat == COI[3]])) 
g
dev.off()
# The graph produced above has labels

pdf(file = paste(doa, toa, min_count, "xy.scat.clean.pdf", sep = "."),
    height = 4, width = 4, useDingbats = F)
g <- ggplot(d9, aes(x = ko_lt, y = wt_lt, color = cat)) +
  geom_point(shape = ".") +
  scale_x_continuous(breaks = c(log10(.1*mf),seq(-6,-1)),
                     labels = c("n.d.",1,10,100,1000,10000,100000),
                     limits = c(log10(.01*mf), .5 + max(d9$ko_lt))) +
  scale_y_continuous(breaks = c(log10(.1*mf),seq(-6,-1)),
                     labels = c("n.d.",1,10,100,1000,10000,100000),
                     limits = c(log10(.01*mf), .5 + max(d9$wt_lt))) +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = .2),
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  geom_hline(yintercept = log10(mf/2), linetype = "dashed", size =.1 ) +
  geom_vline(xintercept =  log10(mf/2), linetype = "dashed", size =.1 ) +
  annotate("text", x = log10(.02*mf), y = .5 + max(d9$wt_lt), color = "royalblue1", label = paste(ds1$clones[ds1$cat == COI[1]])) +
  annotate("text", x = -2.5, y = .5 + max(d9$wt_lt), color = "black", label = paste(ds1$clones[ds1$cat == COI[2]], "'shared' TCRs")) +
  annotate("text", x = .2 + max(d9$ko_lt), y = log10(.05*mf), color = "firebrick1", label = paste(ds1$clones[ds1$cat == COI[3]])) 
g
dev.off()

# For subsequent analyses, write a "recurrent TCRs" clonotype table.
setwd(data_folder)
write.table(dc5, file = "ctsl.recurrent.TCRs.thresh3.txt", sep = "\t", row.names = F, quote = F )

# For subsequent Morisita-Horn Similarity Index calculations
# comparing with the "Ctsl-dependent", "shared", and "newcomer" TCR repertoires, 
# create new types of "Sample" called "wtpri", "shared", and "kopri"
dm <- d8 %>% 
  select(V, CDR3.amino.acid.sequence, J, cat, wtrs, kors) %>%
  mutate(count = wtrs + kors) %>%
  select(V, CDR3.amino.acid.sequence, J, cat, count) %>%
  rename(Sample = cat)
dm$Sample <- factor(dm$Sample, levels = COI)
dm <- arrange(dm, Sample)
setwd(data_folder)
write.table(dm, file = "Ctsl.recurrent.TCRs.3categories.thresh3.txt", sep = "\t", row.names = F, quote = F )
