# TCR diversity analyses on entire repertoires 
# of M2 CD4SP thymocytes 
# from WT mice (n = 3 samples, with each sample pooled from 2 mice) 
# or Ctsl-DeltaTEC mice (n = 3 samples, with each sample pooled from 2 or 3 mice)
# are outputs of this script.

# In RStudio, install and/or load "tidyverse", "iNEXT", and "scales" 
if (!require("tidyverse"))
  install.packages("tidyverse")
library(tidyverse)

if (!require("iNEXT"))
  install.packages("iNEXT")
library(iNEXT)

if (!require("scales"))
  install.packages("scales")
library(scales)
# Set working directory to the folder that contains data
data_folder <- "" # **Instruction: Make "data_folder" the path to the folder containing the file called "fig3.clonotype.table.txt"
data_folder <- "/Users/daleys6/Documents/TCR_seq/ctsl/250103data/" 
# The preceding line is the path to the folder on the author's computer
setwd(data_folder)
df <- read.delim(file = "ctsl.tcra.clonotype.table.txt", header = T)

# For record keeping, set doa (date of analysis) and toa (type of analysis).
doa <- "250106" 
toa <- ".ctsl.inext"
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
df1 <- filter(df, 
              Genotype %in% GOI,
              Cell_type == CTOI,
              Tissue == TOI,
              Count >= min_count) 

# Collapse to clonotype (unique combination of V, J and CDR3aa sequence).
df2 <- df1 %>%
  group_by(V, CDR3.amino.acid.sequence, J, cw, Genotype) %>%
  summarise(count = sum(Count)) %>%
  rename(Sample = cw)

# Set order of factor(Sample)
# By determining the Genotype of each Sample
df2 %>% group_by(Genotype, Sample) %>% summarise(clones = n())
# Samples 45-47 are WT and Samples 48-50 are CtslDeltaTEC
unique(df3$Sample)
SOI <- c("30903b-045", "30903b-046",  "30903b-047", 
         "30903b-048", "30903b-049", "30903b-050" )
df2$Sample <- factor(df2$Sample, levels = SOI)
############################################################################
# iNEXT loop
# Warning: This command takes hours to execute; it's best to run it overnight.
# Let Hill number (q) = 1 to get Shannon Diversity
for (gn in 1:length(SOI)){
  df3 <- df2 %>%
    filter(Sample == SOI[gn]) %>%
    group_by(V, CDR3.amino.acid.sequence, J) %>%
    summarise (abundance = sum(count))
  
  y <- iNEXT(df3$abundance, q=1, datatype = "abundance", knots = 1000, se = TRUE, nboot=5)
  
  write.table(y[1], file = paste("DataInfo", SOI[gn], min_count, "txt",sep = "."), row.names = F, sep = "\t", quote = F)
  write.table(y[2], file = paste("iNextEst", SOI[gn], min_count, "txt",sep = "."), row.names = F, sep = "\t", quote = F)
  write.table(y[3], file = paste("AsyEst", SOI[gn], min_count, "txt",sep = "."), row.names = F, sep = "\t", quote = F)
}

######################################################

######################################################
# make a metadatafile
# open Terminal, navigate to working directory (wd), type ls > fil.txt and press enter.
# This will write filenames in wd into a file called "fil.txt".
rf1 <- read.delim(file = "fil.txt", header =F)
rf2 <- rf1 %>% 
  separate(V1,
           into = c("table", "sample", "min_count", "suffix"),
           sep = "[.]",
           remove = F)
rf2$suffix <- NULL

# Attach genotype 
atg <- rf2 %>% 
  filter(!(V1 == "fil.txt")) %>%
  group_by(sample) %>%
  summarise(rows = n()) %>%
  select(sample)
atg$genotype <- c(rep(GOI[1],3), rep(GOI[2],3))

colnames(rf2)[1] <- "file_name"
rf3 <- rf2 %>% 
  merge(atg) %>%
  filter(!(file_name == "fil.txt")) %>%
  select(file_name, genotype, sample, table, min_count)
write.table(rf3, file = "metadiversity.txt", sep = "\t", row.names = F, quote = F)
#######################################################
# make plots
inventory <- read.delim("metadiversity.txt", sep = "\t", header = TRUE)
metadatafile <- filter(inventory, table == "iNextEst")
filenames <- as.character( metadatafile$file_name )
all.the.data <- lapply (filenames, read.delim, sep = "\t", header=TRUE)
allmystuff <- mapply(`[<-`, all.the.data, 'file_name', value = metadatafile$file_name, SIMPLIFY = FALSE)
tf <- do.call("rbind", allmystuff)
OOI <- 1 # Set Hill number
tf1 <- merge(tf, metadatafile, by = "file_name") %>%
               filter(iNextEst.size_based.Order.q == OOI)
tf1$genotype <- factor(tf1$genotype, levels = GOI)
tf1$sample <- factor(tf1$sample, levels = SOI)

int <- tf1 %>% filter(iNextEst.size_based.Method == "Rarefaction")
obs <- tf1 %>% filter(iNextEst.size_based.Method == "Observed")
ext <- tf1 %>% filter(iNextEst.size_based.Method == "Extrapolation")
shap <- c ( "WT" = 16, "CtslDeltaTEC" = 16)

# plot species accumulation curves
titleS <- c("0" = "Every clonotype", 
            "1" = "Shannon diversity",
            "2" = "Simpson diversity")
pal <- c("royalblue1", "firebrick1")

pdf(file = paste("sac",doa, toa, min_count, OOI, "pdf", sep = "."), width = 3.85, height = 2.5)
g <- ggplot(tf1, aes(x=iNextEst.size_based.m, 
                     y=iNextEst.size_based.qD, 
                     group=sample, 
                     colour=genotype,  shape = genotype, fill = genotype)) +
  theme_bw() + 
 # geom_ribbon(aes(ymin = iNextEst.size_based.qD.LCL, 
  #                  ymax = iNextEst.size_based.qD.UCL), alpha = .2) +
  geom_line(data = int, size = .2) +
  #geom_line(data = ext, size = .2, linetype="dotted") +
  geom_point(data = obs, size = 2, stroke = .5) + 
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = shap) +
  scale_fill_manual(values = pal) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = .5, vjust = 1, size = 16, face = "bold"),
        legend.title = element_blank(),
        panel.grid = element_blank()
        ) +
  labs(y= "Shannon diversity (x 10^3)", x = "TCRa reads (x 10^3)") + 
  scale_y_continuous(breaks = c(1000, 2000, 3000, 4000),
                     labels = c(1:4)) +
  scale_x_continuous(trans = log2_trans(),
                     limits = c(2^10,max(tf1$iNextEst.size_based.m)),
                     breaks = c(4000, 32000, 256000),
                     labels = c(4, 32, 256)
                     )
g
dev.off()
# This plot has labels on the axes and a legend.
###############################################################################

pdf(file = paste("sac",doa, toa, min_count, OOI, "clean.pdf", sep = "."), width = 3, height = 2)
g <- ggplot(tf1, aes(x=iNextEst.size_based.m, 
                     y=iNextEst.size_based.qD, 
                     group=sample, 
                     colour=genotype,  shape = genotype, fill = genotype)) +
  theme_bw() + 
  #geom_ribbon(aes(ymin = iNextEst.qD.LCL, ymax = iNextEst.qD.UCL), alpha = .2) +
  geom_line(data = int, size = .2) +
  #geom_line(data = ext, size = .2, linetype="dotted") +
  geom_point(data = obs, size = 2, stroke = .5) + 
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = shap) +
  scale_fill_manual(values = pal) +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5, vjust = 1, size = 16, face = "bold"),
        legend.position = 'none',
        panel.grid = element_blank()
  ) +
  scale_y_continuous(breaks = c(1000, 2000, 3000, 4000),
                     labels = c(1:4)) +
  scale_x_continuous(trans = log2_trans(),
                     limits = c(2^10,max(tf1$iNextEst.size_based.m)),
                     breaks = c(4000, 32000, 256000),
                     labels = c(4, 32, 256)
  )
g
dev.off()

# This plot lacks labels on the axes and a legend.
# To make the figure for the paper, the axes and legend were annotated manually using Adobe Illustrator.  
###############################################################################################
# Make a file containing source data.
tf_out <- tf1 %>% filter(iNextEst.size_based.Method == "Observed" | iNextEst.size_based.Method == "Rarefaction")
tf_out$sample <- factor(tf_out$sample, levels = SOI)
tf_out <- tf_out %>% arrange(sample)
tfo <- tf_out[,c("genotype", "sample",
                 "iNextEst.size_based.m",
                 "iNextEst.size_based.Method",
                 "iNextEst.size_based.qD",
                 "iNextEst.size_based.SC")]

colnames(tfo) <- c("Genotype", "Sample", "TCRa_reads", 
                   "Observation_type",  "Shannon_diversity",
                   "Coverage")
write.table(tfo, file = paste(doa, toa, "sac_source.txt", sep = "."),
            sep = "\t", row.names = F, quote = F)

#############################################################################
# plot relative diversity 
ca <- tf1 %>%
  filter(iNextEst.size_based.Method == "Extrapolation") %>%
  group_by(sample) %>%
  summarise(
    max_obs = max(iNextEst.size_based.SC)
  )  %>%
  select(max_obs) %>%
  max()

cbb <- tf1 %>%
  filter(iNextEst.size_based.Method == "Observed") %>%
  group_by(sample) %>%
  slice(which.max(iNextEst.size_based.SC)) 

cb <- min(cbb$iNextEst.size_based.SC)

#base <- min(ca,cb)
# As sampling completeness (coverage) approached 1 for each sample (all > 0.96), 
# the “Relative diversity” values are the ratios of the Shannon diversities 
# "observed" in the data.
base <- cb
tf1$delta.c = tf1$iNextEst.size_based.SC - base

tf2 <- tf1 %>%
  filter(iNextEst.size_based.Method == "Observed") 

tf2$d_se <- (tf2$iNextEst.size_based.qD.UCL - tf2$iNextEst.size_based.qD) / 1.96

info <- tf1 %>% 
  select(sample, genotype) %>% 
  distinct(sample, genotype)
library(reshape2)
tf3 <- tf2 %>%
  dcast(sample ~ iNextEst.size_based.Order.q, value.var = "iNextEst.size_based.qD") %>%
  merge(info, by = "sample")


colnames(tf3)[2] <- "diversity"

tf4 <- tf2 %>%
  dcast(sample ~ iNextEst.size_based.Order.q, value.var = "d_se")
colnames(tf4)[ncol(tf4)] <- "stan_err"  

tf5 <- merge(tf3, tf4,by = c('sample'))

denom_div = mean(tf5$diversity[tf5$genotype == GOI[2]])
denom_se = mean(tf5$stan_err[tf5$genotype == GOI[2]])

tf6 <- tf5 %>% 
  mutate(
    div_ratio = diversity / denom_div,
    term_a  = (stan_err)^2 / (diversity)^2,
    term_b = (denom_se)^2 / (denom_div)^2,
    fac_b = (term_a + term_b)^.5,
    ci = div_ratio * fac_b * 1.96,
    dr_lcl=div_ratio - ci, 
    dr_ucl=div_ratio + ci
  ) %>% merge(metadatafile,c("sample", "genotype"))

tf6$genotype <- factor(tf6$genotype, levels = GOI)
tf6$sample <- factor(tf6$sample, levels = SOI)
tf6$xval <- match(tf6$genotype, GOI)
st <- t.test(tf6$div_ratio[tf6$genotype == GOI[1]],
             tf6$div_ratio[tf6$genotype == GOI[2]])

#plot it
pdf(file =paste("reldiv", doa, toa, min_count, OOI, "pdf",sep = "."), width = 2, height = 1.65)
gr <- ggplot(tf6, aes(y=div_ratio, x=xval, colour= genotype, 
                      shape = genotype)) +
  theme_bw() + 
  #geom_errorbar(aes(ymax = dr_ucl, ymin = dr_lcl), size = .2) +
  geom_jitter(size=2, width = .05) + 
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = shap) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(y = "Relative diversity") +
  scale_y_continuous(limits = c(.9 * min(tf6$div_ratio),1.1 * max(tf6$div_ratio))) +
  scale_x_continuous(limits = c(.5,2.5),
                     breaks = seq(1,2),
                     labels = GOI) +
  annotate('text', label = paste ("p =", round(as.numeric(st[3]),4)),
           x = 1.5, y = 2.1, size = 2)
gr
dev.off()
# The above plot has labels.
pdf(file = paste("reldiv", doa, toa, min_count, OOI, "clean.pdf",sep = "."), width = 2, height = 1.65)
gr <- ggplot(tf6, aes(y=div_ratio, x=xval, colour= genotype, 
                      shape = genotype)) +
  theme_bw() + 
  #geom_errorbar(aes(ymax = dr_ucl, ymin = dr_lcl), size = .2) +
  geom_jitter(size=2, width = .05) + 
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = shap) +
  theme(panel.grid  = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") + 
  scale_y_continuous(limits = c(.9 * min(tf6$div_ratio),1.1 * max(tf6$div_ratio))) +
  scale_x_continuous(limits = c(.5,2.5),
                     breaks = seq(1,2),
                     labels = GOI) 
gr
dev.off()
# The above plot is label free.
# Write table for Source Data
sd <- tf6[,c("genotype",
             "sample", 
             "diversity", 
             "stan_err",  
             "div_ratio")]
colnames(sd) <- c("Genotype",
                  "Sample",
                  "Shannon diversity",
                  "Standard error",
                  "Relative diversity")
write.table(sd, file =  paste(doa, toa,"reldiv", min_count, OOI, "txt", sep = "."),row.names = F, quote = F, sep = "\t")
###############################################################################
###########################################################################
