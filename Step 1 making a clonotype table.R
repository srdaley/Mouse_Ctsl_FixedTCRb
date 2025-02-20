# I thank the first author, Elisabetta Petrozziello, and the senior author, Ludger Klein,
# for inviting me to contribute to this exciting research project. 

# Readers should be able to reproduce the figures that  
# report TCR sequence analyses on T cell samples from "Fixed-beta" mice
# using the files in the repository at https://github.com/srdaley/CathepsinL_TCRseq. 

# Download all files in the repository at https://github.com/srdaley/CathepsinL_TCRseq
# Create a path to the folder that contains the downloaded files.
ctsl_folder <- "" # **Instruction: Make "ctsl_folder" the path to the folder that contains the downloaded files.
# ctsl_folder <- "/Volumes/One Touch/QUT/TCR_seq/ctsl/250103data/" 
# The preceding line is the path to the folder on the author's computer

# The text below shows the code used to collate all sample datafiles into one clonotype table.
# To reproduce the results shown in a figure panel(s), e.g. Figs 3c and 3d, 
# use the code in the corresponding file, e.g "Figs 3c and 3d.R".

# In RStudio, install and/or load "tidyverse".
if (!require("tidyverse"))
  install.packages("tidyverse")
library(tidyverse)

# Load the data into one dataframe 
setwd(ctsl_folder)
metadatafile <- read.delim(file = "ctsl.metadata.txt", header = T)
filenames <- as.character( metadatafile$filename )
all.the.data <- lapply (filenames, read.csv, sep = ";", header=TRUE)
allmystuff <- mapply(`[<-`, all.the.data, 'cw', value = metadatafile$cw, SIMPLIFY = FALSE)
d <- do.call("rbind", allmystuff) %>% merge (metadatafile, by = "cw") 

# Exclude rows with non-productive sequences
d2 <- d %>% filter( Productive == "Productive" )

# Exclude rows with a CDR3 amino acid that does not signify one of the 20 common amino acids
d2$B <- str_detect(d2$CDR3.amino.acid.sequence, "B")
d2$Jb <- str_detect(d2$CDR3.amino.acid.sequence, "J")
d2$O <- str_detect(d2$CDR3.amino.acid.sequence, "O")
d2$U <- str_detect(d2$CDR3.amino.acid.sequence, "U")
d2$X <- str_detect(d2$CDR3.amino.acid.sequence, "X")
d2$Z <- str_detect(d2$CDR3.amino.acid.sequence, "Z")
LOI <- c("B", "Jb", "O", "U", "X", "Z")
d2$exc_letters <- rowSums(d2[,LOI])
d3 <- d2 %>%
  filter(exc_letters == 0)
 
# The dataset now contains 1426994 rows.
# Each row is a unique "TCR clone".
# A TCR clone is defined as a unique combination of 
# Sample (called 'cw'), "V", "CDR3.nucleotide.sequence", "J" and "C". 
# View the column names of "d3"
colnames(d3)
# Retain only the needed columns
d4 <- d3[,c("cw", "Mouse_No", "Sex", "Genotype", "Tissue", "Cell_type", 
            "V", 
            "J", 
            "CDR3.nucleotide.sequence", 
            "CDR3.amino.acid.sequence",  
            "Count")]
# d4 is the clonotype table that forms the starting point of all TCRseq analyses in the paper.
write.table(d4, file = "ctsl.tcra.clonotype.table.txt", row.names = F, quote = F, sep = "\t")
# The line above should write a file identical to "ctsl.tcra.clonotype.table.txt" downloadable from the repository on github.   
