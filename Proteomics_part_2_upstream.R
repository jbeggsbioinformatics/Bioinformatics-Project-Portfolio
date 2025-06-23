#Identification of Phosphopeptides in SARS-CoV and SARS-CoV-2 Infected Cells (Workshop 2)

#1.Set working directory

setwd("~/Dropbox/Bioinformatics Msc/LIFE 754 - Proteomics/life_754_assignment_part1")

#2. Read comet output data file
df = read.csv("Animal_20200325_TM_HStmtpro_CoV12_ph2_fr3.txt", skip = 1, header = TRUE, sep = "\t")

#3. Order by e-value
df = df[order(df$e.value),]

#4. Create empty vectors to hold counts
fp_v = vector(mode = "integer", length = nrow(df))
tp_v = vector(mode = "integer", length = nrow(df))
fdr_v = vector(mode = "numeric", length = nrow(df))
isphospho_v = vector(mode = "logical", length = nrow(df))

#5. Generate new columns for FDR thresholds
fdr_01_hits = rep("", nrow(df))
fdr_05_hits = rep("", nrow(df))
fdr_1_hits = rep("", nrow(df))

fdr_01_phospho = rep("", nrow(df))
fdr_05_phospho = rep("", nrow(df))
fdr_1_phospho = rep("", nrow(df))

fp = 0
total = 0 #intialise counter for total peptides and false postive counts

#5. For loop iterates through the rows to calculate FP, TP, FDR, and phospho status
for (i in 1:nrow(df)) {
  prot_acc = df[i, "protein"]
  
  if (substr(prot_acc, 1, 5) == "DECOY") { #searching for decoy peptide sequences
    fp = fp + 1
  } else {
    total = total + 1
  }
  
  pep_mod = df[i, "modifications"]
  isphospho_v[i] = grepl("79.96", pep_mod, fixed = TRUE) #identifies phosphopeptides based on weight of phosphorylation modification (79.96)
  
  fp_v[i] = fp
  tp_v[i] = total - fp_v[i]
  
  #6. Calculation of false discovery rate (FDR) thresholds
  
  fdr_v[i] = fp_v[i] / (tp_v[i] + fp_v[i])
  
  #7. Quantification of peptide spectrum matches (PSMs) within each FDR thresholds using series of conditional statements
  
  if (fdr_v[i] <= 0.01) {
    fdr_01_hits[i] = "1"
    if (isphospho_v[i]) {
      fdr_01_phospho[i] = "1"
    }
  }
  if (fdr_v[i] <= 0.05) {
    fdr_05_hits[i] = "1"
    if (isphospho_v[i]) {
      fdr_05_phospho[i] = "1"
    }
  }
  if (fdr_v[i] <= 0.1) {
    fdr_1_hits[i] = "1"
    if (isphospho_v[i]) {
      fdr_1_phospho[i] = "1"
    }
  }
}

#8. Results added as new columns in the dataframe

df$count_fp = fp_v
df$count_tp = tp_v
df$is_phospho = isphospho_v
df$fdr = fdr_v

#9. PSMs and phospho PSM hits at each FDR threshold added as new columns

df$fdr_01_hits = fdr_01_hits
df$fdr_05_hits = fdr_05_hits
df$fdr_1_hits = fdr_1_hits

df$fdr_01_phospho = fdr_01_phospho
df$fdr_05_phospho = fdr_05_phospho
df$fdr_1_phospho = fdr_1_phospho

#10. Make the top row show counts for hits and phospho counts within the thresholds

df[1, "fdr_01_hits"] = sum(fdr_v <= 0.01)
df[1, "fdr_05_hits"] = sum(fdr_v <= 0.05)
df[1, "fdr_1_hits"] = sum(fdr_v <= 0.1)

df[1, "fdr_01_phospho"] = sum(fdr_v <= 0.01 & isphospho_v)
df[1, "fdr_05_phospho"] = sum(fdr_v <= 0.05 & isphospho_v)
df[1, "fdr_1_phospho"] = sum(fdr_v <= 0.1 & isphospho_v)

#Make the other rows (after first) empty strings

df[2:nrow(df), c("fdr_01_hits", "fdr_05_hits", "fdr_1_hits", "fdr_01_phospho", "fdr_05_phospho", "fdr_1_phospho")] <- ""

#11. The modified data frame is saved as a new CSV for downstream analysis

write.csv(df, "Phosphopeptide_complete.csv", row.names = FALSE)

# Load the ggplot2 package for generation of scatter plot
library(ggplot2)

#12. Creation of the scatter plot showing True Positive Count vs FDR

ggplot(df, aes(x = tp_v, y = fdr_v)) +
  geom_point(color = "blue") + 
  labs(title = "Scatter Plot of True Positive Count vs FDR",
       x = "True Positive Count",
       y = "FDR") +
  theme_minimal()

  
  