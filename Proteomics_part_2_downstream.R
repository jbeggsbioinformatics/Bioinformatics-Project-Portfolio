##Identification of Phosphopeptides in SARS-CoV and SARS-CoV-2 Infected Cells (Workshop 2 - Downstream Analysis)

#1.Set working directory

setwd("~/Dropbox/Bioinformatics Msc/LIFE 754 - Proteomics/life_754_assignment_part1/Life_754_uploadfiles")
library(ggplot2)

#2. Create dataframe of upstream output file 

psm <- read.csv("phosphopeptide_complete.csv")

#3. PSMs filtered at 1%, 5%, 10% FDR for both phospho and non-phospho peptides

# FDR <= 0.01 (1%), FDR <= 0.05 (5%), FDR <= 0.10 (10%)

# 1% FDR
psm_fdr_01 <- psm[psm$fdr <= 0.01,]

# 5% FDR
psm_fdr_05 <- psm[psm$fdr <= 0.05,]

# 10% FDR
psm_fdr_10 <- psm[psm$fdr <= 0.10,]

#4. Filter for phosphopeptides (where is_phospho is TRUE)
psm_fdr_01_phospho <- psm_fdr_01[psm_fdr_01$is_phospho == TRUE,]
psm_fdr_05_phospho <- psm_fdr_05[psm_fdr_05$is_phospho == TRUE,]
psm_fdr_10_phospho <- psm_fdr_10[psm_fdr_10$is_phospho == TRUE,]

#5. Filter non-phosphopeptides (where is_phospho is FALSE)
psm_fdr_01_nonphospho <- psm_fdr_01[psm_fdr_01$is_phospho == FALSE,]
psm_fdr_05_nonphospho <- psm_fdr_05[psm_fdr_05$is_phospho == FALSE,]
psm_fdr_10_nonphospho <- psm_fdr_10[psm_fdr_10$is_phospho == FALSE,]

#  isphospho_v[i] = grepl("79.96"), meaning its phosphorylated

#6. Creation of summary table for phospho and non-phospho counts at each FDR threshold
summary_table <- data.frame(
  FDR_Threshold = c("1% FDR", "5% FDR", "10% FDR"),
  Phospho_Peptides = c(nrow(psm_fdr_01_phospho), nrow(psm_fdr_05_phospho), nrow(psm_fdr_10_phospho)),
  Non_Phospho_Peptides = c(nrow(psm_fdr_01_nonphospho), nrow(psm_fdr_05_nonphospho), nrow(psm_fdr_10_nonphospho))
)

#7. Table reordered for FDR_Threshold factor levels to reverse the order, increasing in stringency
summary_table$FDR_Threshold <- factor(summary_table$FDR_Threshold, levels = c("10% FDR", "5% FDR", "1% FDR"))

# Reshaping of the table for plotting efficiency

summary_table_long <- reshape(summary_table, 
                              varying = c("Phospho_Peptides", "Non_Phospho_Peptides"), 
                              v.names = "Count", 
                              timevar = "Peptide_Type", 
                              times = c("Phospho", "Non-Phospho"),
                              direction = "long")

#8. Bar Chart Visualisation of the data, generated using ggplot

ggplot(summary_table_long, aes(x = FDR_Threshold, y = Count, fill = Peptide_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Counts of PSMs at 1%, 5% and 10% False Discovery Rate (FDR)", 
       x = "FDR Threshold", 
       y = "Count of PSMs") +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.8), 
            vjust = -0.2, size = 4) # Add PSM count as text on top of bars



