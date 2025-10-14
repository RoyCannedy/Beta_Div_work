## Frequency tables for all corrupted datasets (0–100% by 10%)

# Loop over saved corrupted files
for (p in seq(0, 100, by = 10)) {
  cat("\n==============================\n")
  cat(sprintf("Corruption Level: %d%%\n", p))
  cat("==============================\n")
  
  chds_corrupted <- read.csv(sprintf("CHDS_corrupted/CHDS_corrupted_%dperc.csv", p),
                             stringsAsFactors = TRUE)
  
  # Frequency table of the full combination
  combo_counts <- as.data.frame(table(chds_corrupted$Social,
                                      chds_corrupted$Economic,
                                      chds_corrupted$Events,
                                      chds_corrupted$Admission))
  
  # Rename columns
  names(combo_counts) <- c("Social", "Economic", "Events", "Admission", "Count")
  
  # View only nonzero combinations
  print(subset(combo_counts, Count > 0))
}


cat("\n==============================\n")
cat("Original (0% baseline)\n")
cat("==============================\n")


combo_counts <- as.data.frame(table(chds_orig$Social,
                                    chds_orig$Economic,
                                    chds_orig$Events,
                                    chds_orig$Admission))

# Rename columns
names(combo_counts) <- c("Social", "Economic", "Events", "Admission", "Count")

# View only nonzero combinations
subset(combo_counts, Count > 0)
