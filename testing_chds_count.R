# Frequency table of the full combination
combo_counts <- as.data.frame(table(chds_corrupted$Social,
                                    chds_corrupted$Economic,
                                    chds_corrupted$Events,
                                    chds_corrupted$Admission))

# Rename columns
names(combo_counts) <- c("Social", "Economic", "Events", "Admission", "Count")

# View only nonzero combinations
subset(combo_counts, Count > 0)

combo_counts <- as.data.frame(table(chds_orig$Social,
                                    chds_orig$Economic,
                                    chds_orig$Events,
                                    chds_orig$Admission))

# Rename columns
names(combo_counts) <- c("Social", "Economic", "Events", "Admission", "Count")

# View only nonzero combinations
subset(combo_counts, Count > 0)