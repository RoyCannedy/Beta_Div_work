
#Social Economic  Events Admission Count
#6     Low     High    High        No     8
#14    Low     High Average       Yes     2
#18    Low     High    High       Yes     3
#22    Low     High     Low       Yes     3


# Copy the original
chds_orig <- read.csv("Desktop/betadiv_roy//CHDS.latentexample1.csv", stringsAsFactors = TRUE)

# Define combos to inflate
boosts <- list(
  list(Social="Low", Economic="High", Events="High",    Admission="No"),
  list(Social="Low", Economic="High", Events="Average", Admission="Yes"),
  list(Social="Low", Economic="High", Events="High",    Admission="Yes"),
  list(Social="Low", Economic="High", Events="Low",     Admission="Yes")
)

# Directory to save all corrupted versions
dir.create("CHDS_corrupted", showWarnings = FALSE)

# Loop over corruption percentages: 0%, 10%, ..., 100%
for (p in seq(0, 1, by = 0.1)) {
  chds_corrupted <- chds_orig
  
  for (b in boosts) {
    rows <- chds_orig[
      chds_orig$Social    == b$Social &
        chds_orig$Economic  == b$Economic &
        chds_orig$Events    == b$Events &
        chds_orig$Admission == b$Admission, ]
    
    # Calculate how many extra rows to add based on percentage of original group size
    if (nrow(rows) > 0) {
      n_to_add <- ceiling(nrow(rows) * p)
      chds_corrupted <- rbind(chds_corrupted, rows[rep(1:nrow(rows), length.out = n_to_add), ])
    }
  }
  
  out_csv <- sprintf("CHDS_corrupted/CHDS_corrupted_%dperc.csv", round(p * 100))
  write.csv(chds_corrupted, out_csv, row.names = FALSE)
  
  cat(sprintf("Saved %s | corruption level: %.0f%% | total rows: %d\n",
              out_csv, p * 100, nrow(chds_corrupted)))
}
