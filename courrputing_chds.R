'''


can you increase these values 
Social Economic  Events Admission Count
6     Low     High    High        No     8
14    Low     High Average       Yes     2
18    Low     High    High       Yes     3
22    Low     High     Low       Yes     3
'''

# Copy the original
chds_recorrupted <- chds_orig

# Define the combos to increase
boosts <- list(
  list(Social="Low", Economic="High", Events="High", Admission="No", n=8),
  list(Social="Low", Economic="High", Events="Average", Admission="Yes", n=2),
  list(Social="Low", Economic="High", Events="High", Admission="Yes", n=3),
  list(Social="Low", Economic="High", Events="Low", Admission="Yes", n=3)
)
for (b in boosts) {
  rows <- chds_orig[
    chds_orig$Social    == b$Social &
      chds_orig$Economic  == b$Economic &
      chds_orig$Events    == b$Events &
      chds_orig$Admission == b$Admission, ]
  
  # Duplicate 'n' times and add to dataset
  if (nrow(rows) > 0) {
    chds_recorrupted <- rbind(chds_recorrupted, rows[rep(1, b$n), ])
  }
}
