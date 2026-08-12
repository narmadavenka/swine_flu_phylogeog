# Set up logging
sink("subsample_debug.log")
on.exit(sink())

trees_file <- "1990_USflu_30k.trees"
output_file <- "H3_1990-2025_thinned_1000_final.trees"

# Read file line by line to avoid loading everything into memory
cat("Reading trees file...\n")

# Extract posteriors and line numbers
posteriors <- c()
tree_lines <- c()
line_num <- 0

con <- file(trees_file, "r")
while (length(line <- readLines(con, n = 1)) > 0) {
  line_num <- line_num + 1
  
  if (grepl("^tree STATE", line)) {
    # Extract posterior probability (adjust regex to your format)
    posterior <- as.numeric(gsub(".*\\[&posterior=([0-9.]+)\\].*", "\\1", line))
    
    if (!is.na(posterior)) {
      posteriors <- c(posteriors, posterior)
      tree_lines <- c(tree_lines, line_num)
    }
  }
}
close(con)

cat("Total trees:", length(posteriors), "\n")

# Find top 10%
threshold <- quantile(posteriors, 0.90)
top_idx <- which(posteriors >= threshold)
cat("Top 10% threshold:", threshold, "\n")
cat("Trees in top 10%:", length(top_idx), "\n")

# Randomly sample 1000
set.seed(42)
final_idx <- sample(top_idx, size = min(1000, length(top_idx)), replace = FALSE)
final_lines <- tree_lines[final_idx]

cat("Sampled", length(final_lines), "trees\n")

# Write output
cat("Writing subsampled trees...\n")

# Copy header
con_in <- file(trees_file, "r")
con_out <- file(output_file, "w")

line_num <- 0
while (length(line <- readLines(con_in, n = 1)) > 0) {
  line_num <- line_num + 1
  
  # Copy header/metadata
  if (!grepl("^tree STATE", line)) {
    writeLines(line, con_out)
  } else if (line_num %in% final_lines) {
    # Write selected trees
    writeLines(line, con_out)
  }
}

close(con_in)
close(con_out)