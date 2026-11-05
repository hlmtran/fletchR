# Directory containing your R packages
path = .libPaths()[1]

# Find R package directories
descriptions = list.files(
  path,
  pattern = "^DESCRIPTION$",
  recursive = TRUE,
  full.names = TRUE
)

packages = dirname(descriptions)

# Get the latest modification time for each package
lastModified = sapply(packages, function(pkg) {
  files = list.files(pkg, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  files = files[!grepl("(^|/)\\.git(/|$)", files)]
  
  max(file.info(files)$mtime, na.rm = TRUE)
})

# Sort packages by most recent modification
result = data.frame(
  package = basename(packages),
  lastModified = as.POSIXct(lastModified, origin = "1970-01-01")
)

result = result[order(result$lastModified, decreasing = TRUE), ]

# Show packages modified in the last 7 days
subset(result, lastModified >= Sys.time() - 7 * 24 * 60 * 60)
