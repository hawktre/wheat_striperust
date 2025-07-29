## load up the packages we will need:  (uncomment as required)
user_lib <- "~/R/x86_64-pc-linux-gnu-library/4.4"
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(user_lib)

# Set CRAN repository
# This is necessary for the Rscript to run on the HPC cluster
options(repos = c(CRAN = "https://cloud.r-project.org")) 

# Install required packages if not already installed
pkgs <- c("data.table", "here", "dplyr", "purrr", "stringr", "tidyr")

install.packages(pkgs, dependencies = TRUE)