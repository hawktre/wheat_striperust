# Set custom user library path
user_lib <- file.path(Sys.getenv("HOME"), "R_libs/4.4")

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

# Use only the custom user library
.libPaths(user_lib)

# Set CRAN repository (important for non-interactive jobs)
options(repos = c(CRAN = "https://cloud.r-project.org")) 

# Install required packages if not already installed
pkgs <- c("data.table", "here", "dplyr", "purrr")

install.packages(pkgs, dependencies = TRUE, Ncpus = 4)