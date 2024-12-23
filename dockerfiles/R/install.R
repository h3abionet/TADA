## install script for R pkgs modified from https://github.com/davismcc/r-tidybioc-img/blob/master/install.R
install.packages("BiocManager")
## update installed packages
BiocManager::install()

pkgs <- c(
  "devtools",
  "RCurl",
  "tidyverse",
  "pander",
  "phangorn",
  "dplyr",
  # "dada2",
  "DECIPHER",
  "digest",
  "biomformat",
  "optparse",
  "yaml"
)

# developer (Github) URLs; dev branch only for now
dev <- c("benjjneb/dada2")

# check that desired packages are available
ap.db <- available.packages(contrib.url(BiocManager::repositories()))
ap <- rownames(ap.db)
pkgs_to_install <- pkgs[pkgs %in% ap]

# do not reinstall packages that are already installed in the image
ip.db <- installed.packages()
ip <- rownames(ip.db)
pkgs_to_install <- pkgs_to_install[!(pkgs_to_install %in% ip)]

BiocManager::install(pkgs_to_install)

if (length(dev) > 0) { lapply(dev, devtools::install_github, upgrade = "never") }

## just in case there were warnings, we want to see them
## without having to scroll up:
warnings()

if (!is.null(warnings()))
{
  w <- capture.output(warnings())
  if (length(grep("is not available|had non-zero exit status", w)))
    quit("no", 1L)
}
