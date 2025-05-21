setwd("/home/ycen/proj/rpkg/scDesignPop")
library(scDesignPop)
# Run once
# Remove docs/ from gitignore to ensure it is checked into git.
usethis::use_pkgdown()
# Run everytime you want to update your site
pkgdown::build_site()

# TODO:Missing alt-text in README.md and scDesignPop.rmd
pkgdown::clean_cache()
