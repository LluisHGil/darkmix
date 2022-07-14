library(knitr)
opts_knit$set(base.dir = '_static') # Change the base dir where to save figures
knit("docs/source/darkmix_steps.Rmd")
