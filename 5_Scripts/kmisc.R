# Install
remotes::install_github("stan-dev/cmdstanr")
options(repos = c(
  yihui = 'https://yihui.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'
))

install.packages('xfun')

# Load pakages and source file
library(xfun)      # for file_ext() function
library(ggplot2)   # for fggplot() function
library(scales)    # for scientific_format() function 
source("C:/Users/sima/Desktop/4b_MrBayes/2_mb_runs/plot_mrb.R")

# Run
plot_mrb()
