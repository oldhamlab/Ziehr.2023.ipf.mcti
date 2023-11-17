# targets.R

# setup -------------------------------------------------------------------

#pipeline packages
suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
})

# target options
tar_option_set(
  packages = c(
    'tidyverse',
    'patchwork',
    'scales',
    'SummarizedExperiment'
  ),
  format = 'rds',
  controller = crew::crew_controller_local(workers = 4)
)

# source functions
tar_source()


# targets -----------------------------------------------------------------

list(
  NULL
)
