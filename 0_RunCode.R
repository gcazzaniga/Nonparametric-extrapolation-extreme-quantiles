### Load packages ##############################################################

packages <- c("climod", "kdensity", "rstudioapi", "rtop", "reshape2")

for (package in packages) {
  if (package %in% rownames(installed.packages())) {
    do.call('library', list(package))
  } else {
    stop(
      sprintf(
        "Package %s is not installed.
                 See https://github.com/gcazzaniga/Nonparametric-extrapolation-extreme-quantiles.git for the list of the required packages",
        package
      )
    )
  }
}

### Set current working direction ##############################################

setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path))

### Variables ##################################################################

input <- read.table("input.txt", fill = TRUE, stringsAsFactors = FALSE)

distribution <- input[1, !is.na(input[1, ])]
colnames(distribution)[1] <- "name"
for (n_par in 2:length(distribution)) {
  colnames(distribution)[n_par] <- sprintf("par%d", n_par - 1)
}
n_sim_max <- as.numeric(input[2, !is.na(input[2, ])])
sample_len <- as.numeric(input[3, !is.na(input[3, ])])
return_period <- as.numeric(input[4, !is.na(input[4, ])])

### Run main code #############################################################

source("1_main.R")
