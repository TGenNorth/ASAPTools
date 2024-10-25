#' Process ASAP DF for SMOR discards in parallel
#' @description
#' Function will use a user specified dataframe from read.ASAP function and parse the SMOR discards data into a df that can be easily plotted.
#' @param read.ASAP.df imported df
#' @param num_cores number of cores to use for parallel processing
#' @return An object with data from the read.ASAP dataframe.
#' @export ASAP.get.smor.discards
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_replace
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel

ASAP.get.smor.discards <- function(read.ASAP.df, num_cores = 1) {
  library(dplyr)
  library(stringr)
  library(foreach)
  library(doParallel)

  # Register parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  smor_discards_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run <- read.ASAP.df$run[[i]]
    name <- read.ASAP.df$name[[i]]
    assay_name <- read.ASAP.df$assay_name[[i]]
    avg_breadth <- read.ASAP.df$breadth[[i]]
    smor_discards <- read.ASAP.df$smor_discards[[i]]

    # Remove any extra quotes and split by comma
    smor_discards <- unlist(strsplit(smor_discards, ","))
    smor_discards <- as.numeric(smor_discards)
    position <- seq(1, length(smor_discards))

    temp <- data.frame(run, name, assay_name, position, smor_discards, avg_breadth)

    print(paste("Processing complete for:", name))
    return(temp)
  }

  # Stop the parallel cluster
  stopCluster(cl)

  smor_discards_out$position <- as.numeric(as.character(smor_discards_out$position))
  smor_discards_out$smor_discards <- as.numeric(as.character(smor_discards_out$smor_discards))
  smor_discards_out$avg_breadth <- as.numeric(as.character(smor_discards_out$avg_breadth))

  return(smor_discards_out)
}
