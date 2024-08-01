#' Process ASAP DF for Proportions Data
#' @description
#' Function will use a user specified dataframe from read.ASAP function and parse the proportions data into a df that can be easily plotted.
#' @param read.ASAP.df imported df
#' @param num_cores number of cores to use for parallel processing
#' @return An object with data from the read.ASAP dataframe.
#' @export ASAP.get.proportions
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_replace
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel

ASAP.get.proportions <- function(read.ASAP.df, num_cores = 1) {
  library(dplyr)
  library(stringr)
  library(foreach)
  library(doParallel)

  # Register parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  proportions_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run <- read.ASAP.df$run[[i]]
    name <- read.ASAP.df$name[[i]]
    assay_name <- read.ASAP.df$assay_name[[i]]
    avg_breadth <- read.ASAP.df$breadth[[i]]
    proportions <- read.ASAP.df$proportions[[i]]

    # Remove any extra quotes and split by comma
    proportions <- unlist(strsplit(proportions, ","))
    proportions <- as.numeric(proportions)
    position <- seq(1, length(proportions))

    temp <- data.frame(run, name, assay_name, position, proportions, avg_breadth)

    print(paste("Processing complete for:", name))
    return(temp)
  }

  # Stop the parallel cluster
  stopCluster(cl)

  proportions_out$position <- as.numeric(as.character(proportions_out$position))
  proportions_out$proportions <- as.numeric(as.character(proportions_out$proportions))
  proportions_out$avg_breadth <- as.numeric(as.character(proportions_out$avg_breadth))

  return(proportions_out)
}
