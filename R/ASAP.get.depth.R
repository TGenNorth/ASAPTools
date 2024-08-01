#' Process ASAP DF for Depth Data
#' @description
#' Function will use a user specified dataframe from read.ASAP function and parse the depth data into a df that can be easily plotted.
#' @param read.ASAP.df imported df
#' @param num_cores number of cores to use for parallel processing
#' @return An object with data from the read.ASAP dataframe.
#' @export ASAP.get.depth
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_replace
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel

ASAP.get.depth <- function(read.ASAP.df, num_cores = 1) {
  library(dplyr)
  library(stringr)
  library(foreach)
  library(doParallel)

  # Register parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  depth_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run <- read.ASAP.df$run[[i]]
    name <- read.ASAP.df$name[[i]]
    assay_name <- read.ASAP.df$assay_name[[i]]
    avg_breadth <- read.ASAP.df$breadth[[i]]
    depth <- read.ASAP.df$depths[[i]]

    # Remove any extra quotes and split by comma
    depth <- unlist(strsplit(depth, ","))
    depth <- as.numeric(depth)
    position <- seq(1, length(depth))

    temp <- data.frame(run, name, assay_name, position, depth, avg_breadth)

    print(paste("Processing complete for:", name))
    return(temp)
  }

  # Stop the parallel cluster
  stopCluster(cl)

  depth_out$position <- as.numeric(as.character(depth_out$position))
  depth_out$depth <- as.numeric(as.character(depth_out$depth))
  depth_out$avg_breadth <- as.numeric(as.character(depth_out$avg_breadth))

  return(depth_out)
}
