#' Process ASAP DF for Depth Data
#' @description
#' Function will use a user specified dataframe from read.ASAP function and parse the depth data into a df that can be easily plotted.
#' @param read.ASAP imported df
#' @return An object with data from the read.ASAP dataframe.
#' @export ASAP_get_depth
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_replace

ASAP_get_depth <- function(read.ASAP.df){

  depth_out <- data.frame()

  for (i in 1:nrow(read.ASAP.df)) {

    run = read.ASAP.df$run[[i]]

    name = read.ASAP.df$name[[i]]

    assay_name = read.ASAP.df$assay_name[[i]]

    avg_breadth <- read.ASAP.df$breadth[[i]]

    depth <- read.ASAP.df$depths[[i]]

    # Depth <- stringr::str_replace_all(Depth, "\"", "")

    depth <- unlist(strsplit(depth, ","))

    depth <- as.numeric(depth)

    position <- seq(1:length(depth))

    temp <- as.data.frame(cbind(run, name, assay_name, position, depth, avg_breadth))

    depth_out <- dplyr::bind_rows(depth_out, temp)

    print(paste("Processing complete for:", name))
  }

  depth_out$position <- as.numeric(as.character(depth_out$position))
  depth_out$depth <- as.numeric(as.character(depth_out$depth))

  depth_out$avg_breadth <- as.numeric(as.character(depth_out$avg_breadth))

  return(depth_out)
}
