#' Process ASAP DF for Proportions Data
#' @description
#' Function will use a user specified dataframe from read.ASAP function and parse the proportions data into a df that can be easily plotted.
#' @param read.ASAP imported df
#' @return An object with data from the read.ASAP dataframe.
#' @export ASAP_get_proportions
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_replace

ASAP_get_proportions <- function(read.ASAP.df){

  proportions_out <- data.frame()

  for (i in 1:nrow(read.ASAP.df)) {

    run = read.ASAP.df$run[[i]]

    name = read.ASAP.df$name[[i]]

    assay_name = read.ASAP.df$assay_name[[i]]

    avg_breadth <- read.ASAP.df$breadth[[i]]

    proportions <- read.ASAP.df$proportions[[i]]

    # proportions <- stringr::str_replace_all(proportions, "\"", "")

    proportions <- unlist(strsplit(proportions, ","))

    proportions <- as.numeric(proportions)

    position <- seq(1:length(proportions))

    temp <- as.data.frame(cbind(run, name, assay_name, position, proportions, avg_breadth))

    proportions_out <- dplyr::bind_rows(proportions_out, temp)

    print(paste("Processing complete for:", name))
  }

  proportions_out$position <- as.numeric(as.character(proportions_out$position))
  proportions_out$proportions <- as.numeric(as.character(proportions_out$proportions))

  proportions_out$avg_breadth <- as.numeric(as.character(proportions_out$avg_breadth))

  return(proportions_out)
}
