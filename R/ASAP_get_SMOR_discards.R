#' Process ASAP DF for SMOR discards
#' @description
#' Function will use a user specified dataframe from read.ASAP function and parse the SMOR discards data into a df that can be easily plotted.
#' @param read.ASAP imported df
#' @return An object with data from the read.ASAP dataframe.
#' @export ASAP_get_smor_discards
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_replace

ASAP_get_smor_discards <- function(read.ASAP.df){

  smor_discards_out <- data.frame()

  for (i in 1:nrow(read.ASAP.df)) {

    run = read.ASAP.df$run[[i]]

    name = read.ASAP.df$name[[i]]

    assay_name = read.ASAP.df$assay_name[[i]]

    avg_breadth <- read.ASAP.df$breadth[[i]]

    smor_discards <- read.ASAP.df$smor_discards[[i]]

    # discards <- stringr::str_replace_all(discards, "\"", "")

    smor_discards <- unlist(strsplit(smor_discards, ","))

    smor_discards <- as.numeric(smor_discards)

    position <- seq(1:length(smor_discards))

    temp <- as.data.frame(cbind(run, name, assay_name, position, smor_discards, avg_breadth))

    smor_discards_out <- dplyr::bind_rows(smor_discards_out, temp)

    print(paste("Processing complete for:", name))
  }

  smor_discards_out$position <- as.numeric(as.character(smor_discards_out$position))
  smor_discards_out$smor_discards <- as.numeric(as.character(smor_discards_out$smor_discards))

  smor_discards_out$avg_breadth <- as.numeric(as.character(smor_discards_out$avg_breadth))

  return(smor_discards_out)
}
