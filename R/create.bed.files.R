#' @title Create an ASAP formatted bed file
#' @description
#' This function creates a bed file for ASAP from find.primers output.
#' @param find.primers.df A dataframe from the find.primers function.
#' @param save.path A path to automatically save the bed file to in a .bed format.
#' @returns A dataframe formatted as a bed file.
#' @export create.asap.bed.file
#' @import Biostrings
#' @import data.table
#' @import tidyverse
#' @import foreach

create.asap.bed.file <- function(find.primers.df, save.path = NA){

  # Use find.primers.df and filter to complete records (removes primers that were not found)
  Temp <- find.primers.df %>%
    na.omit()

  # Select columns needed for asap bed file input. https://stackoverflow.com/questions/17108191/how-to-export-proper-tsv
  Temp <- Temp %>%
    select(assay, primer, direction, start, end)

  # If save.path was specified save it to that directory.
  if(!is.na(save.path)){
    write.table(Temp, file= save.path, sep='\t')
  }

  # Return the asap bed file to R env.
  return(Temp)
}
