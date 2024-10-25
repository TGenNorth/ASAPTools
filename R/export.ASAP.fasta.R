#' @title Export ASAP Consensus Fasta
#'
#' @description
#' This function exports the consensus sequence from read.ASAP dataframe.
#' @param read.ASAP.df from read.ASAP with added column "Sequence_Name"
#' @param savepath path to save the file to.
#' @return Statement of export completion
#' @export export.ASAP.fasta
#' @import tidyverse
#'

export.ASAP.fasta <- function(read.ASAP.df, savepath) {

  # Initialize an empty character vector `fastaLines` to store the FASTA formatted lines.
  fastaLines = c()

  # Loop through each row in the data frame `read.ASAP.df`
  for (rowNum in 1:nrow(read.ASAP.df)) {

    # Append the FASTA header line to `fastaLines`.
    # The ">" symbol is followed by the "Sequence_Name" from the current row.
    fastaLines = c(fastaLines, as.character(paste(">", read.ASAP.df[rowNum, "Sequence_Name"], sep = "")))

    # Append the sequence data line (consensus sequence) to `fastaLines`.
    fastaLines = c(fastaLines, as.character(read.ASAP.df[rowNum, "consensus_seq"]))
  }

  # Open a connection to the file specified by `savepath` for writing.
  fileConn <- file(savepath)

  # Write each element in `fastaLines` to the file as a separate line.
  writeLines(fastaLines, fileConn)

  # Close the file connection to save changes and free up resources.
  close(fileConn)

  cat(paste("Consensus sequences saved to:", savepath))
}
