#' @title Find primers in reference sequence
#'
#' @description
#' This function uses reference sequence and a list of primers to identify and search for the locations, returning a dataframe with the reference name, primer name, primer direction, primer sequence, searched direction, start location, end location, matching width, and returned sequence. Be aware that this function searches both the original and reverse complimented sequences.
#' @param fasta.path The path to the reference fasta file.
#' @param primer.names A list of primer names.
#' @param primer.direction A list of primer directions.
#' @param primer.list A list of primers to align to the reference file.
#' @param max.mismatch A number of mismatches to allow in search. Note: that the script records ambigous characters as mismatches, thus, ensure this number is greater than the number of ambigous characters.
#' @returns An dataframe with the reference name, primer name, primer direction, primer sequence, searched direction, start location, end location, matching width, and returned sequence.
#' @export find.primers
#' @import Biostrings
#' @import data.table
#' @import tidyverse
#' @import foreach

find.primers <- function(fasta.path, primer.names, primer.direction, primer.list, max.mismatch){

  # Read the FASTA file
  sequence_file <- Biostrings::readDNAStringSet(fasta.path)

  # Iterate over each sequence in the FASTA file
  bed_data <- foreach(i = 1:length(primer.list), .combine = rbind) %do% {

      PRIMER <- primer.names[[i]]
      PRIMER_DIRECTION <- primer.direction[[i]]
      PRIMER_SEQUENCE <- primer.list[[i]]

      Combine_Out <- data.frame()

      for (GENE in 1:length(sequence_file)) {

        sequences <- sequence_file[[GENE]]

        ASSAY <- names(sequence_file)[[GENE]]

        Temp <- matchPattern(pattern = PRIMER_SEQUENCE, subject = as.character(sequences),
                       max.mismatch=max.mismatch, with.indels=T, fixed=TRUE,
                       algorithm="auto")

        if(nrow(data.frame(Temp)) > 0) {

          Out <- cbind(assay = ASSAY,
                       primer = PRIMER,
                       direction = PRIMER_DIRECTION,
                       primer.sequence = PRIMER_SEQUENCE,
                       primer.search = PRIMER_SEQUENCE,
                       searched.direction = "Forward",
                       data.frame(Temp))

        }else{
            Out <- cbind(assay = ASSAY,
                         primer = PRIMER,
                         direction = PRIMER_DIRECTION,
                         primer.sequence = PRIMER_SEQUENCE,
                         primer.search = PRIMER_SEQUENCE,
                         searched.direction = "Forward",
                         start = NA,
                         end = NA,
                         width = NA,
                         seq = NA)
        }

        ## Do reverse compliment search
        Temp <- matchPattern(pattern = as.character(reverseComplement(DNAString(PRIMER_SEQUENCE))), subject = as.character(sequences),
                             max.mismatch=max.mismatch, with.indels=T, fixed=TRUE,
                             algorithm="auto")

        if(nrow(data.frame(Temp)) > 0) {

          Out <- rbind(Out, cbind(assay = ASSAY,
                       primer = PRIMER,
                       direction = PRIMER_DIRECTION,
                       primer.sequence = PRIMER_SEQUENCE,
                       primer.search = as.character(reverseComplement(DNAString(PRIMER_SEQUENCE))),
                       searched.direction = "Reverse",
                       data.frame(Temp)))

        }else{
          Out <- rbind(Out, cbind(assay = ASSAY,
                       primer = PRIMER,
                       direction = PRIMER_DIRECTION,
                       primer.sequence = PRIMER_SEQUENCE,
                       primer.search = as.character(reverseComplement(DNAString(PRIMER_SEQUENCE))),
                       searched.direction = "Reverse",
                       start = NA,
                       end = NA,
                       width = NA,
                       seq = NA))
        }
        Combine_Out <- rbind(Combine_Out, Out)
      }

      Combine_Out
  }

  bed_data$Mismatches <- NA

  bed_data[complete.cases(bed_data)]

  which(complete.cases(bed_data$seq))

  for (ROW in which(complete.cases(bed_data$seq))) {
    bed_data$Mismatches[ROW] <- find.diff.in.seq(bed_data$primer.search[ROW], bed_data$seq[ROW])
  }

  return(bed_data)
}

###################################
# Testing
###################################

# fasta.path = "/scratch/tporter/TB_20240726_AmpSeq_Bedfile/reference.fasta"
# Primers <- read_csv("/scratch/tporter/TB_20240726_AmpSeq_Bedfile/temp_primer_list_ed.csv")
# primer.names = Primers$Names
# primer.direction = Primers$Direction
# primer.list = Primers$Primer
# max.mismatch = 2


# fasta.path <- "/scratch/tporter/RSV_20230402_SequencingMethodCompare_SQL/KT992094.fasta"
#
# primer.names <- c("RSVA-F1",
#                   "RSVA-R1",
#                   "RSVA-F2",
#                   "RSVA-R2",
#                   "RSVA-F3",
#                   "RSVA-R3",
#                   "RSVA-F4",
#                   "RSVA-R4",
#                   "RSVA-F5",
#                   "RSVA-R5",
#                   "RSVA-F6",
#                   "RSVA-R6",
#                   "RSVA-F7",
#                   "RSVA-R7",
#                   "RSVA-F8",
#                   "RSVA-R8",
#                   "RSVA-F9",
#                   "RSVA-R9",
#                   "RSVA-F10",
#                   "RSVA-R10")
#
# primer.direction <- rep(c("F", "R"), 10)
#
# primer.list <- c("ACGSGAAAAAATGCGTACAAC",
#                  "GAAGATTGTGCTATACCAAAATGAACA",
#                  "ACAGGCATGACTCTCCTGAT",
#                  "TTGGGTGTGGATATTTGTTTCAC",
#                  "GGGCAAATATGGAAACATACGTG",
#                  "GTTTGCYGAGGCTATGAATATGAT",
#                  "ACCTGGGACACTCTCAATCA",
#                  "GACATGATAGAGTAACTTTGCTGTCT",
#                  "GAACAACAGACTACTAGAGATTACCAG",
#                  "AGGAGTTTGCTCATGGCAA",
#                  "GTCACGAAGGAATCCTTGCA",
#                  "CCCTCTACCTCTTTTATTATGTAGAACC",
#                  "AGCTTAGGCTTAAGATGYGGA",
#                  "TGAGTTTGACCTTCCATGAGT",
#                  "GGTGTACAATCTCTATTTTCCTGGT",
#                  "CGATTAATAGGGCTAGTATCAAAGTG",
#                  "GGGTTGGTTCATCTACACAAGAG",
#                  "CGCAAYAATAAATTCCCTGCTCC",
#                  "CGTCTACAATGATTAGAACCAATTAC",
#                  "ACGAGAAAAAAAGTGTCAAAAACTAA")
#
# max.mismatch <- 2
#
# find.primers.df <- find.primers(fasta.path, primer.names, primer.direction, primer.list, max.mismatch)
