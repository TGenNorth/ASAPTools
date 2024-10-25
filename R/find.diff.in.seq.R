#' @title Find Differences between two sequence
#'
#' @description
#' Find differences between two sequences and return the differences separated by ",".
#' Note: This function needs to be adjusted to deal with IUPAC ambiguity
#' @param a A character reference string.
#' @param b A character string that is compared to the reference
#' @export find.diff.in.seq
#' @return A character string of differences within the sequence
#' @import tidyverse

find.diff.in.seq<- function(a, b){

  # a = "GTTGACCTACACAGGTGgCA"
  # b = "GTTGACCTACACAGGTGCCA"

  a <- toupper(a)
  b <- toupper(b)

  seq.a <- unlist(strsplit(a,split=""))
  seq.b <- unlist(strsplit(b,split=""))
  diff.d <- rbind(seq.a,seq.b)
  only.diff <- diff.d[,diff.d[1,]!=diff.d[2,]]
  pos <- which(diff.d[1,]!=diff.d[2,])

  Out <- as.data.frame(t(rbind(pos, as.data.frame(only.diff))))

  if(nrow(Out) == 0)
    {Out = "No Difference"} else {
    Out <- mutate(Out, SNP = paste(seq.a, `1`, seq.b, sep = "-"))
    Out <- paste(Out$SNP, collapse = ", ")
    }

  return(Out)
}
