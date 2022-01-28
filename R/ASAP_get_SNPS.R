#' Process ASAP XML List for SNPS
#'
#' @param List of all XML objects to extract SNP data from R parsed XML.
#' @return An object with data regarding SNPS from the ASAP xml.
#' @export ASAP_get_SNPS

ASAP_get_SNPS <- function(FinalList){

  SNPS_Out <- data.frame()

  for (i in 1:length(FinalList)) {

    Run = FinalList[[i]]$XML

    Sequence_Name = FinalList[[i]]$Sequence_Name

    Average_Breadth <- FinalList[[i]]$Breadth

    SNPS <- FinalList[[i]]$SNPS

    SNPS <- cbind(Run, Sequence_Name, Average_Breadth, SNPS)

    SNPS_Out <- rbind(SNPS_Out, SNPS)

    print(paste("Processing complete for:", Sequence_Name))
  }

  SNPS_Out$Position <- as.numeric(as.character(SNPS_Out$Position))
  SNPS_Out$Depth <- as.numeric(as.character(SNPS_Out$Depth))
  SNPS_Out$Average_Breadth <- as.numeric(as.character(SNPS_Out$Average_Breadth))
  SNPS_Out$Count <- as.numeric(as.character(SNPS_Out$Count))
  SNPS_Out$Percent <- as.numeric(as.character(SNPS_Out$Percent))

  return(SNPS_Out)
}
