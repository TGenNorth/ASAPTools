#' @title Import ASAP XML Results
#'
#' @description
#' This function uses a user provided ASAP results XML path and a user provided XML name to read in the XML into a matrix. Several other functions can then be used to parse the matrix.
#' @param XML_List A list of XML path.
#' @param XML_Names A list of XML names.
#' @return An object with data from the ASAP xml.
#' @export ASAP_read_sars_xml
#' @importFrom xml2 read_xml
#' @importFrom rvest html_nodes
#' @importFrom stringr str_split
#' @importFrom stringr str_replace
#' @importFrom tidyr separate

ASAP_read_sars_xml <- function(XML_List, XML_Names){

  Final_Out <- NULL

  for (n in 1:length(XML_List)) {

    XML <- XML_List[n]

    XML_Name <- XML_Names[n]

    Temp <- as.list(xml2::read_xml(XML, options = "HUGE"))

    TempNodes <- rvest::html_nodes(Temp, "sample")

    FinalList<-vector("list",length(TempNodes))

    for( File in 1:length(TempNodes)) {

      as.character(TempNodes[File])

      ### Pull Sample Data

      Temp_Info <- as.data.frame(stringr::str_split(as.character(TempNodes[File]), " "))

      Name <- gsub('^.{6}|.{1}$', '', Temp_Info[8,])

      Mapped_Reads <- gsub('^.{14}|.{1}$', '', Temp_Info[6,])

      Unassigned_Reads <- gsub('^.{18}|.{1}$', '', Temp_Info[10,])

      Unmapped_Reads <- gsub('^.{16}|.{3}$', '', Temp_Info[11,])

      Bam_File <- gsub('^.{10}|.{1}$', '', Temp_Info[2,])

      Average_Depth <- gsub('^.{15}|.{16}$', '', as.character(rvest::html_nodes(TempNodes[File], "average_depth")))

      Depths <- gsub('^.{8}|.{9}$', '', as.character(rvest::html_nodes(TempNodes[File], "depths")))

      Consensus_Sequence <- gsub('^.{20}|.{21}$', '', as.character(rvest::html_nodes(TempNodes[File], "consensus_sequence")))

      Breadth <- gsub('^.{9}|.{10}$', '', as.character(rvest::html_nodes(TempNodes[File], "breadth")))

      Gap_Filled_Consensus <- gsub('^.{30}|.{31}$', '',as.character(rvest::html_nodes(TempNodes[File], "gapfilled_consensus_sequence")))

      Proportions <- gsub('^.{13}|.{14}$', '',as.character(rvest::html_nodes(TempNodes[File], "proportions")))

      ###### Extracts SNPS into DF
      SNPS <- rvest::html_nodes(TempNodes[File], "snp")
      Out <- data.frame()
      for (i in 1:length(SNPS)) {

        SNP_Temp <- t(as.data.frame(stringr::str_split(as.character(SNPS)[i], pattern = c("<snp depth=|name=|position=|reference=|<snp_call count=|percent=|<base_distribution"))))[2:8]

        as.character(SNPS)[i]

        SNP_Temp <- as.data.frame(t(stringr::str_replace(SNP_Temp, "\\n|\"|\\\"", "")))
        names(SNP_Temp) <- c("Depth", "SNP_Name", "Position", "Reference", "Count", "Percent_Call", "Base_Distribution")
        SNP_Temp <- tidyr::separate(SNP_Temp, Percent_Call, into = c("Percent", "SNP"), sep = "\">")

        SNP_Temp

        SNP_Temp <- data.frame(lapply(SNP_Temp, function(x) {
          gsub("\"", "", x)
        }))

        SNP_Temp <- data.frame(lapply(SNP_Temp, function(x) {
          gsub(">\n", "", x)
        }))


        SNP_Temp <- data.frame(lapply(SNP_Temp, function(x) {
          gsub("</snp_call", "", x)
        }))

        SNP_Temp <- data.frame(lapply(SNP_Temp, function(x) {
          gsub("/</snp>", "", x)
        }))

        SNP_Temp <- data.frame(lapply(SNP_Temp, function(x) {
          trimws(x)
        }))

        names(SNP_Temp) <- c("Depth", "SNP_Name", "Position", "Reference", "Count", "Percent", "Call", "Base_Distribution")

        SNP_Temp

        Out <- rbind(Out, SNP_Temp)
      }

      Out

      Out_List <- list(XML = XML_Name, Sequence_Name = Name, Bam_File = Bam_File, Mapped_Reads = Mapped_Reads, Unassigned_reads = Unassigned_Reads, Unmapped_Reads = Unmapped_Reads, Breadth=Breadth, Average_Depth = Average_Depth, SNPS = Out, Depths = Depths, Proportions=Proportions, Consensus_Sequence = Consensus_Sequence, Gap_Filled_Consensus = Gap_Filled_Consensus)

      FinalList[[File]] <- Out_List
    }

    Final_Out <- c(Final_Out, FinalList)

  }
  return(Final_Out)
}



