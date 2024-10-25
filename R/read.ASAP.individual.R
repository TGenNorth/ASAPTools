#' @title Read ASAP Results for an Individual xml file
#'
#' @description
#' This function uses a user provided ASAP results XML path to parse sample and assay information from the ASAP output for an individual sample. SNPs will not be included in this function and the function read.ASAP.snps should be used instead.
#' @param XML A path with the individual ASAP XML file to parse.
#' @return An object with data from the ASAP xml.
#' @export read.ASAP.individual
#' @import xml2
#' @import tidyverse
#'

read.ASAP.individual <- function(XML){

  Out <- data.frame()

  xml_data <- read_xml(XML, options = "HUGE")

  Run_Info = cbind(run = "Individual_XML_processing")

  #SAMPLE = 5
  #Subset xml_data
  Sample_Node <- xml_data

  #xml_attrs(Sample_Node) Get the attributes in the node

    #Pull out the attributes of interest
  Sample_Info <- cbind(
      name = xml_attr(Sample_Node, "name"),
      mapped_reads = xml_attr(Sample_Node, "mapped_reads"),
      unassigned_reads = xml_attr(Sample_Node, "unassigned_reads"),
      unmapped_reads = xml_attr(Sample_Node, "unmapped_reads")
    )

    #For each assay within a sample
  for (ASSAY in 1:xml_length(Sample_Node)) {

      #ASSAY = 1
      #subset Sample node to assay
      Assay_Node <- xml_child(xml_data)

      #xml_attrs(Assay_Node) list the attrs for the assay

      #Pull out the attributes of interest
      Assay_Info <- cbind(
        assay_function = xml_attr(Assay_Node, "function"),
        assay_gene = xml_attr(Assay_Node, "gene"),
        assay_name = xml_attr(Assay_Node, "name"),
        assay_type = xml_attr(Assay_Node, "type")
      )

      #For each amplicon in assay:
      for(AMPLICON in 1:xml_length(Assay_Node)){

        #AMPLICON = 1
        #Get a single amplicon node
        Amplicon_Node <- xml_child(Assay_Node, AMPLICON)

        xml_attrs(Amplicon_Node) #show the amplicon node attributes

        Breadth = xml_contents(xml_child(Amplicon_Node, "breadth"))

        Amplicon_Info <- cbind(
          #Add the amplicon number to the amp_info
          amplicon_number = AMPLICON,
          #Parse out number of reads in amplicon
          amplicon_reads = xml_attrs(Amplicon_Node, "reads")[[1]],
          #Parse out the Amplicon Variant info
          amplicon_variant = ifelse(length(xml_attrs(Amplicon_Node, "reads")) == 2, xml_attrs(Amplicon_Node, "reads")[[2]], "No variant"),
          #Parse out the breadth of coverage
          breadth = ifelse(length(as.character(xml_contents(xml_child(Amplicon_Node, "breadth")))) == 0, "No Breadth", as.character(xml_contents(xml_child(Amplicon_Node, "breadth")))),
          #Parse out the average breadth
          avg_depth = ifelse(length(as.character(xml_contents(xml_child(Amplicon_Node, "average_depth")))) == 0, "No Average Depth", as.character(xml_contents(xml_child(Amplicon_Node, "average_depth")))),
          #Parse the consensus sequence
          consensus_seq = ifelse(length(as.character(xml_contents(xml_child(Amplicon_Node, "consensus_sequence")))) == 0, "No Consensus Sequence", as.character(xml_contents(xml_child(Amplicon_Node, "consensus_sequence")))),
          #Parse out the depths
          depths = ifelse(length(as.character(xml_contents(xml_child(Amplicon_Node, "depths")))) == 0, "No Depth", as.character(xml_contents(xml_child(Amplicon_Node, "depths")))),
          #Parse out the proportions
          proportions = ifelse(length(as.character(xml_contents(xml_child(Amplicon_Node, "proportions")))) == 0, "No Proportions", as.character(xml_contents(xml_child(Amplicon_Node, "proportions")))),
          #Parse out the discarded reads
          smor_discards = ifelse(length(as.character(xml_contents(xml_child(Amplicon_Node, "discards")))) == 0, "No SMOR Analysis", as.character(xml_contents(xml_child(Amplicon_Node, "discards"))))
          )

        Temp <- cbind(Run_Info, Sample_Info, Assay_Info, Amplicon_Info)

        Out <- rbind(Out, Temp)

    }
    print(paste("Processing complete for:", xml_attr(Sample_Node, "name")))
  }

  row.names(Out) <- 1:nrow(Out)

  Out[c(3:5,11,13,14)] <- lapply(Out[c(3:5,11,13,14)],as.numeric)

  return(Out)

}


