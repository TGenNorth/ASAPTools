#' @title Read ASAP Significance Messages
#' @description
#' This function uses a user provided ASAP results XML path to parse sample and assay information from the ASAP output and display the significance.
#' @param XML A path with the ASAP XML file to parse.
#' @return An object with data from the ASAP xml.
#' @export read.ASAP.significance
#' @import xml2
#' @import tidyverse

library(xml2)
library(tidyverse)

XML <- "/scratch/djasso-selles/NGARD/Su_Tuberculosis_project2021/17Aug2022_TannerASAP/17Aug2022_TannerASAP_analysis.xml"

read.ASAP.significance <- function(XML){

  Out <- data.frame()

  xml_data <- read_xml(XML)

  Run_Info = cbind(run = xml_attrs(xml_data))

  #for each sample create a node list with specific child
  for (SAMPLE in 7:8) { #1:xml_length(xml_data)

    SAMPLE <- 7

    #Subset xml_data
    Sample_Node <- xml_child(xml_data, SAMPLE)

    #xml_attrs(Sample_Node) Get the attributes in the node

    #Pull out the attributes of interest
    Sample_Info <- cbind(
      name = xml_attr(Sample_Node, "name")
      #Mapped_Reads = xml_attr(Sample_Node, "mapped_reads"),
      #unassigned_reads = xml_attr(Sample_Node, "unassigned_reads"),
      #unmapped_reads = xml_attr(Sample_Node, "unassigned_reads")
    )

    #For each assay within a sample
    for (ASSAY in 1:xml_length(Sample_Node)) {

      ASSAY = 4
      #subset Sample node to assay
      Assay_Node <- xml_child(Sample_Node, ASSAY)

      #xml_attrs(Assay_Node) list the attrs for the assay

      #Pull out the attributes of interest
      Assay_Info <- cbind(
        assay_function = xml_attr(Assay_Node, "function"),
        assay_gene = xml_attr(Assay_Node, "gene"),
        assay_name = xml_attr(Assay_Node, "name"),
        assay_type = xml_attr(Assay_Node, "type")
      )

      as.character(xml_child(Amplicon_Node, "significance"))

      xml_attrs(Assay_Node, attr)

      #For each amplicon in assay:
      for(AMPLICON in 1:xml_length(Assay_Node)){

        AMPLICON = 1
        #Get a single amplicon node
        Amplicon_Node <- xml_child(Assay_Node, AMPLICON)

        #xml_attrs(Amplicon_Node) #show the amplicon node attributes

        Breadth = xml_contents(xml_child(Amplicon_Node, "breadth"))

        Amplicon_Info <- cbind(
          #Add the amplicon number to the amp_info
          amplicon_number = AMPLICON,
          amplicon_significance = NA,
          amplicon_significance_flag = NA
        )

        for (SNP in 1:xml_length(Amplicon_Node)) {
          #SNP = 1
          #Get a single amplicon node
          SNP_Node <- xml_child(Amplicon_Node, SNP)

          if (xml_name(SNP_Node) == "significance"){

            Amplicon_Info[1,2] <- ifelse(length(xml_text(SNP_Node)) > 0, length(xml_text(SNP_Node)) > 0, NA)

            Amplicon_Info[1,3] <- ifelse(length(as.character(xml_attrs(SNP_Node))) > 0, as.character(xml_attrs(SNP_Node)), NA)
            length(xml_text(SNP_Node)) > 0

          }

          if (xml_name(SNP_Node) == "snp"){

            xml_attrs(SNP_Node)

            xml_children(SNP_Node)

            if(xml_length(SNP_Node) == 2){
              #Get and format the SNP_Distirbution Data
              SNP_Dist = as.character(xml_child(SNP_Node, 2))
              SNP_Dist = str_remove(SNP_Dist, "<base_distribution ")
              SNP_Dist = str_remove_all(SNP_Dist, "\"")
              SNP_Dist = str_remove_all(SNP_Dist, "/>")
            }else{
              SNP_Dist = "No reads matching SNP"
            }
            #Get and format SNP information
            SNP_call = as.character(xml_child(SNP_Node, 1))
            SNP_call = str_remove(SNP_call, "<snp_call count=")
            SNP_call = str_remove_all(SNP_call, "\"")
            SNP_call = str_remove_all(SNP_call, "percent=")
            SNP_call = str_replace(SNP_call, ">", " ")
            SNP_call = str_remove_all(SNP_call, "</snp_call>")
            SNP_call = as.list(str_split(SNP_call, " "))

            #Create object with values of interest
            SNP_Info <- cbind(
              location_depth = xml_attr(SNP_Node, "depth"),
              snp_name = xml_attr(SNP_Node, "name"),
              snp_position = xml_attr(SNP_Node, "position"),
              snp_reference = xml_attr(SNP_Node, "reference"),
              snp_depth = SNP_call[[1]][[1]],
              snp_proportion = SNP_call[[1]][[2]],
              snp_call = SNP_call[[1]][[3]],
              snp_distribution = SNP_Dist
            )

            Temp <- cbind(Run_Info, Sample_Info, Assay_Info, Amplicon_Info, SNP_Info)

            Out <- rbind(Out, Temp)
          }


        }



      }
    }
    print(paste("Processing complete for:", xml_attr(Sample_Node, "name")))
  }

  row.names(Out) <- 1:nrow(Out)

  return(Out)

}
