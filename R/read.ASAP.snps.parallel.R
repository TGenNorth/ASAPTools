#' @title Read ASAP SNP Results in Parallel
#'
#' @description
#' This function uses a user provided ASAP results XML path to parse sample and assay information from the ASAP output and display the SNPs.
#' @param XML A path with the ASAP XML file to parse.
#' @return An object with data from the ASAP xml.
#' @export read.ASAP.snps.parallel
#' @import xml2
#' @import tidyverse
#' @import foreach
#' @import doParallel
#' @import parallelly

read.ASAP.snps.parallel <- function(XML, cores){

  library(doParallel)
  library(foreach)

  cores <- cores
  cl <- makeCluster(cores) #not to overload your computer
  registerDoParallel(cl)

  Out <- data.frame()

  xml_data <- read_xml(XML)

  Run_Info = cbind(run = xml_attrs(xml_data))

  Iterations <- xml_length(xml_data)

  rm(xml_data)

  gc()

  #for each sample create a node list with specific child
  Out <- foreach(SAMPLE = 1:Iterations, .combine = rbind) %dopar% {

    library(xml2)
    library(tidyverse)

    xml_data <- read_xml(XML, options = "HUGE")

    #SAMPLE = 24
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

    Temp_Out <- data.frame()

    #For each assay within a sample
    for (ASSAY in 1:xml_length(Sample_Node)) {

      #ASSAY = 1
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

      #For each amplicon in assay:
      for(AMPLICON in 1:xml_length(Assay_Node)){

        #AMPLICON = 1
        #Get a single amplicon node
        Amplicon_Node <- xml_child(Assay_Node, AMPLICON)

        #xml_attrs(Amplicon_Node) #show the amplicon node attributes

        Breadth = xml_contents(xml_child(Amplicon_Node, "breadth"))

        Amplicon_Info <- cbind(
          #Add the amplicon number to the amp_info
          amplicon_number = AMPLICON
          #Parse out number of reads in amplicon
          #Amplicon_reads = xml_attrs(Amplicon_Node, "reads"),
          #Parse out the breadth of coverage
          #Breadth = as.character(xml_contents(xml_child(Amplicon_Node, "breadth"))),
          #Parse out the average breadth
          #Avg_Depth = as.character(xml_contents(xml_child(Amplicon_Node, "average_depth")))
        )

        for (SNP in 1:xml_length(Amplicon_Node)) {
          #SNP = 1
          #Get a single amplicon node
          SNP_Node <- xml_child(Amplicon_Node, SNP)

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
            Temp_Out <- rbind(Temp_Out, Temp)
          }
        }
      }
      gc()
      }
    Temp_Out
  }
  stopCluster(cl)
  return(Out)
}
