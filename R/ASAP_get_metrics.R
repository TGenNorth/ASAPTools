#' Process ASAP XML List for General Metrics
#'
#' @param List of all XML objects to extract depth data from.
#' @return An object with general metrics by sample.

ASAP_get_metrics <- function(FinalList){

Out <- data.frame()
for (i in 1:length(FinalList)) {

  Temp <- cbind(Run = FinalList[[i]]$XML,
                Sequence_Name = FinalList[[i]]$Sequence_Name,
                Average_Depth = ifelse(is_empty(FinalList[[i]]$Average_Depth), 0, FinalList[[i]]$Average_Depth),
                Breadth = ifelse(is_empty(FinalList[[i]]$Breadth), 0, FinalList[[i]]$Breadth),
                Mapped = FinalList[[i]]$Mapped_Reads,
                Unmapped_Reads = FinalList[[i]]$Unmapped_Reads,
                Unassigned = FinalList[[i]]$Unassigned_reads)

  Out <- rbind(Out, Temp)

  #print(paste("Processing complete for:", i))
}

General_QC <- Out



General_QC$Breadth <- as.numeric(as.character(General_QC$Breadth))
General_QC$Average_Depth <- as.numeric(as.character(General_QC$Average_Depth))
General_QC$Mapped <- as.numeric(as.character(General_QC$Mapped))
General_QC$Unassigned <- as.numeric(as.character(General_QC$Unassigned))
General_QC$Unmapped_Reads <- as.numeric(as.character(General_QC$Unmapped_Reads))

return(General_QC)
}
