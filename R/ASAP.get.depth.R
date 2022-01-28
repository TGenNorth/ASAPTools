#' Process ASAP XML List for Depth Data
#'
#' @param List of all XML objects to extract depth data from.
#' @return An object with data from the ASAP xml.

ASAP.get.depth <- function(FinalList){

Depth_Out <- data.frame()

for (i in 1:length(FinalList)) {

  Run = FinalList[[i]]$XML

  Sequence_Name = FinalList[[i]]$Sequence_Name

  Average_Breadth <- FinalList[[i]]$Breadth

  Depth <- FinalList[[i]]$Depths

  Depth <- stringr::str_replace_all(Depth, "\"", "")

  Depth <- unlist(strsplit(Depth, ","))

  Depth <- as.numeric(Depth)

  Position <- seq(1:length(Depth))

  Temp <- as.data.frame(cbind(Run, Sequence_Name, Position, Depth, Average_Breadth))

  Depth_Out <- dplyr::bind_rows(Depth_Out, Temp)

  print(paste("Processing complete for:", Sequence_Name))
}

Depth_Out$Position <- as.numeric(as.character(Depth_Out$Position))
Depth_Out$Depth <- as.numeric(as.character(Depth_Out$Depth))

Depth_Out$Average_Breadth <- as.numeric(as.character(Depth_Out$Average_Breadth))

return(Depth_Out)
}
