#' @title Plot SARS Coverage Depth
#' @name ASAP_plot_sars_depth
#' @description
#' This function creates plots to display the coverage across the genome of select ASAP XMLs.
#' @param XML_List A list of XML path.
#' @param XML_Names A list of XML names.
#' @param Out_Dir A directory of where the output should be saved.
#' @return Objects should be saved to specified directory.
#' @export ASAP_plot_sars_depth
#' @import tidyr
#' @import tidyverse

# library(tidyr)
# library(tidyverse)

ASAP_plot_sars_depth <- function(XML_List, XML_Names, Out_Dir){

    XML <- ASAPTools::ASAP_read_sars_xml(XML_List, XML_Names)

    Metrics <- ASAPTools::ASAP_get_metrics(XML)

    Depth_Data <- ASAPTools::ASAP_get_depth(XML)

    Metrics <- Metrics %>%
        mutate(Coverage = ifelse(Breadth > 90, ">90%", ifelse(Breadth > 80, "80-90%", ifelse(Breadth > 70, "70-80%", "< 70%"))))

    Metrics$Coverage <- factor(Metrics$Coverage, levels = c(">90%", "80-90%", "70-80%", "< 70%"))

    Depth_Data <- left_join(Depth_Data, Metrics, by = c("Run", "Sequence_Name"))

    Object <- Depth_Data %>%
      group_by(Run, Coverage, Position) %>%
      summarise(Median = median(Depth),
                Min = min(Depth),
                Max = max(Depth)) %>%
      ggplot(aes(x = Position))+
        geom_ribbon(aes(ymin = Min, ymax = Max, fill = Run), alpha = 0.25)+
        geom_line(aes(y = Median, color = Run))+
        facet_wrap(~Coverage, ncol = 1)+
        theme_bw()+
        scale_y_continuous(trans = "log10")+
        theme(legend.position = "bottom")+
        xlab("Genome position") +
        ylab("Coverage")+
        ggtitle("Median coverage across genomes", subtitle = "Shaded region represents the min max across sample groups.")

    ggsave(paste(Out_Dir, "SARS_Breadth_Coverage_Plot", XML_Names, ".png", sep = ""), plot = Object, height =10, width = 15)

    print("Analysis complete")

    return(Object)
}

