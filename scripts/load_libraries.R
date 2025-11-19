## This script will load all necessary libraries for all 
## analyses in this code base

library(patchwork)
library(metR)
library(cmocean)
library(rnaturalearthdata)
library(rnaturalearth)
library(reticulate)
library(lmodel2)
library(pheatmap)
library(viridis)
library(dunn.test)
library(purrr)
library(RColorBrewer)
library(ggthemes)
library(broom)
library(tidyverse)
library(lubridate)
library(data.table)

## Setting up a plotting theme for transcriptome-related figures

###Theme code
theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = 'black'),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major.x  = element_blank(),
            panel.grid.major.y  = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_blank(),
            strip.text = element_text(face="bold")
    ))
  
}

