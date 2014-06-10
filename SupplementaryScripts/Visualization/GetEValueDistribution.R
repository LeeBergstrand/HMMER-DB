#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand 
# Descript: Creates plot showing E-value distribution of HMMHits from HMMER-DB
#             
# Requirements: 
#-----------------------------------------------------------------------------
#Imports:
library(ggplot2)
library(scales)
library(RSQLite)

#Functions:
reverselog_trans = function(base = exp(1)) {
  trans = function(x) -log(x, base)
  inv = function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

sqlite = dbDriver("SQLite")
HMMDB = dbConnect(sqlite, "/Users/lee/Data/SteriodHMMs/HMMDB.sqlite")

data = dbGetQuery(HMMDB, "SELECT
                            HMM_Hits.HMM_E_Value,
	                          HMM_Hits.HMM_Model,
                          	HMM_Data.HMM_Family
                          FROM
	                          HMM_Hits,
	                          HMM_Data
                          WHERE
	                          HMM_Hits.HMM_Model = HMM_Data.HMM_Model
                            AND HMM_Hits.HMM_Coverage > 0.50")

plotObj = ggplot(data, aes(x = HMM_E_Value, y = HMM_Model, colour = factor(HMM_Family)))
plotObj + geom_point(alpha = 1/10) + scale_x_continuous(trans = reverselog_trans(10)) + 
          facet_grid(HMM_Family ~ ., scales = "free") + theme_bw() + 
          theme(legend.position = "none", plot.background = element_rect(fill = "black"), 
                axis.text = element_text(colour = "white"),
                axis.text.y = element_text(angle = -45, hjust = 1),
                axis.text.x = element_text(angle = -90, hjust = 1),
                axis.ticks = element_line(colour = "white"))