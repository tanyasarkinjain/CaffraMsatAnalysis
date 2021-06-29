require(gstudio)
require(ggplot2)
require(ggmap)
require(popgraph)
library(igraph)
require(maps)
require(raster)
require(fields)
library(dplyr)

#Reading Coordinate data
SADataCoords<-read.csv("/Users/tanyajain/Desktop/SouthernAfricaResearch/Wogan_etal_Cossypha_caffra_SamplingInfoTanya.csv",header=TRUE)
SADataCoords

#length of unique samples in Coords: 339
length(unique(SADataCoords$Sample))
#length of unique rows in Coords: 354
nrow(unique(SADataCoords))
SADataCoords[duplicated(SADataCoords$Sample, fromLast = TRUE) == TRUE,]
#deleting duplicates in Coords
SADataCoords <- SADataCoords[duplicated(SADataCoords$Sample, fromLast = TRUE) == FALSE,]
SADataCoords

#Sorts data by Taxon + Population to create the .gen file with more ease
SortedSADataCoords <- SADataCoords[
  with(SADataCoords, order(Taxon, Population)),
]
SortedSADataCoords

write.csv(SortedSADataCoords,"Tanya_Wogan_etal_Cossypha_caffra_SamplingInfoSorted.csv", row.names = FALSE)

#joining the data with the data in the igraph file
SADataMsat<-read_population(path="/Users/tanyajain/Desktop/SouthernAfricaResearch/cc2_Msats_GenePop.gen", type="genepop", locus.columns = 3:16, header=TRUE) 
SADataMsat
SADataMsat <- as.data.frame(SADataMsat)

#created 
SortedSADataMsat = inner_join(SortedSADataCoords, SADataMsat, by = c("Sample" = 'ID'), copy = FALSE, suffix = c(".x", ".y"),)
SortedSADataMsat
SortedSADataMsat <- SortedSADataMsat[, c(1, 40:53)]
SortedSADataMsat


SortedSADataMsat$Sample <- paste(SortedSADataMsat$Sample, ', ', sep=" ")

SortedSADataMsat[] <- lapply(SortedSADataMsat, function(x) sub("$^", "000:000", x))
SortedSADataMsat[]
SortedSADataMsat[] <- lapply(SortedSADataMsat, function(x) gsub(":", "", x))
SortedSADataMsat


#in genepop format
write.table(SortedSADataMsat, file = "SortedSADataMsat2.txt", sep = " ",
            row.names = F, col.names = T, quote = FALSE)


#Relay



