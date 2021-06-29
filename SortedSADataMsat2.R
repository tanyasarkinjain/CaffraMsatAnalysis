require(gstudio)
require(ggplot2)
require(ggmap)
require(popgraph)
library(igraph)
require(maps)
require(raster)
require(fields)
require(comprehenr)
require(topoDistance)
require(raster)
require(rgdal)
library(raster)
library(sp)
library(rgdal)
require(gdistance)

#https://wiki.osgeo.org/wiki/South_African_Geodata
#https://www.jpl.nasa.gov/images/srtm-data-release-for-africa-colored-height

#------------------------------------------------------------

#expearimenting with topoDistance
#NEED: Digital elevation model, and a Raster layer for conductence?
wd <- ("/Users/tanyajain/Desktop/SouthernAfricaResearch/")
setwd(wd)
DEM <- raster(paste0(wd, "PIA04965.tif"))
DEM@crs
plot(DEM, col = terrain.colors(99))

PopulationCentroids2

xy <- matrix(ncol = 2, byrow = TRUE,c(-33.56606, 20.87105, -29.52967, 19.66842))
colnames(xy) <- c("latitude", "longitude")
xy

#tdist <- topoDist(DEM, xy, paths = TRUE)
#tdist

plot(DEM, main="r", xlab="Longitude (degrees)", ylab="Latitude (degrees)")

#---------------------------------------
#Reading MSAT data
SADataMsat2<-read_population(path="/Users/tanyajain/Desktop/SouthernAfricaResearch/SortedSADataMsat2removed.gen", type="genepop", locus.columns = 2:16, header=TRUE) 
SADataMsat2[3]
SADataMsat2
class(SADataMsat)
SADataMsat2
SADataMsat2$Population

SADataMsat2
SADataMsat211 <- SADataMsat2[SADataMsat2$Population == 'Pop-11', ]
SADataMsat211


#reading in CSV file 
SADataCoords2<-read.csv("/Users/tanyajain/Desktop/SouthernAfricaResearch/Tanya_Wogan_etal_Cossypha_caffra_SamplingInfoSortedremoved.csv",header=TRUE)
SADataCoords2


#Convert to MSAT data to multivariate: review: what exactly does this do?
SADataMsat2.mv <- to_mv(SADataMsat2, ploidy = 2, alleles=NA, drop.allele = FALSE, leave.as.na = FALSE)
SADataMsat2.mv
class(SADataMsat2.mv)



#Transform the multivariate data into population graph object

#removed populations with less than 3 individuals:Cossypha caffra iolaema-Tanzania_Iringa (1 sample), and Cossypha caffra iolaema	Mozambique_iolaema (2 individuals)
#removed: Cossypha dichroa	dichroa, Cossypha natalensis	natalensis
#Warning: Pop 11 and Pop 5 have fewer than 4 individuals (3). Choosing to procceed.

pops2 <- SADataMsat2$Population
pops2
SAGraphMsat2 <- popgraph(x = SADataMsat2.mv, groups = pops2)
SAGraphMsat2
V(SAGraphMsat2)

#Make the vertex names and group meaningful
#Cossypha natalensis population samples were taken both for malwai and South Africa: What should be thee region
length(unique(SADataCoords2$Population))
unique(SADataCoords2$Population)
V(SAGraphMsat2)$name <- c("SouthAfrica_EasternCape", "Malawi_Nyika", "DRC_Kivu", "Namibia_namaq", "SouthAfrica_namaq", "SouthAfrica_FreeState", "SouthAfrica_Gauteng","SouthAfrica_KwazuluNatal", "SouthAfrica_Limpopo", "SouthAfrica_Mpumalanga", "SouthAfrica_NorthernCape", "SouthAfrica_WesternCape", "Malawi_MtMulanje")

V(SAGraphMsat2)$group <- c('South Africa',           'Malawi',      'DRC',        'Namibia',       'South Africa',      'South Africa',         'South Africa',        'South Africa',            'South Africa',        'South Africa',           'South Africa',              'South Africa',           'Malawi')


#Checking vertex names, and edges 
V(SAGraphMsat2)$name
V(SAGraphMsat2)$group
E(SAGraphMsat2)

#DADataMsat Rows: 332, SADataCoords Rows: 332
nrow(SADataMsat2)
nrow(SADataCoords2)

#length of unique samples in Coords: 332
length(unique(SADataCoords2$Sample))
#length of unique rows in Coords: 332
nrow(unique(SADataCoords2))

#Selecting relevent columns from SADataCoords2 (not the loci data)
#Verified that coordinate data was correct
library(dplyr)
SADataCoordsR2 <- select(SADataCoords2, 'Sample', 'Latitude', 'Longitude', 'Genbank')
SADataCoordsR2 

#centroid data is WRONG
#Merge the Msat Data with the coordinate data
#verified
MsatCoordsDF2 <- merge(SADataMsat2, SADataCoordsR2, by.x = 'ID', by.y = 'Sample', sort = F)
MsatCoordsDF2

#find the centroids and pass on the data that has no latitude longitude data (No null latitude/longitude data)
PopulationCentroids2 <- aggregate(cbind(Latitude,Longitude) ~ Population, data = MsatCoordsDF2, FUN=mean, na.action=na.pass)
PopulationCentroids2

#set the population of the centroids to the main populations
unique(SADataCoords2$Population)

PopulationCentroids2$Population <- c('SouthAfrica_EasternCape', "Malawi_Nyika", "DRC_Kivu", "Namibia_namaq", "SouthAfrica_namaq", "SouthAfrica_FreeState", "SouthAfrica_Gauteng", "SouthAfrica_KwazuluNatal",
                                     "SouthAfrica_Limpopo", "SouthAfrica_Mpumalanga", "SouthAfrica_NorthernCape", "SouthAfrica_WesternCape" 
                                     ,"Malawi_MtMulanje" )
PopulationCentroids2

#Adding centroid spatial data based to graph (decorated based on Populations Centroids which coorespond to the stratum Population)
#lumps populations together but may be worth plotting individuals as well.


SAGraphMsatSP2 <- decorate_graph(SAGraphMsat2, PopulationCentroids2, stratum = "Population")
SAGraphMsatSP2

register_google(key = 'AIzaSyDQIIitnLwt_DfIJugW2ZrnEYYLIW6l2yQ')

#spatially plot the popgraph
OverallCentroid2 <- c(mean(V(SAGraphMsatSP2)$Longitude), mean(V(SAGraphMsatSP2)$Latitude))
OverallCentroid2
PopulationCentroids2
SAPopMap2 <-get_map(OverallCentroid2, maptype = "satellite", zoom = 4)
SAPopMap2<- ggmap(SAPopMap2) + geom_edgeset( aes(x=Longitude,y=Latitude), SAGraphMsatSP2, color="white") +  geom_nodeset( aes(x=Longitude, y=Latitude, color = group, size=3), SAGraphMsatSP2) 
SAPopMap2 + xlab("Latitude") + ylab("Longitude")

#plot the graph by itself without spatial data to understand cGD
plot(SAGraphMsat2, edge.color="black", vertex.label.color="darkred", vertex.color="#cccccc", vertex.label.dist=1)

#plotting using fruchterman.reingold layout
layout <- layout.fruchterman.reingold(SAGraphMsat2)
plot(SAGraphMsat2, layout=layout, edge.color="black", vertex.label.color="darkred", vertex.color="#cccccc", vertex.label.dist=4, size = .2)

#color by region 
#creating the graph from scratch: edges + nodes colored by region
p <- ggplot() + geom_edgeset( aes(x=Longitude,y=Latitude), SAGraphMsatSP2, color="darkgrey" )
p <- p + geom_nodeset( aes(x=Longitude, y=Latitude, color=group, size=size), SAGraphMsatSP2) 
p <- p + xlab("Longitude") + ylab("Latitude") 
p

#Matrix representation of the population graph object 
SAGraphMsatSP2
cGD2 <- to_matrix(SAGraphMsatSP2, mode="shortest path")
require(fields)
pDist2 <- rdist.earth( cbind(V(SAGraphMsatSP2)$Longitude, V(SAGraphMsatSP2)$Latitude ))
#creates a table of physical & cGD distance of the edges  between the vetricies 
dfSA2 <- data.frame( cGD=cGD2[upper.tri(cGD2)], Phys=pDist2[upper.tri(pDist2)])
#Despite a determined 8 connections, cGD and Phys contain 15 distances for all possible connections (including the ones deemed not present)
dfSA2$cGD
dfSA2$Phys



#Spearmans Coorelation test: Sees if there is a pattern of isolation by distance
'Spearmans Coorelation is a measure of the strength of a monotonic relationship between two variables. (rather than linear)
It is appropriate in this case, to examine, if, as, phsyical distance increase, 
cGD edge weight (interpopulation variance) increases. 
We expect, in an isolation by distance scenerio, monitically increasing graph.

Analysis:
The rho value of 0.5621088 indicates a "Strong" monitonically increasing
relationship between the phsyical and cGD distance.

The p-value, 1.492e-07 < 0.05, We reject the null hypothesis and conclude
that there is a statistically significant monotonic relationship between
the physical distance and cGD distance.

"A second issue potentially limiting the utility of pairwise genetic distance for landscape genetic analysis is
that the theoretical relationship between genetic distance (including linearized FST) and the spatial separation of populations is linear only when gene flow across
the landscape is a homogeneous function of spatial distance. This assumed homogeneity conflicts with a basic
tenet of landscape genetics that the movement of
migrants is influenced by ecological variables (e.g. suitable habitat, dispersal corridors, topography) whose
spatial distributions are decidedly heterogeneous." - Dyer et al 2010

Across the entire Population Graph,
cGD is estimated as the length of the shortest (geodesic)
path connecting pairs of populations.
'
cor.test( dfSA2$Phys, dfSA2$cGD, method="spearman")

#Graph of physical VS Conditional Genetic Distance using loess
qplot(Phys, cGD, geom="point", data=dfSA2) + stat_smooth(method="loess") + xlab("Physical Distance") + ylab("Conditional Genetic Distance")

#Graph of physical VS Conditional Genetic Distance us lm
qplot(Phys, cGD, geom="point", data=dfSA2) + stat_smooth(method = lm) + xlab("Physical Distance") + ylab("Conditional Genetic Distance")


#creating a 'loess' object to see which population connections fall above and below the loess fitted line
loessMod50 <- loess( cGD ~ Phys, data=dfSA2, span=0.50) # 50% smoothing span
dfSA2.edges <- as.data.frame(get.edgelist(SAGraphMsatSP2))
dfSA2.edges

#predicted cGD for each physical distance for each population
dfSA2.predicted <- dfSA2
dfSA2.predicted$PredictedcGD <- predict(loessMod50, dfSA2$Phys)
dfSA2.predicted

#https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-294X.2010.04748.x?casa_token=-lu4w1Mg6zYAAAAA:uI4JO_uJ5g1w6DnxpjPHlPjNJNdoDWIm7CfcBHJo7x6BcJVyd2Xb2IjPf6n7elVchfxNlBN1DhCYK2vQ
#Above provides an explanation of cGD

#Data Frame of edges where the actual cGD is less than predicted
dfSA2.lower <- dfSA2.predicted[dfSA2.predicted$cGD < dfSA2.predicted$PredictedcGD,]
dfSA2.lower$difference <- dfSA2.lower$PredictedcGD - dfSA2.lower$cGD
#attaches the edge connections to each row of the data frame according to the actual cGD2 (cGD2 contains a list of edge )
dfSA2.lower$edges <- to_vec(for(x in dfSA2.lower$cGD) paste((V(SAGraphMsatSP2)$name[which(cGD2 == x,arr.ind = TRUE)[1,][1]]), (V(SAGraphMsatSP2)$name[which(cGD2 == x, arr.ind = TRUE)[1,][2]]), sep = "_"))
dfSA2.lower

#Data Frame of Edges where the actual cGD is more than predicted
dfSA2.greater <- dfSA2.predicted[dfSA2.predicted$cGD > dfSA2.predicted$PredictedcGD,]
dfSA2.greater$difference <- dfSA2.greater$cGD - dfSA2.greater$PredictedcGD 
dfSA2.greater$edges <- to_vec(for(x in dfSA2.greater$cGD) paste((V(SAGraphMsatSP2)$name[which(cGD2 == x,arr.ind = TRUE)[1,][1]]), (V(SAGraphMsatSP2)$name[which(cGD2 == x, arr.ind = TRUE)[1,][2]]), sep = "_"))
dfSA2.greater

#Bootstrapping:
#ask about colinear variables

all_edge_vec <- to_vec(for(x in dfSA2$cGD) paste((V(SAGraphMsatSP2)$name[which(cGD2 == x,arr.ind = TRUE)[1,][1]]), (V(SAGraphMsatSP2)$name[which(cGD2 == x, arr.ind = TRUE)[1,][2]]), sep = "--"))
all_edge_vec


edge_list = list()

resamples <- lapply(1:100, function(i) MsatCoordsDF2[sample(nrow(MsatCoordsDF2), replace = TRUE), ])
resample1 <- resamples[[1]]
resample1_loci <- resample1[3:16]
class(resample1)
resample1.mv <- to_mv(resample1_loci)
resample1$Population

#original frequency dataframe
edge_freq <- data.frame(all_edge_vec, rep(0, 78))
edge_freq

resample1graph <- popgraph(resample1.mv, groups = resample1$Population)
V(resample1graph)$name <- c("SouthAfrica_EasternCape", "Malawi_Nyika", "DRC_Kivu", "Namibia_namaq", "SouthAfrica_namaq", "SouthAfrica_FreeState", "SouthAfrica_Gauteng","SouthAfrica_KwazuluNatal", "SouthAfrica_Limpopo", "SouthAfrica_Mpumalanga", "SouthAfrica_NorthernCape", "SouthAfrica_WesternCape", "Malawi_MtMulanje")
V(resample1graph)$group <- c('South Africa',           'Malawi',      'DRC',        'Namibia',       'South Africa',      'South Africa',         'South Africa',        'South Africa',            'South Africa',        'South Africa',           'South Africa',              'South Africa',           'Malawi')
V(resample1graph)

#keep track of the frequencys of the edges
#need a list of all possible edges


plot(resample1graph)

resample1graphSP <- decorate_graph(resample1graph, PopulationCentroids2, stratum = "Population")

#color by region 
#creating the graph from scratch: edges + nodes colored by region
p <- ggplot() + geom_edgeset( aes(x=Longitude,y=Latitude), resample1graphSP, color="darkgrey" )
p <- p + geom_nodeset( aes(x=Longitude, y=Latitude, color=group, size=size), resample1graphSP) 
p <- p + xlab("Longitude") + ylab("Latitude") 
p

#Get edge list and transform into -- format
e<-as_edgelist(resample1graph, names = TRUE)
e <- 
#to_vec(for(x in e) paste((V(SAGraphMsatSP2)$name[which(cGD2 == x,arr.ind = TRUE)[1,][1]]), (V(SAGraphMsatSP2)$name[which(cGD2 == x, arr.ind = TRUE)[1,][2]]), sep = "_"))
#update frequency dataframe with existant edges
#calculting edge 'stability': first need to obtain all possible edges
