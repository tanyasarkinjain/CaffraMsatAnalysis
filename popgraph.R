#Notes: probably a good idea to join the spatial data via a dataframe at the begining if that can be preseverd
#for bootstrapping this would be immensly helpful

#loading all requirements

'TODO: read about graph theory....
TODO: find a way to measure variation in graphs
TODO: look by region
'

#https://reader.elsevier.com/reader/sd/pii/S1055790314001808?token=F74B2BC878A446D3DA1D8F966E0B8C44E275B829C837CE38DE4A9F2D864A8A93FFDE88F4E82D7FBE2DB8C64F0A52EECC&originRegion=us-east-1&originCreation=20210617062521


require(gstudio)
require(ggplot2)
require(ggmap)
require(popgraph)
library(igraph)
require(maps)
require(raster)
require(fields)

#Reading MSAT data
SADataMsat<-read_population(path="/Users/tanyajain/Desktop/SouthernAfricaResearch/cc2_Msats_GenePop.gen", type="genepop", locus.columns = 3:16, header=TRUE) 
SADataMsat[3]
class(SADataMsat)

#Reading Coordinate data
SADataCoords<-read.csv("/Users/tanyajain/Desktop/SouthernAfricaResearch/Wogan_etal_Cossypha_caffra_SamplingInfoTanya.csv",header=TRUE)
SADataCoords

#Convert to MSAT data to multivariate: review: what exactly does this do?
SADataMsat.mv <- to_mv(SADataMsat, ploidy = 2, alleles=NA, drop.allele = FALSE, leave.as.na = FALSE)
SADataMsat.mv
class(SADataMsat.mv)

#Transform the multivariate data into population graph object
#Warnings: Pop-3, Pop-5, Pop-7 have fewer than 4 individuals
pops <- SADataMsat$Population
pops
length(pops)
SAGraphMsat <- popgraph(x = SADataMsat.mv, groups = pops)
SAGraphMsat
V(SAGraphMsat)

#Make the vertex names and group meaningful
#Cossypha natalensis population samples were taken both for malwai and South Africa: What should be thee region
V(SAGraphMsat)$name <- c("Cossypha caffra caffra","Cossypha caffra iolaema","Cossypha caffra kivuensis","Cossypha caffra namaquensis", "Cossypha natalensis", "Cossypha dichroa")
V(SAGraphMsat)$group <- c("South Africa", "Malawi", "Democratic Republic Congo", "Namibia","Malawi/SA", "South Africa" )

#Checking vertex names, and edges 
V(SAGraphMsat)$name
V(SAGraphMsat)$group
E(SAGraphMsat)

#DADataMsat Rows: 339, SADataCoords Rows: 354
nrow(SADataMsat)
nrow(SADataCoords)

'repeated samples in SADataCoords: 
Cossypha caffra namaquensi
Index: 90:224
Difference between duplicates: values in Citaiton
Chose to delete first occurences 
'
#length of unique samples in Coords: 339
length(unique(SADataCoords$Sample))
#length of unique rows in Coords: 354
nrow(unique(SADataCoords))
SADataCoords[duplicated(SADataCoords$Sample, fromLast = TRUE) == TRUE,]
#deleting duplicates in Coords
SADataCoords <- SADataCoords[duplicated(SADataCoords$Sample, fromLast = TRUE) == FALSE,]
SADataCoords

#Selecting relevent columns from SADataCoords
library(dplyr)
SADataCoordsR <- select(SADataCoords, 'Sample', 'Latitude', 'Longitude', 'Genbank')

#joining spatial dataframe with MSAT dataframe by ID/Sample
MsatCoordsDF <- merge(SADataMsat, SADataCoordsR, by.x = 'ID', by.y = 'Sample', sort = F)

#finding the centroid for each population using MSATCoordsDF
PopulationCentroids <- aggregate(cbind(Latitude,Longitude) ~ Population, data = MsatCoordsDF, FUN=mean, na.action=na.pass)
PopulationCentroids$Population <- c("Cossypha caffra caffra","Cossypha caffra iolaema","Cossypha caffra kivuensis","Cossypha caffra namaquensis", "Cossypha natalensis", "Cossypha dichroa")

#Adding centroid spatial data based to graph
SAGraphMsatSP <- decorate_graph(SAGraphMsat, PopulationCentroids, stratum = "Population")
SAGraphMsatSP

#check that the spatial vectors are present
V(SAGraphMsatSP)$Latitude
V(SAGraphMsatSP)$Longitude

register_google(key = '##')

#spatially plot the popgraph
OverallCentroid <- c(mean(V(SAGraphMsatSP)$Longitude), mean(V(SAGraphMsatSP)$Latitude))
OverallCentroid
PopulationCentroids
SAPopMap <-get_map(OverallCentroid, maptype = "satellite", zoom = 4)
SAPopMap<- ggmap(SAPopMap) + geom_edgeset( aes(x=Longitude,y=Latitude), SAGraphMsatSP, color="white") +  geom_nodeset( aes(x=Longitude, y=Latitude, color = group, size=3), SAGraphMsatSP) 
SAPopMap + xlab("Latitude") + ylab("Longitude")


#Matrix representation of the population graph object 
SAGraphMsatSP
cGD <- to_matrix(SAGraphMsatSP, mode="shortest path")
require(fields)
pDist <- rdist.earth( cbind(V(SAGraphMsatSP)$Longitude, V(SAGraphMsatSP)$Latitude ) )
#creates a table of physical & cGD distance of the edges  between the vetricies 
dfSA <- data.frame( cGD=cGD[upper.tri(cGD)], Phys=pDist[upper.tri(pDist)])
#Despite a determined 8 connections, cGD and Phys contain 15 distances for all possible connections (including the ones deemed not present)
dfSA$cGD
dfSA$Phys



#Spearmans Coorelation test: Sees if there is a pattern of isolation by distance
'Spearmans Coorelation is a measure of the strength of a monotonic relationship between two variables. (rather than linear)
It is appropriate in this case, to examine, if, as, phsyical distance increase, 
cGD edge weight (interpopulation variance) increases. 
We expect, in an isolation by distance scenerio, monitically increasing graph.

Analysis:
The rho value of 0.7107143 indicates a "Strong" monitonically increasing
relationship between the phsyical and cGD distance.

The p-value, 0.004055 < 0.05, We reject the null hypothesis and conclude
that there is a statistically significant monotonic relationship between
the physical distance and cGD distance.
'
cor.test( dfSA$Phys, dfSA$cGD, method="spearman")

#Graph of physical VS Conditional Genetic Distance using loess
qplot(Phys, cGD, geom="point", data=dfSA) + stat_smooth(method="loess") + xlab("Physical Distance") + ylab("Conditional Genetic Distance")

#Graph of physical VS Conditional Genetic Distance us lm
qplot(Phys, cGD, geom="point", data=dfSA) + stat_smooth(method = lm) + xlab("Physical Distance") + ylab("Conditional Genetic Distance")

#creating a 'loess' object to see which population connections fall above and below the loess fitted line
loessMod50 <- loess( cGD ~ Phys, data=dfSA, span=0.50) # 50% smoothing span

dfSA.edges <- as.data.frame(get.edgelist(SAGraphMsatSP))
dfSA.edges


#predicted cGD for each physical distance for each population
dfSA.predicted <- dfSA
dfSA.predicted$PredictedcGD <- predict(loessMod50, dfSA$Phys)

#Data Frame of edges where the actual cGD is less than predicted
dfSA.lower <- dfSA.predicted[dfSA.predicted$cGD < dfSA.predicted$PredictedcGD,]
dfSA.lower$difference <- dfSA.lower$PredictedcGD - dfSA.lower$cGD
dfSA.lower$edges <- to_vec(for(x in dfSA.lower$cGD) paste((V(SAGraphMsatSP)$name[which(cGD == x,arr.ind = TRUE)[1,][1]]), (V(SAGraphMsatSP)$name[which(cGD == x, arr.ind = TRUE)[1,][2]]), sep = "_"))
dfSA.lower

#Data Frame of Edges where the actual cGD is more than predicted
dfSA.greater <- dfSA.predicted[dfSA.predicted$cGD > dfSA.predicted$PredictedcGD,]
dfSA.greater$difference <- dfSA.greater$cGD - dfSA.greater$PredictedcGD 
dfSA.greater$edges <- to_vec(for(x in dfSA.greater$cGD) paste((V(SAGraphMsatSP)$name[which(cGD == x,arr.ind = TRUE)[1,][1]]), (V(SAGraphMsatSP)$name[which(cGD == x, arr.ind = TRUE)[1,][2]]), sep = "_"))
dfSA.greater

'Stratified Bootstrapping Procedure:
1. Removing rows for populations that have less than 4 samples as genepop is not reccomended in these cases.


Why not resample the largest data set by drawing samples of the size of the smallest one and compare the chosen statistics (mean for instance) : from the smallest set we have an empirical value, from the largest set we have a resampled distribution and under the null hypothsis, the empirical value should not be too atypical as compared to the distribution

TODO: look into how many samples you need for bootstrapping and stratified bootstrapping
read paper referenced in paper


Seqboot from PHYLLP package
Created something called a \'bootstrap support'.
' Stability of edges among geographic groupswas assessed 
using a bootstrap approach with 200 bootstrap pseu-doreplicates, 
which were generated using seqboot from the PHYLIPpackage (Felsenstein, 1989) 
and analysed as the original data set.The proportion of replicates where a
certain edge is found consti-tutes its bootstrap support. Edges with bootstrap
support of 50%or more were considered stable (Escobar García et al., 2012)
'

MsatCoordsDF

#Merge loci and spatial dataframe with datframe with population frequency informaion 
#https://www.statisticshowto.com/unequal-sample-sizes/


PopulationFrequency <- table(MsatCoordsDF$Population)
PopulationFrequncyDF <- as.data.frame(PopulationFrequency)
MsatCoordsDFFreq <- merge(MsatCoordsDF, PopulationFrequncyDF, by.x = 'Population', by.y = 'Var1', all.x = TRUE)
#removing rows that coorespond to a population with less than 10 samples
MsatCoordsDF10 <- MsatCoordsDFFreq[MsatCoordsDFFreq$Freq >= 3,]
#Checking populations have been removed
pops2 = unique(MsatCoordsDF10$Population)
MsatCoordsDF10


#resampling 1000 times
MsatCoordsDF10

#resamples are correct
resamples <- lapply(1:10, function(i) MsatCoordsDF10[sample(nrow(MsatCoordsDF10), replace = TRUE), ])
column_class(resamples[[1]], class = "locus")
resamples[[1]]

resamples.mv <- lapply((1:10), function(i) to_mv(resamples[[i]][3:16], ploidy = 2, alleles=NA, drop.allele = FALSE, leave.as.na = FALSE))

resamples.graphs <- lapply((1:10), function(i) popgraph(x = resamples.mv[[i]], groups = resamples[[i]]$Population))
plot(resamples.graphs[[1]])
#look at 'edge support'

