#installing package from github so updates are included. Installing from CRAN does not work.
'install.packages("githubinstall")
install.packages("ggplot2"")
library(githubinstall)

install.packages( c("RgoogleMaps",
                    "geosphere",
                    "proto",
                    "sampling",
                    "seqinr",
                    "spacetime",
                    "spdep"), 
                  dependencies=TRUE )

library(devtools)
githubinstall("dyerlab/popgraph")
githubinstall("dyerlab/gstudio")'


#loading package
require(gstudio)
#
#info about package
help(package = "gstudio")
vignette("gstudio")
#system("open http://dyerlab.github.io/gstudio/")

#Genetic Data: creating a locus object
x <- locus(c(1, 2))
class(x)
x

#different types of loci
loc0 <- locus()
loc1 <- locus(1)

#can put in actual allele nitrogenous base
loc2 <- locus(c("C", "A"), phased = TRUE)
loc3 <- locus(c("C", "A"))
loc4 <- locus(12, type = "zyme")
loc5 <- locus(1:4)

#polyploid: more than 2 sets of chromosomes so more alleles at a specific loci
loc6 <- locus("A:C:C:G", type = "separated")
loc7 <- locus(1:1)

#c combines in R
loci <- c(loc0, loc1, loc2, loc3, loc4, loc5, loc6, loc7)
loci

#creating a data-frame of loci
df <- data.frame(ID = 0:7, Loci = loci)
df

#can do normal processing: this checks if any of the loci are 'empty'
is.na(df$Loci)

#checks if the loci are heterozygous
is_heterozygote(df$Loci)

#importing data
'notes about how data should be organized:
1. First line has info on data that is ignored by read_population()
2: names of loci are listed next
3: single line with word \'pop\' represents a new population for the popgraph.
4: identification number, comma, genotypes: per line for an indiviidual in a pop: loci are represented with 3 digits: 003005 (3:5) for diploid organism (i assumer the same order in which the loci names are listed?)
- ex) 1234, 003005'

#reading in different types of input
#[1] "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/gstudio/extdata"
#above path indicates where example files are placed on my computer

#alleles for a single locus (column)
#default, homozygous: single digit, 1.Heterozygous: seperated by colon, 1:2
file <- system.file("extdata", "data_2_column.csv", package = "gstudio")
data <- read_population(file, type = "column", locus.columns = 4:7)
data

#gametic phase: a particular association of alleles at different loci on the same chromosome. Gametic phase is influenced by genetic linkage.[2] If: arranged alphanumerically, linkage may not be clear: a:A b:B in reality may be A:a b:B where A and b are linked.
#captures linkage data
#exactly the same as reading in the data as above,however, there is the additional paramater 'phased.' Maintains order
file <- system.file("extdata", "data_2_column.csv", package = "gstudio")
data <- read_population(file, type = "column", locus.columns = 4:7, phased = TRUE)
data

#ALFP-like data: Amplified fragment length polymorphism: 0/1 indicates the presence/absense of a particular band after electrophoretic separation  after PCR of fragments.
file <- system.file("extdata", "data_aflp.csv", package = "gstudio")
data <- read_population(file, type = "aflp", locus.columns = c(4, 5))
data

#SNP Minor Allele Data: minor allele refers to the frequency in which the 2nd most common allele occurs in a population. B is the minor one
#alleles are 0 1 or 2 indicating number of minor alleles 
file <- system.file("extdata", "data_snp.csv", package = "gstudio")
data <- read_population(file, type = "snp", locus.columns = 4:7)
data

#zyme-like data: allozyme: allelic variants of enzymes encoded by structural genes, determined via electrophoresis. Loci can be determined
file <- system.file("extdata", "data_zymelike.csv", package = "gstudio")
data <- read_population(file, type = "zyme", locus.columns = 4:7)
data

#pre-seperated data for higher ploidy?
file <- system.file("extdata", "data_separated.csv", package = "gstudio")
data <- read_population(file, type = "separated", locus.columns = c(4, 5))
data


#Saving Data
save(df, file = "MyData.rda")
load("MyData.rda")
df
ls()

#saving as text
write_population(df, file = "~/Desktop/MyData.csv")

#https://rdrr.io/github/dyerlab/gstudio/src/R/write_population.R
#write_population(df, file = "~/Desktop/MyData.txt", mode = "genepop",  stratum="locus")
#write_population(df, file = "~/Desktop/MyData.str", mode = "structure",  stratum="locus")

#example data from sonaron beetle
#LTRS, WNT,EN, EF refer to locuses, and the values the specific alleles along with the frequency of those alleles
data(arapat)
summary(arapat)

#convience functions

#identifies locus column names
column_class(arapat, class = "locus")
#identifies locus column index
column_class(arapat, class = "locus", mode = "index")
#partitioning: returns data frame partioned by stratum
names(arapat)
clades <- partition(arapat, stratum = "Species")
names(clades)

#different number of individuals but same # of columns
lapply(clades, dim)

#example of using partitioning and lapply

#how to get frequency table to contain 0: modify either table or unlist
pops <- partition(arapat, stratum = "Population")

#____$____ gives the dimensions of the frequencies of the 2nd input categories for every 1st inputs
counts <- lapply(pops, function(x) {
  return(table(x$Species))
})

m <- matrix(unlist(counts), ncol = 3, byrow = TRUE)
rownames(m) <- names(pops)[1:17]
colnames(m) <- levels(arapat$Species)
m[1:10, ]

#plotting populations: Did not work
require(ggplot2)
plot(arapat)

#a different mechanism for plotting
#loading in the package
require(ggmap)
#chooses all rows, but specific columns 2,3,6 and 5, assigns to data
data <- unique(arapat[, c(2, 3, 6, 5)])
#gets the longitude and latitude columns and gets the mean for the columns (1 for row 2 for columns) 
#(delinated by column 2 titled cluster)
centroid <- apply(data[, 3:4], 2, mean)
centroid


#gstudio population_map did not work as intended: used ggmap directly to generate map in gstudio documentation
#https://www.nceas.ucsb.edu/sites/default/files/2020-04/ggmapCheatsheet.pdf

maptype = c("roadmap", "terrain", "satellite", "hybrid")
register_google(key = '##')

myLocation = c(lon = -111.39669, lat = 25.65789)
typeof(myLocation)

#https://journal.r-project.org/archive/2013-1/kahle-wickham.pdf
my_map <- get_map(location = myLocation, source = 'google', maptype = 'roadmap', crop = FALSE, zoom = 7) 
ggmap(my_map)+  geom_point(aes(x = Longitude, y = Latitude, color = Cluster), data = data, size = 3)


#Allele Frequencies
#the function frequencies constructs a table automatically

#single locus (gives the allele frequencies at that locus)
freqs <- frequencies(arapat$EF)
class(freqs)
freqs

#multilocus
freqs.loci <- frequencies(arapat)
freqs.loci[1:10, ]

#substrata and Allele Frequencies
#seperateas further by the stratum provided
freqs.strata <- frequencies(arapat, stratum = "Cluster")
freqs.strata[1:10, ]

#plotting allele frequencies
#plots the frequency for the allele specified
#arapat$MP20 : gives you a vector of alleles that coorespond to the locus MP20 
plot(arapat$MP20)
plot(arapat$MP20, mode = 'pie')

#color dictated by cluster. (aes = aesthetic mappings)
#ggplot intitializes a ggplot object
ggplot() + geom_locus(aes(x = MP20, fill = Cluster), data = arapat)

#freqs.strata was the  allele frequrncy by cluster and loci. Here, finds the frequency of a locus where locus must belong to MP20 and AML
#compares those two loci
f <- freqs.strata[freqs.strata$Locus %in% c("MP20", "AML"), ]
f
#the summary
'length: shows the number of individuals/rows
Frequency: Breaks down the distribution of allele frequencies'
summary(f)

ggplot(f) + geom_frequencies(f) + facet_grid(Stratum ~ .) + theme(legend.position = "none")

#frequency gradients: useful for observing how allele frequencies change in relation to some variable other than stratum
#observing all the the alleles for which the [] statement is true
baja <- arapat[arapat$Species != "Mainland", ]
plot(baja$EN)

#examing the 01 allele
freq <- frequencies(baja, stratum = "Population", loci = 'EN')
freq.01 <- freq[freqs$Allele == "01", ]
freq.01

#merging coordinate data with frequency data for allele 01 at loci EN
coords <- strata_coordinates(baja)
df <- merge(freq.01, coords)
df[1:10, ]

#plotting allele frequency vs. coordinate Data: could be useful for msat data?
ggplot(df, aes(x = Latitude, y = Frequency)) + geom_line(linetype = 2) + geom_point(size = 2)

#plotting allele frequency vs. coordinate Data but separating based on cluster

#adding a column where the species is 'Baja'
df$Species <- 'Baja'
#looks specifically for the unique populations in the baja dataframe where the cluster is "SCBP-A
pops.with.scbp <- as.character(unique(baja$Population[baja$Cluster == "SCBP-A"]))
#selects the species in the original dataframe, df, and set all that have scbp-A to "Cape"
df$Species[df$Stratum %in% pops.with.scbp] <- "Cape"

#graph of species as determined by the presence or lack thereof of SCBP-A
#determined based on the thought that those in this cluster belong to a different species?
ggplot(df, aes(x = Latitude, y = Frequency)) + geom_line(linetype = 2) + geom_point(size = 5, 
                                                                                    aes(color = Species))

#Spatial Frequency Plots
#this plot shows the frequency of allele 01 at Locus 01 at location
map <- population_map(baja)
ggmap(map) + geom_point(aes(x = Longitude, y = Latitude, size = Frequency), 
                        data = df)

#another way to do it using get_map directly 
centroid2 <- apply(df[, 5:6], 2, mean)
my_map <- get_map(location = centroid2, source = 'google', maptype = 'roadmap', crop = FALSE, zoom = 7) 
ggmap(my_map)+  geom_point(aes(x = Longitude, y = Latitude, size = Frequency), data = df)

#Multivariate Analogs for Loci
to_mv(arapat$WNT[1:10])
to_mv(arapat$WNT[1:10], drop.allele = TRUE)

'TODO: review SVD and principle component analysis'
x <- to_mv(arapat, drop.allele = TRUE)
arapat
x
fit.pca <- princomp(x, cor = TRUE)
#summarizes the proportion of variance in alleles captured by each principle component
#what exactly does fit.pca contain?
summary(fit.pca)

fit.pca
pred <- predict(fit.pca)
pred
#forming a 'feature' table with the principle components and species and population data
df <- data.frame(PC1 = pred[, 1], PC2 = pred[, 2], Species = arapat$Species, 
                 Clade = arapat$Cluster, Pop = arapat$Population)

#graph of the principle components: what variation in species and clade is captured by the first 2 principle components?
ggplot(df) + geom_point(aes(x = PC1, y = PC2, shape = Species, color = Clade), 
                        size = 3, alpha = 0.75)

#examining the peninsula species in particular and looking at the different clades
baja <- pred[df$Species == "Peninsula", ]

'TODO: review clustering tecniques such as agglomerate clustering'
h <- hclust(dist(baja), method = "single")
#is this like a decision tree?
plot(h, main = "Main Baja California Clade", xlab = "")

#http://dyerlab.github.io/popgraph/
#______________________________________________________________________________________________________________
#Creating Population Graphs
#install.packages('popgraph')
require(popgraph)

#can create a popgraph object from scratch
A <- matrix(0, nrow=5, ncol=5)
A[1,2] <- A[2,3] <- A[1,3] <- A[3,4] <- A[4,5] <- 1
A <- A + t(A)
A
g <- as.popgraph( A )
g

#can name the names of the vertexes
library(igraph)
V(g)$name <- c("Olympia","Bellingham","St. Louis","Ames","Richmond")
V(g)$group <- c("West","West", "Central","Central","East")
V(g)$color <- "#cca160"

#can access both the vertex names and the edge connections
list.vertex.attributes( g )
E(g)

E(g)$color <- c("red","red", "red", "blue","dark green")
list.edge.attributes( g )

#accessing the built in data set
data(lopho)
class(lopho)
lopho

#decorate_graph() enables the user to add additional information from an external source
#it does so by combining data from an external source (in this case a dataframe object)
data(baja)
summary(baja)

#adds the additional information from baja to lopo, according to the shared values in Population
lopho <- decorate_graph(lopho, baja, stratum="Population")
lopho

#using igraph routines for visualizing
plot(g)
plot(g, edge.color="black", vertex.label.color="darkred", vertex.color="#cccccc", vertex.label.dist=1)
layout <- layout.circle( g )
plot( g, layout=layout)

layout <- layout.fruchterman.reingold(g)
plot( g, layout=layout)

#plotting using ggplot
#structures the edges & connections, made by looking at genetic variation, on a topographic plane
require(ggplot2)
#creating the edges in the grapj
p <- ggplot() 
p <- p + geom_edgeset( aes(x=Longitude,y=Latitude), lopho) 
p

#adding nodes to the graph
p <- p +  geom_nodeset( aes(x=Longitude, y=Latitude), lopho, size=4)
p

#creating the graph from scratch: edges + nodes colored by region
p <- ggplot() + geom_edgeset( aes(x=Longitude,y=Latitude), lopho, color="darkgrey" )
p <- p + geom_nodeset( aes(x=Longitude, y=Latitude, color=Region, size=size), lopho) 
p <- p + xlab("Longitude") + ylab("Latitude") 
p

#Fruchterman-Reingold algorithm can be used to fid a layout of 'edges'that minimizes the amount of 'energy' so visualization is clear
#https://en.wikipedia.org/wiki/Force-directed_graph_drawing
#easier visualization of the connections: but does not include topographic data?

#applying the Fruchterman-Reingold algorithm to the lopho
c <- layout.fruchterman.reingold( lopho )
#what does this do?: replaces x and y values of lopho
V(lopho)$x <- c[,1]
V(lopho)$y <- c[,2]
p <- ggplot() + geom_edgeset( aes(x,y), lopho, color="darkgrey" )
p <- p + geom_nodeset( aes(x, y, color=Region, size=size), lopho) 
p

#saving a population graph:
save( lopho, file="MyLophoGraph.rda")

#writing the graph is not working
#write.popgraph(lopho,file="~/Desktop/Cactus.pgraph", format="pgraph")

#Spatial Population Graphs
#mapping the nodes and edges onto real space

#Genera & Quick maps
V(g)$Latitude <- c( 47.15, 48.75,38.81, 42.26, 37.74 )
V(g)$Longitude <- c(-122.89,-122.49,-89.98, -93.47, -77.16 )

require(maps)
map( "state" )
#overlap popgraph is not functioning at this time
#Error in overlay_popgraph(g) : Not functioning at this time, I'm in Europe and will fix when I return
#overlay_popgraph(g)

#integrating google and ggplot2 for plotting
require(ggmap)
citation("ggmap")
#this is the centroid: the mean of the longitude and the mean of the latitude
location <- c( mean(V(lopho)$Longitude), mean(V(lopho)$Latitude))
map <- get_map(location,maptype="satellite", zoom=6)
dim(map)
map[1:4,1:4]
#this map object contains raster images? which can be passed onto ggmap() replacing ggplot() for initialization of the graph
p <- ggmap( map ) 
p <- p + geom_edgeset( aes(x=Longitude,y=Latitude), lopho, color="white" ) 
p <- p + geom_nodeset( aes(x=Longitude, y=Latitude, color=Region, size=size), lopho) 
p + xlab("Longitude") + ylab("Latitude")

#Integrating Raster Maps: Raster maps can capture temprature, topography/elevation (might want to use a climate map for rastering?: Ask Guin)
#install.packages('raster')
require(raster)
data(alt)
plot(alt)

#converting spatial nodes and lines into objects that interact with rasters
lopho.nodes <- to_SpatialPoints(lopho)
lopho.nodes
lopho.edges <- to_SpatialLines(lopho)
lopho.edges

plot( alt )
#look into how these paramaters work
plot( lopho.edges, add=TRUE, col="#555555" )
plot( lopho.nodes, add=TRUE, col="black", cex=1.5 )
plot( lopho.nodes, add=TRUE, col=V(lopho)$color, pch=16, cex=1.5 )

#creates a DataFrame with node information: name, lattitude, and longitude: this comes directly from the lopho table
df.nodes <- data.frame(Pop=V(lopho)$name, Latitude=V(lopho)$Latitude, Longitude=V(lopho)$Longitude)
df.nodes
summary(df.nodes)

#rhis creates a 'column' elevation thatextracts the altitide from the raster alt at the points specified by lopho nodes
df.nodes$Elevation <- raster::extract(alt, lopho.nodes )
summary(df.nodes)

#creating a data frame edge with the edges of the lopho and the weights of those edges
df.edge <- data.frame(Weight=E(lopho)$weight )
summary(df.edge)

#also need to extract the 'elevation profile' for the edges 
plot(alt)
#looking at a specific transect
plot(lopho.edges[3],add=TRUE,col="red",lwd=3)
require(ggplot2)

#returns a 'list' of elevations
Elevation <- extract( alt, lopho.edges[3] )[[1]]
#what does extent do?: shows the latitudes and longitudes the transect spans
e <- extent( lopho.edges[3] )
e

#This creates a sequence that cooresponds to the span of latitudes. The sequence goes from the min lattitude to the max and the 
#line cooresponds to the elevation

Latitude <- seq(ymin(e),ymax(e),length.out=length(Elevation))
qplot( Latitude, Elevation, geom="line" )

#Extracting Graph-Theoretic Parameters
'popgraph is based upon igraph: functionalities that can be applied to  igraph and sna packages can be used often for popgraph objects'

to_matrix( g, mode="adjacency")
to_matrix( g, mode="edge weight")
cGD <- to_matrix( lopho, mode="shortest path")
cGD[1:5,1:5]


#can see physical distances according to topographic data using backage fields and function rdist.earth
#install.packages('fields')
require(fields)

pDist <- rdist.earth( cbind( V(lopho)$Longitude, V(lopho)$Latitude ) )

#we can plot the distance data to examine if there is a 
#the matrixes will be symmetric as the matrix shows distance so we extract the upper triangle of the matrix
#the data in the data frame the shortest distance (cGD). the physical distance from the upper triangle are values
df <- data.frame( cGD=cGD[upper.tri(cGD)], Phys=pDist[upper.tri(pDist)])
#what does this value tell us?: tells us the correlation between shortest distance according to genetics? and physicall distance

#https://statistics.laerd.com/statistical-guides/spearmans-rank-order-correlation-statistical-guide.php
'The Spearmans rank-order correlation is the nonparametric version of the Pearson product-moment 
correlation. Spearmans correlation coefficient,(ρ, also signified by rs) measures the strength and
direction of association between two ranked variables'

#this will be very useful graph

cor.test( df$Phys, df$cGD, method="spearman")
qplot( Phys, cGD, geom="point", data=df) + stat_smooth(method="loess") + xlab("Physical Distance") + ylab("Conditional Genetic Distance")

#Node specific paramaters : useful for identifying sink populations
#https://fishbio.com/field-notes/inside-fishbio/source-to-sink#:~:text=Sink%20populations%20exist%20in%20low,source%20population%2C%20would%20become%20extinct.
'sink population: exist in low quality habitat patches: 
not able to support a population in isolation and without contribution of individuals from a source population would be extinct'

#properties of df.nodes: assigning columns in df.nodes to data extracted from lopho
#an analysis of what these values mean can be understood by reading a bit on graph theory
df.nodes$closeness <- closeness(lopho)
df.nodes$closeness
df.nodes$betweenness <- betweenness(lopho)
df.nodes$betweenness
df.nodes$degree <- degree( lopho )
df.nodes$degree
df.nodes$eigenCent <- evcent( lopho )$vector
df.nodes$eigenCent
df.nodes$Region <- factor(V(lopho)$Region)
df.nodes$Region
summary(df.nodes,color="Region")

#install.packages('GGally')
require(GGally)
ggpairs(df.nodes,columns=2:9, aes(color = Region))

df.edge$betweenness <- edge.betweenness(lopho)
df.edge$Region <- rep("Baja",52)
df.edge$Region[36:52] <- "Sonora"
df.edge$Region[c(11,24,27,35)] <- "Cortez"
ggpairs(df.edge, color="Region")

#Edge specific paramaters: assigning the correct regions to the entries
#why do the diagonals come up? Don't they give 0 information
df.edge$betweenness <- edge.betweenness(lopho)
df.edge$Region <- rep("Baja",52)
df.edge$Region[36:52] <- "Sonora"
df.edge$Region[c(11,24,27,35)] <- "Cortez"
ggpairs(df.edge, aes(color = Region))

#Testing for topical congruence
#two different sets of collected data on organisms belonging to the same 2 regions: two sets are the lopho and the upiga
#What exactly is topological congruence? Could this be useful for understanding the topological congruence of the bird species studied in the last paper
data(upiga)
upiga <- decorate_graph(upiga,baja,stratum="Population")
upiga.nodes <- to_SpatialPoints(upiga)
upiga.edges <- to_SpatialLines(upiga)

par(mfrow=c(1,2))
plot(lopho)
plot(upiga)

cong <- congruence_topology(lopho,upiga)
plot(cong)

cong <- decorate_graph( cong, baja )
cong.nodes <- to_SpatialPoints(cong)
cong.edges <- to_SpatialLines(cong)
plot(alt)
plot(cong.edges,add=T)
plot(cong.nodes,add=T, pch=16, col="red")

#are the nodes in the 2 graphs close
test_congruence(lopho,upiga,method="distance")
test_congruence(lopho,upiga, method="combinatorial")



#Questions:
'what are multivariate encodings? This is involved in transforming genotype 
information into information that can be synthesized using popgraphs'

V(lopho)$Longitude


