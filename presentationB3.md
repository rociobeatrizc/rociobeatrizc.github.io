
# Introduction 
Defining where species are found is challenging because biodiversity is complex, and human sampling errors make it even harder. This short project aims to tackle a common sampling mistake by using **virtual species** distributions and **ecological niches**.

The use of empirical data when assessing an error is problematic because each dataset has many confounding factors, which preclude generalisation. A valuable alternative is the simulation of virtual species distributions, because underlying mechanisms that generate such distribution patterns are known. 


A sampling bias in a species often leads to a sample of occurrences that doesn't represent its ecological niche, which is the role and conditions in which that species survives. Comparing the realized ecological niche of the species to that derived from the 'biased' subsample gives us an idea of the discrepancy.

# Input data
## Region of interest
The upload of the area of interest in the form of a **shapefile** is required. It is important to be in the correct folder, within which all the extensions must be present: .shp, .dbf, .shx, .prj.

Our area of interest is the province of Parma, located in Emilia Romagna, in northern Italy

``` R
# Loading required package
library(sf)

# Set working directory 
setwd("C:/chelsa")

# Upload shapefile
aoi <- st_read("parma.shp")
aoi <- aoi$geometry
plot(aoi)
``` 

## Data downloading from CHELSA 

The input data for virtual species is environmental spatial data (**raster data**). 


**CHELSA** (Climatologies at high resolution for the earth’s land surface areas) is a very high resolution (30 arc sec, ~1km) global downscaled climate data set: it is built to provide free access to high resolution climate data for research and application, and is constantly updated and refined.

Bioclimatic variables are derived variables developed for species distribution modeling and related ecological applications: <https://chelsa-climate.org/bioclim>

You can directly download CHELSA data into R. The download creates a series of folders: it is necessary to move to the correct one to work on the downloaded files.


``` R
# Loading required package
library(ClimDatDownloadR) 

# Alternatively, a direct link with GitHub can be created
if(!require(devtools)) install.packages("devtools")
library(devtools)
devtools::install_github("HelgeJentsch/ClimDatDownloadR")

# Download bioclimatic variables from CHELSA 
Chelsa.Clim.download(
  
  # Starting from the workind directory, specify the path 
  save.location = "strade_parma",
  
  # 'bio' contains all the bioclimatic variables
  parameter = "bio",
  
  # Some variables are chosen from the 19 available
  bio.var = c(1, 7, 13, 12, 14),
  
  # Version
  version.var = "2.1",
  
  # Cropping along the area of interest
  clipping = TRUE,
  clip.shapefile = aoi,
  
  # Insert the coordinates of the area of interest (bounding box)
  clip.extent = c(9.439404, 44.34708, 10.50532, 45.04535),
  
  # Buffer, if needed
  # buffer = 3,
  
  # Other commands
  convert.files.to.asc = FALSE,
  stacking.data = TRUE,
  combine.raw.zip = FALSE,
  delete.raw.data = FALSE,
  save.bib.file = TRUE
)
``` 
The required format is a ```RasterStack``` containing all the environmental variables with which you want to generate a virtual species distribution. 
``` R
# Loading required package
library(raster)
library(viridis)

# String containing the names of raster files
rastlist <- list.files(path ="strade_parma/bio/ChelsaV2.1Climatologies/clipped_2024-02-28_09-42-52", pattern = "CHELSA")

# Using the list of names, all the files are imported into a single raster package
mydata <- stack(rastlist)

# Change data names
names(mydata) <- c("mean annual T", "annual range air T", "annual precip", "amount of prec. wettest month", "amount prec. driest month")

# bio1
plot(mydata[[1]], col =  magma(500, alpha = 1, begin = 0, end = 1, direction = 1))
```
<p align="center">
  <img src="magma_t.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>bio1 is the mean annual air temperature</em>
</p>


### Correlation Matrix 
In studying species distribution based on climate variables, conducting a correlation analysis is crucial. 

This analysis helps identify which climate variables are strongly associated with each other and with species distribution patterns. Avoiding highly correlated variables is important because they provide redundant information, potentially leading to multicollinearity issues in statistical models. 

By focusing on uncorrelated variables, models become more interpretable and accurately capture the relationships between climate and species distribution. This enhances predictive accuracy and conserves computational resources by excluding unnecessary variables.
``` R
# Loading required package
library(corrplot)

# Subsample 10% of pixels and calculate pairwise correlations
r1 <- mydata$mean.annual.T
cor <- cor(sampleRandom(mydata, size= ncell(r1) * 0.30 ), method = "pearson")

# Plot correlation matrix
df <- corrplot(cor, method = "number")
```
<p align="center">
  <img src="corrplot.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>Correlation matrix</em>
</p>

# Random species generation
Generating random species from known environmental data allows controlling the factors that can influence the distribution of real data. To create a series of occurrence points for a species, it is necessary to go through 3 steps
## First step: suitability
By intersecting bioclimatic data, the first output obtained is a suitability map. In creating this map, various parameters can be set, allowing for closer approximation to a real-world scenario. 
- ```species.type```. There are two main possibilities to combine response functions to calculate the environmental suitability:
**additive** (the response functions are added) and **multiplicative** (the response functions are multiplied, default).
- ```approach``` By default, the function ```generateRandomSp``` uses a **PCA** approach if you have 6 or more variables, and a ‘**response functions**’ approach if you have less than 6 variables. 
- ```relations```.  Four relations are possible for a random generation of virtual species: **gaussian**, **linear**, **logistic** and **quadratic** relations. By default, all the relation types are used. You can choose to use any combination inside the argument.
``` R
# Loading required packages
library(virtualspecies)
library(ggplot2)
library(tidyverse)

# Suitability map generation
random.sp <- generateRandomSp(raster.stack = mydata,
                              convert.to.PA = FALSE,
                              
                              # How to combine response functions
                              species.type = "multiplicative",

                              # Random approach between PCA and response function
                              approach = "random",
                              
                              # Response function
                              relations = "gaussian",

                              # Realistic species
                              realistic.sp = TRUE,
                              plot = FALSE)

plot(random.sp$suitab.raster, col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))                             
``` 
<p align="center">
  <img src="suitability_markdown.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>Suitability map</em>
</p>

You can visualize how each variable contributed to generating the suitability map using the ```plotResponse``` command

``` R
plotResponse(random.sp)
``` 
<p align="center">
  <img src="plotresponse.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>Response functions</em>
</p>


## Second step: presence/absence
In the second step, the suitability map is converted into a binary **presence/absence map** through a probability function (logisatic curve) that associates the suitability value with the probability of finding the virtual species for each pixel. 

This subset of the environmental
niche that is actually occupied by the species corresponds to the
**realized niche** (*Hutchinson, 1957*). 

Customizing the function is also possible in this step by the parameters ```alpha``` and ```beta``` which determine the shape of the logistic curve. 
- ```alpha``` drives the ‘slope’ of the curve. 
- ```beta``` controls the inflexion point. A higher beta will decrease the probability of finding suitable conditions (smaller distribution range).
- ```species.prevalence```. The species prevalence is the number of places (here, pixels) actually occupied by the species out of the total number of places (pixels) available.
``` R
# Presence/Absence: requires defining the parameters alpha, beta, and species prevalence
new.pres <-convertToPA(random.sp,
                       beta = "random",
                       alpha = -0.05, plot = FALSE, 
                       species.prevalence = 0.1)

plot(new.pres$pa.raster)
``` 
<p align="center">
  <img src="realized_niche.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>Realized niche</em>
</p>

## Third step: occurrences
The third step consists of generating, within the presence/absence raster, a series of occurrence points. This step involves many parameters that are important to test.
- ```n```: number of occurrences
- ```type```. The type of occurrences can be 'presence-absence' or 'presence only'. In the case of choosing the second type, points will be generated within the realized niche, so the maximum number of points allowed is equal to the quantity of pixels present in the presence-absence raster with a value of 1. You can easily find this number (ncell) with the following command:
``` R
 > new.pres$pa.raster %>% raster()
class      : RasterLayer 
dimensions : 83, 128, 10624  (nrow, ncol, ncell)
...
```
- ```sample.prevalence```: the sample prevalence is the proportion of samples in which the species has been found.
- ```error.probability```: the probability of an erroneous detection of our virtual species. It's a human error and could depend, for example, on the sampler's ability to correctly recognize a species: the error is common in the case of plants. The error probability is useless in a ‘presence only’ sampling scheme, because the sampling strictly occurs within the boundaries of the species distribution range. Nevertheless, if you still want to introduce errors, then make a sampling scheme ‘presence-absence’, and use only the ‘presence’ points obtained. 
- ```detection.probability```: the detection probability of our virtual species, ranging from 0 (species cannot be detected) to 1 (species is always detected). This factor, instead, depends on the species: there are species more likely to be sampled.
- ```correct.by.suitability``` a parameter you can set to complexify the detection probability by making it dependent on the environmental suitability. In this case, cells will be weighted by the environmental suitability: less suitable cells will have a lesser chance of detecting the species. 
``` R
  presence.points <- sampleOccurrences(new.pres,
                                       n = 50, 
                                       type = "presence only",
                                       sample.prevalence = 0.9,
                                       error.probability = 0,
                                       detection.probability = 1,
                                       correct.by.suitability = TRUE,
                                       plot = FALSE)

# Plot
plot(mydata[[1]], col =  magma(500, alpha = 1, begin = 0, end = 1, direction = 1))
points(presence.points$sample.points, col = "black", pch = 19)
  ``` 

<p align="center">
  <img src="slippery.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>Occurrences of the random species</em>
</p>

# Preliminary Steps for Niche Analysis
The following are a series of necessary steps to associate the corresponding bioclimatic variables with the realized niche (step 2) and occurrence points (step 3). 
``` R
# Raster with presence/absence points (step 2): used to create the realized niche 
# It must be a RasterLayer object
raster01 <- new.pres$pa.raster %>% raster()

# The bioclimatic layers are extracted one by one from the stack
r1 <- mydata$mean.annual.T
r2 <- mydata$annual.range.air.T
r3 <- mydata$annual.precip
r4 <- mydata$amount.of.prec..wettest.month
r5 <- mydata$amount.prec..driest.month

# A stack is created containing the bioclimatic variables and the raster of presence/absence (realized niche)
stack_pa <- brick(r1, r2, r3, r4, r5, raster01)

# The values are extracted from the stack, converted into a dataframe, and only the presence pixels (1) are retained
values <- stack_pa %>% rasterToPoints() %>% as.data.frame()
filtered_pa <- values %>% filter(., lyr.1 == 1) %>% as.data.frame()
```
This way, ```filtered_pa``` is obtained, a dataframe containing the points of the realized niche and their corresponding bioclimatic variables. 


``` R
# Conversion from a dataframe to a raster
raster_pa <- filtered_pa %>% .[,-8] %>% rasterFromXYZ()
plot(raster_pa)
```
<p align="center">
  <img src="raster_pa.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>Bioclimatic variables only corresponding to the realized niche</em>
</p>

The same steps will be performed for the occurrence points of the species
```R
# The raster of occurrences is transformed into a dataset, from which the rows satisfying both conditions Real = 1 and Observed = 1 are preserved
raster_occurences <- presence.points$sample.points %>% as.data.frame() %>% .[.$Real == 1 & .$Observed == 1, ]

# The environmental variables are associated with the occurrences using their coordinates
stack_occ <- brick(r1, r2, r3, r4, r5)
values_occ <- stack_occ %>%  rasterToPoints() %>% as.data.frame()
filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))
```



```filtered_occ``` contains the coordinates of occurrence points and their corresponding environmental variables


# Introducing roadside bias
One of the most recognized forms of bias in distributional data is the high concentration of observations (or collection sites) along road: models based on data collected along roads should lead to inaccurate predictions when applied to larger area (*Kadmon et al., 2004*)

The main mechanism by which roadside bias can affect the accuracy of bioclimatic models is through climatic bias in the distribution of the road network.

Our goal is to associate the occurrence points with the probability of being sampled, knowing that as one moves away from roads, the probability decreases: the first step involves calculating the distance of the points from the nearest road. In our case, the road network is that of the province of Parma.
``` R
# Loading required package
library(terra)

# The occurrences (data.frame) are transformed into a SpatVector object
coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")

# For calculating distances in meters, it is necessary to change the CRS from 4326 to 3857. Keeping WGS84 will result in distances in degrees
coord_occ_dist <- terra::project(coord_occ, "EPSG:3857")

# Upload the shapefile containing the road network and convert it into a SpatVector object
setwd("C:/chelsa")
roads <- "roads_parma.shp" %>% st_read() %>%  .$geometry %>% terra::vect() 

# Here, the CRS also needs to be changed
roads <- terra::project(roads, "EPSG:3857")

# Visualization
plot(roads)
points(coord_occ_dist)
```
<p align="center">
  <img src="roads_points.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>Occurrence points in relation to the road network</em>
</p>

```R
# With this function, the distance of each point from all roads is obtained
dist <- distance(coord_occ_dist, roads, unit ="m")

# Distance of each point from the nearest road
min_dist <- apply(dist, 1, min)
```
With a cumulative density function, you can observe the distribution of points relative to roads: it saturates quickly
``` R
# The distances are extracted and associated with the dataset of occurrences
coord_occ$distance <- min_dist

# To create the function, distances are extracted from the dataset of occurrences in the form of a dataframe
occ_hist <- coord_occ %>% .$distance %>% as.data.frame()


# Empirical cumulative distribution function
ecdf_fun <- ecdf(occ_hist$occ_hist)

# Plot
plot(ecdf_fun, verticals = TRUE,
     main = "Cumulative Frequency",
     xlab = "Distance (m)", 
     ylab = "Relative Frequency",
     xlim = c(0, 33000), pch = 20, col = "dodgerblue3")

```

<p align="center">
  <img src="cumulative.png" alt="Image Description" style="display: block; margin: auto;" />
</p>
<p align="center">
  <em>Cumulative frequency</em>
</p>

## Sampling probability
Dato che la probabilità di campionare diminuisce con la distanza dalla strada, possiamo immaginare una curva logaritmica decrescente.
Il parametro c opera una trasformazione sulle distanze, modulando la pendenza della curva. Valori prossimi allo 0 di c abbassano il valore della distanza, associando anche a grandi distanze una buona probabilità di trovare campioni.
``` R
########## Probabilità di campionamento
# c <- "sampling effort": 1-(((log(c*distanza_minima))/(log(c*max(distanza_minima)))))
prob_campionamento <- 1-(((log(min_dist))/(log(max(min_dist)))))

#### Funzione sperimentale
# options(scipen=999)
# prob_campionamento <- (0.9999/((max(distanza_minima))^0.9999))*((1/distanza_minima)^(1-0.9999))
hist_prob <- bind_cols(prob_camp = prob_campionamento, distance = min_dist) %>% 
  ggplot(aes(x = distance, y = prob_camp)) +
  ylim(0, 1) +
  geom_point()

hist_prob
```

Comunque venga costruita, la curva significa che non tutti i punti verranno campionati allo stesso modo: vicino alle strade sono campionati di più, mentre man mano che ci si allontana, con la diminuzione della probabilità, diminuisce anche il numero effettivo di punti raccolti. 
Vediamo innanzitutto come si distribuiscono i punti rispetto alla distanza
```R
 # In breaks introduco il numero di intervalli
# Secondo quale criterio? Intervalli uguali, ma devo definirli a mano
coord_occ %>% as.data.frame() %>% ggplot(., aes(x = distance)) + geom_histogram(binwidth = 500)

summary(min_dist)
```
Risulta utile stabilire delle classi di distanza, suddividendo le distanze in intervalli: ogni classe avrà il suo numero di punti, la sua probabilità media di campionamento e la quantità di punti che si ottengono se si considera il roadside bias. 
```R
#Calcola i limiti degli intervalli in base alle statistiche
lim <- c(0, 500, 1000, 5000, 10000, 15000, 25000, 57000)  

# Definisci le nuove etichette per gli intervalli
# nuove_etichette <- c("0-50", "50-100", "100-300", "300-500", "500-1000", "1000-1500", "1500-3000", "3000-5000", "5000-10000", "10000-15000", "15000-25000", "25000-33000")
?bind_cols
distance_class <- bind_cols(prob_camp = prob_campionamento,  coord_occ %>% as.data.frame(geom="XY"))   %>% 
  add_column(class = cut(.$distance, breaks = lim, labels = F))


# Punti per classe
count(distance_class, class)
nrow(distance_class)

groups <- distance_class %>% 
  group_by(class) %>% 
  summarise(plotn = n(), meanprob = mean(prob_camp)) %>% 
  add_column(plot_sampl = .$plotn * .$meanprob)

groups


whole_data <- merge(distance_class, groups, by="class")
whole_data
nrow(whole_data) 

```


## Dataset resampled based on probability
Ricampionando il dataset originale attraverso la probabilità, che opera da filtro, si ottiene un dataset ridotto. 

``` R

### Punti filtrati per probabilità
by_prob <- whole_data %>% 
  group_by(class) %>% 
  group_split() %>% 
  map(function(z){
    set.seed(1234)
    z %>% 
      sample_frac(size = mean(z$prob_camp))
  }) %>% 
  do.call(bind_rows, .)


count(by_prob, class)
view(by_prob)
by_prob
nrow(by_prob)

plot(strade)
points(whole_data %>% terra::vect(., geom = c("x","y"), crs="epsg:4326") %>% 
         terra::project(., "EPSG:3857"), col="blue")
points(by_prob %>% terra::vect(., geom = c("x","y"), crs="epsg:4326") %>% 
         terra::project(., "EPSG:3857"), col="red")

########## Filtra solo le righe con la classe 3
points_prob <- by_prob %>% filter(class==1) %>% terra::vect(., geom = c("x","y"), crs="epsg:4326")
# Per proiettarli vicino alle strade: %>%  terra::project(., "EPSG:3857"), col="red")
nrow(points_prob)
#################################################
# Estrai le coordinate x e y
coordinate_prob <- points_prob %>%  st_as_sf() %>% sf::st_coordinates()
coordinate_prob
# Crea un dataframe con le coordinate x e y
coordinate_prob <- data.frame(x = coordinate_prob[, "X"], y = coordinate_prob[, "Y"])
print(coordinate_prob)

```


# Niche analysis 

Kadmon
Since data on species distribution are often biased
and incomplete, distribution maps derived directly
from such data may reflect patterns of sampling efforts
rather than true patterns of species distribution. One
possible approach to overcome this problem is to apply
"ecological niche models" (sensu Peterson et al. 2002)
as a tool for mapping distribution ranges. According
to this idea, the data available on the distribution of
the relevant taxon are used to identify its distribution
in an ecological, rather than geographical space. 

How much does resource use by species A overlap with
that of species B? Recent concern over the effects of global
change on species distributions has emphasized the need to
quantify differences among species in their environmental
requirements in a geographical context and at an extent comparable
to that of species ranges.


Measurement of niche overlap
The comparison of zij between two entities can be used to
calculate niche overlap using the D metric (Schoener, 1970;
reviewed in Warren et al., 2008).
This metric varies between 0 (no overlap) and 1 (complete
overlap)

This is possible in R through ecospat
The aim of the ecospat package is to make available novel tools and methods to support spatial analyses and modeling
of species niches and distributions in a coherent workfl ow. Th e package is written in the R language (R Development Core
Team) and contains several features, unique in their implementation, that are complementary to other existing R packages.
Pre-modeling analyses include species niche overlap, quantifi cations and comparisons between distinct ranges or time periods,
measures of phylogenetic diversity, and other data exploration functionalities. 

Riprendiamo il dataset che contiene le coordinate dei punti campionati secondo la probabilità: scegliamo, ad esempio, di visualizzare la nicchia ecologica generata dai punti che appartengono alla classe di distanza 3, da tot metri a tot metri: quanto la nicchia di questo sottocampione si discosta dalla nicchia realizzata?
A questi punti bisogna associare le variabili bioclimatiche che si trovano in corrispondenza.

``` R
########## Filtra solo le righe con la classe 3
points_prob <- by_prob %>% filter(class==1) %>% terra::vect(., geom = c("x","y"), crs="epsg:4326")
# Per proiettarli vicino alle strade: %>%  terra::project(., "EPSG:3857"), col="red")
nrow(points_prob)
#################################################
# Estrai le coordinate x e y
coordinate_prob <- points_prob %>%  st_as_sf() %>% sf::st_coordinates()
coordinate_prob
# Crea un dataframe con le coordinate x e y
coordinate_prob <- data.frame(x = coordinate_prob[, "X"], y = coordinate_prob[, "Y"])
print(coordinate_prob)
###############   Nicchia
# Uniamo i due dataframe (occorrenze originali e punti a meno di 1 km) basandoci sulle colonne "x" e "y" 
new_points <- merge(filtered_occ, coordinate_prob, by = c("x", "y"))

new_points


## Passo da dataframe a raster
new_points <- new_points[,-(8:9)]
new_points
raster_occ <- rasterFromXYZ(new_points)
```
For the niche quantification, a Principal Component Analysis (PCA) of the environmental data is carried out.

``` R

### A questo punto ho le occorrenze e le variabili climatiche ad esse legate 
### ecospat niche overlap
### https://rstudio-pubs-static.s3.amazonaws.com/1044769_cfe9b0a3704844309f1884a38bc037ed.html
### Measuring ecological niche overlap from occurrence and spatial environmental

library(ecospat)
library(ade4)
# For the niche quantification, a matrix with the background environmental variables from both ranges, 
# as well as the global environment are needed. 
# After cropping the rasters, getValuesis applied to convert them to a data frame.

# Dati bioclimatici dello step 2 sottoforma di dataframe
env_pa <- getValues(raster_pa)
env_pa

# Dati bioclimatici dello step 3 sottoforma di dataframe
env_occ <- getValues(raster_occ)
env_occ

# remove missing values
env_occ <- env_occ[complete.cases(env_occ), ]
env_occ <- as.data.frame(env_occ)

env_pa <- env_pa[complete.cases(env_pa), ]
env_pa

# produce global environmental background data
globalEnvM <- rbind(env_pa, env_occ)
globalEnvM

## PCA
pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)

global.scores <- pca.clim$li

pa.scores <-
  suprow(pca.clim,
         data.frame(filtered_pa)[, colnames(globalEnvM)])$li   

occ.scores <-
  suprow(pca.clim,
         data.frame(filtered_occ)[, colnames(globalEnvM)])$li

pa.scores1 <- suprow(pca.clim, env_pa)$li
occ.scores1<- suprow(pca.clim, env_occ)$li

data.frame(filtered_pa)[, colnames(globalEnvM)]


# calculate the Occurrence Density Grid for both native and invasive species
pagrid <- ecospat.grid.clim.dyn(global.scores,
                                pa.scores,
                                pa.scores1)

occgrid <- ecospat.grid.clim.dyn(global.scores,
                                 occ.scores, 
                                 occ.scores1)
dev.off()
ecospat.plot.niche.dyn(pagrid, occgrid, quant = 0.1, interest = 2, name.axis1 = "PC1", name.axis2 = "PC2")
?ecospat.plot.niche

# plot variable contributions
ecospat.plot.contrib(contrib=pca.clim$co, eigen=pca.clim$eig)

# dynamics index
ecospat.niche.dyn.index(pagrid, occgrid, intersection = 0.1)$dynamic.index.w

# calculate niche overlap
ecospat.niche.overlap(pagrid, occgrid, cor=T)

# perform the Niche Equivalency Test
#eq.test <- ecospat.niche.equivalency.test(pagrid, occgrid, rep = 100, ncores = 2)
# perform the Niche Similarity Test
# sim.test <- ecospat.niche.similarity.test(pagrid, occgrid, rep = 100, rand.type = 2, ncores = 2)
# plot Equivalency and Similarity Test
# par(mfrow=c(1,2))
# ecospat.plot.overlap.test(eq.test, "D", "Equivalency") 
# ecospat.plot.overlap.test(sim.test, "D", "Similarity")
```

# Cube (alla fine)
Aggiungendo una colonna contenente il tempo, ovvero un anno casuale all'interno dell'intervallo climatico (1980-2010) si ottiene un cubo di dati: specie, x, y, distanza. Quest'ultimo fattore, come abbiamo visto, influenza il campionamento.

Il codice finora genera soltanton una specie: per ottenere un dataframe contenente un numero a scelta di specie casuali, occorre iterare il procedimento come segue.
``` R
# Definisci il numero di iterazioni
num_iterazioni <- 3

# Crea una lista vuota per memorizzare i dataframe
dataframe_list_prima <- list()
dataframe_list_dopo <- list()

# Ciclo for per generare specie e calcolare le distanze per ciascuna iterazione
for (i in 1:num_iterazioni) {
  set.seed(i)
  random.sp <- generateRandomSp(raster.stack = mydata,
                                convert.to.PA = FALSE,
                                species.type = "multiplicative",
                                approach = "random",
                                relations = c("gaussian", "logistic", "quadratic"),
                                realistic.sp = TRUE,
                                plot = FALSE)
  random_beta <- runif(1)
  new.pres <- convertToPA(random.sp,
                          beta = random_beta,
                          alpha = -0.05, plot = FALSE)
  presence.points <- sampleOccurrences(new.pres,
                                       n = 50, 
                                       type = "presence only",
                                       sample.prevalence = 0.9,
                                       error.probability = 0,
                                       detection.probability = 1,
                                       correct.by.suitability = TRUE,
                                       plot = FALSE)
  raster01 <- new.pres$pa.raster
  raster01 <- raster01 %>% raster()
  
  
  r1 <- mydata$mean.annual.T
  r2 <- mydata$annual.range.air.T
  r3 <- mydata$annual.precip
  r4 <- mydata$amount.of.prec..wettest.month
  r5 <- mydata$amount.prec..driest.month
  
  stack_pa <- brick(r1, r2, r3, r4, r5, raster01)
  values <- rasterToPoints(stack_pa)
  values <- as.data.frame(values)
  
  filtered_pa <- filter(values, lyr.1 == 1)
  filtered_pa <- as.data.frame(filtered_pa)
  filtered_pa <- filtered_pa[,-8]
  
  raster_pa <- rasterFromXYZ(filtered_pa)
  raster_occurences <- presence.points$sample.points
  raster_occurences <- as.data.frame(raster_occurences)
  
  filtered_occ <- raster_occurences[raster_occurences$Real == 1 & raster_occurences$Observed == 1, ]
  filtered_occ <- filtered_occ[,-(3:4)]
  
  stack_occ <- brick(r1, r2, r3, r4, r5)
  
  values_occ <- rasterToPoints(stack_occ)
  values_occ <- as.data.frame(values_occ)
  
  filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))
  
  coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")
  coord_occ_dist <- terra::project(coord_occ, "EPSG:3857")
  
  setwd("C:/chelsa")
  aoi <- st_read("parma.shp")
  aoi <- aoi$geometry

  strade <- st_read("strade_parma.shp")
  strade <- strade$geometry
  
  strade <- terra::vect(strade)
  
  aoi <- terra::vect(aoi)
  aoi <- terra::project(aoi,"EPSG:3857")
  strade <- terra::project(strade, "EPSG:3857")
  dist <- distance(coord_occ_dist, strade, unit ="m")
  
  distanza_minima <- apply(dist, 1, min)
  coord_occ$distanze <- distanza_minima
  
  prob_campionamento <- 1-(((log(distanza_minima))/(log(max(distanza_minima)))))
  distance_class <- bind_cols(prob_camp = prob_campionamento,  coord_occ %>% as.data.frame(geom="XY"))   %>% 
    add_column(class = cut(.$distanze, breaks = lim, labels = F))
  
  groups <- distance_class %>% 
    group_by(class) %>% 
    summarise(plotn = n(), meanprob = mean(prob_camp)) %>% 
    add_column(plot_sampl = .$plotn * .$meanprob)
  
  whole_data <- merge(distance_class, groups, by="class")
  set.seed(i)
  whole_data$tempo <- sample(1980:2010, nrow(whole_data), replace = TRUE)
  whole_data
  
  by_prob <- whole_data %>% 
    group_by(class) %>% 
    group_split() %>% 
    map(function(z){
      set.seed(1234)
      z %>% 
        sample_frac(size = mean(z$prob_camp))
    }) %>% 
    do.call(bind_rows, .)
  
  coord_occ <- as.data.frame(whole_data)
  dataframe_prima <- whole_data[,-(1:9)]
  dataframe_prima <- dataframe_prima[,-(4:6)]
  
  # Aggiungi una colonna "specie" con il valore dell'iterazione corrente
  dataframe_prima$specie <- paste("specie", i)
  
  # Aggiungi il dataframe ottenuto alla lista
  dataframe_list[[i]] <- dataframe_prima
  by_prob <- as.data.frame(by_prob)
  by_prob
  dataframe_dopo <- by_prob[,-(1:9)]
  dataframe_dopo <- dataframe_dopo[,-(4:6)]
  
  # Aggiungi una colonna "specie" con il valore dell'iterazione corrente
  dataframe_prima$specie <- paste("specie", i)
  dataframe_dopo$specie <- paste("specie", i)
  
  
  # Aggiungi il dataframe ottenuto alla lista
  dataframe_list_prima[[i]] <- dataframe_prima
  dataframe_list_dopo[[i]] <- dataframe_dopo
  
}




# Combinare tutti i dataframe ottenuti in un unico dataframe
dataframe_before <- do.call(rbind, dataframe_list_prima)
dataframe_after <- do.call(rbind, dataframe_list_dopo)

dataframe_before <- dataframe_before[, c("specie", "x", "y", "tempo", "distanze")]
dataframe_after <- dataframe_after[, c("specie", "x", "y", "tempo", "distanze")]
```

