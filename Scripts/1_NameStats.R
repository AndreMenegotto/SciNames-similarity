
# ------------------------------------------------------------
# Script:   1_NameStats.R
# Purpose:  Check similarity among valid names of marine species
# Updated:  15-06-2026
# ------------------------------------------------------------


#. Load packages ----
library(stringdist)
library(dplyr)



#. Loading and preparing the data ----

## Read the data set
# To access WoRMS data it is necessary to apply a request at 'https://www.marinespecies.org/usersrequest.php'
patch <- "./Data/raw/WoRMS_download_2025-04-01/"
taxon <- data.table::fread(paste(patch,"taxon.txt",sep=""), na.strings=c("","NA"))
taxon <- as.data.frame(taxon)


## Select only valid species names
species <- taxon %>% filter(taxonRank %in% "Species" & taxonomicStatus %in% "accepted")
species <- species %>% select(scientificName, phylum) %>% distinct(scientificName, .keep_all = T)
any(duplicated(species$scientificName)) #FALSE

# Remove hybrid or unidentified species
rmSP <- grep("\\[|\\]|[0-9]|×", species$scientificName)
species <- species[-rmSP,]
species$scientificName <- gsub("['\"]", "", species$scientificName)

# Sort names and convert to lower case
species <- species[order(species$scientificName),]
species$scientificName <- tolower(species$scientificName)

speciesList <- strsplit(species$scientificName, " ")
rmSP <- which(sapply(speciesList, function(x) any(nchar(x) == 1)))
speciesList <- speciesList[-rmSP]
species <- species[-rmSP,]

table(lengths(speciesList))
# 3 = subgenus
# 4 = incertae sedis
# >4 = virus description

pos <- which(lengths(speciesList) < 4)
speciesList <- speciesList[pos]
length(speciesList) #229,578

species <- species[pos,]
length(table(species$phylum)) #85

# Select only genus and epithet (ignore subgenus)
genus <- sapply(speciesList,"[[",1)
epithet <- sapply(speciesList,tail,1L)

tmpSpecies <- paste(genus, epithet, sep=" ")
length(tmpSpecies) #229,578
length(unique(tmpSpecies)) #229,508 (removing subgenus created duplicates)

# Select unique names
pRM <- which(!duplicated(tmpSpecies))

genus <- genus[pRM]
epithet <- epithet[pRM]
tmpSpecies <- tmpSpecies[pRM]
species <- species[pRM,]
speciesList <- strsplit(tmpSpecies, " ")
#save(genus, epithet, tmpSpecies, species, speciesList, file = "./Data/tmpRData/CheckedNames.RData")
#load("./Data/tmpRData/CheckedNames.RData")



#. Comparing similarity (general names) ----

## Get the lowest number of edits among all valid names

# Output to save the lowest number of edits and the closest name
nearNames <- data.frame(vec=tmpSpecies, closest=NA, min_dist=NA)
pb <- txtProgressBar(min = 0, max = nrow(nearNames), style = 3)

# Loop
vec <- nearNames$vec
for(i in 1:nrow(nearNames))
{
  allOthers <- vec[-i]
  
  # Calculate the distance between target and all other names
  dis <- stringdist(vec[i], allOthers, method = "dl")
  posD <- which.min(dis)
  
  # Add the information in the output data.frame
  nearNames$closest[i] <- allOthers[posD]
  nearNames$min_dist[i] <- dis[posD]
  
  setTxtProgressBar(pb, i)
}
#save(nearNames, file = "./Data/tmpRData/nearNames.RData")
#load("./Data/tmpRData/nearNames.RData")


## Check results
sum(nearNames$min_dist==1) #1,026
sum(nearNames$min_dist==2) #8,063
sum(nearNames$min_dist==3) #28,564
round(c(1026,8063,28564)/nrow(nearNames)*100, 2) #0.45 3.51 12.45


## Histogram (Fig. 1A)
barplot(table(nearNames$min_dist)/nrow(nearNames), space=0, col="grey80", ylim=c(0,.25), xlim=c(0,19), las=2, axisnames=F, cex.axis=1.25)
axis(side=1, at=seq(1,19,2)-.5, labels=NA, tck=-.035)
axis(side=1, at=seq(2,18,2)-.5, labels=NA, tck=-.02, lwd=0, lwd.ticks = 1)
axis(side=1, at=seq(1,19,2)-.5, labels=seq(1,19,2), tick=F, cex.axis=1.25)


## Checking for multiple matches for those at DL <= 4

# Get only the names at 4 or less DL distance from other names 
x <- nearNames$vec[nearNames$min_dist<5]

# Output to save the lowest number of edits and all the closest names
nearNamesAll4 <- data.frame(vec=character(), closest=character(), min_dist=integer())
pb <- txtProgressBar(min = 0, max = length(x), style = 3)

# Loop
vec <- x
for(i in 1:length(x))
{
  allOthers <- vec[-i]
  
  # Calculate the distance between target and all other names
  dis <- stringdist(vec[i], allOthers, method = "dl")
  posD <- which.min(dis)
  posD <- which(dis == dis[posD])
  
  # Add the information in the output data.frame
  tmp <- data.frame(vec=vec[i], closest=allOthers[posD], min_dist=dis[posD])
  nearNamesAll4 <- rbind(nearNamesAll4, tmp)
  
  setTxtProgressBar(pb, i)
}
#save(nearNamesAll4, file = "./Data/tmpRData/nearNamesAll4.RData")
#load("./Data/tmpRData/nearNamesAll4.RData")

round(table(table(nearNamesAll4$vec[nearNamesAll4$min_dist==1]))/length(unique(nearNamesAll4$vec[nearNamesAll4$min_dist==1]))*100, 2)
round(table(table(nearNamesAll4$vec[nearNamesAll4$min_dist==2]))/length(unique(nearNamesAll4$vec[nearNamesAll4$min_dist==2]))*100, 2)
round(table(table(nearNamesAll4$vec[nearNamesAll4$min_dist==3]))/length(unique(nearNamesAll4$vec[nearNamesAll4$min_dist==3]))*100, 2)
round(table(table(nearNamesAll4$vec[nearNamesAll4$min_dist==4]))/length(unique(nearNamesAll4$vec[nearNamesAll4$min_dist==4]))*100, 2)


## Get names with three edits or less
edt3 <- nearNamesAll4[nearNamesAll4$min_dist<4,]

speciesList1 <- strsplit(edt3$vec, " ")
speciesList2 <- strsplit(edt3$closest, " ")

# Split genus and epithet
genus1 <- sapply(speciesList1,"[[",1)
epithet1 <- sapply(speciesList1,tail,1L)
genus2 <- sapply(speciesList2,"[[",1)
epithet2 <- sapply(speciesList2,tail,1L)

# Get the number of edits for each component
dlGE <- stringdist(genus1, genus2, method = "dl")
dlEP <- stringdist(epithet1, epithet2, method = "dl")

# Check where the edits are concentrated
#View(cbind(edt3, dlGE, dlEP))
round(sum(dlGE==0 | dlEP==0)/nrow(edt3)*100, 2) #99.20%
round(sum(dlGE==0)/nrow(edt3)*100, 2) #81.17%
round(sum(dlEP==0)/nrow(edt3)*100, 2) #18.03%



#. Applying post-matching filters ----
source("./Scripts/0_postFilters.R")


## Post-Filter from Taxatmatch
resFiltered <- numeric()
for(i in 1:4)
{
  x <- nearNamesAll4[nearNamesAll4$min_dist==i,]
  
  genus1 <- sapply(strsplit(x$vec, " "),"[[",1)
  epithet1 <- sapply(strsplit(x$vec, " "),tail,1L)
  
  genus2 <- sapply(strsplit(x$closest, " "),"[[",1)
  epithet2 <- sapply(strsplit(x$closest, " "),tail,1L)
  
  dlGe <- mapply(stringdist, genus1, genus2, method = "dl")
  dlEp <- mapply(stringdist, epithet1, epithet2, method = "dl")
  
  x2 <- data.frame(genus1, epithet1, genus2, epithet2, dlGe, dlEp)
  
  postGe <- mapply(genus_postFilter, genus=x2$genus1, preGenera=x2$genus2, GeneraDist=x2$dlGe)
  x2 <- x2[postGe,]
  
  postSP <- mapply(epithet_postFilter, epithet = x2$epithet1, preEpithet = x2$epithet2, EpithetDist = x2$dlEp, GeneraDist = x2$dlGe)
  x2 <- x2[postSP,]
  
  resFiltered <- rbind(resFiltered, x2)
  
  print(round(nrow(x2)/nrow(x)*100, 2))
}


## How much remain after the filters
zz <- data.frame(sciname=paste(resFiltered$genus1,resFiltered$epithet1), DL=rowSums(resFiltered[, c("dlGe","dlEp")]))
zz <- zz %>% distinct(sciname, .keep_all = T)

table(zz$DL) #1026 4070 9060 2622
round(table(zz$DL)/nrow(nearNames)*100, 2) #0.45 1.77 3.95 1.14

# How much was effectively removed
round(table(zz$DL)/table(nearNames$min_dist[nearNames$min_dist<5])*100, 2)
#100.00  50.48  31.72   4.92


## Check where the edits are concentrated in the post-filter remaining matches
aa <- rowSums(resFiltered[, c("dlGe","dlEp")])

gg <- resFiltered[aa<4,"dlGe"]
round(sum(aa[aa<4]==gg)/sum(aa<4)*100, 2) #21.16

ee <- resFiltered[aa<4,"dlEp"]
round(sum(aa[aa<4]==ee)/sum(aa<4)*100, 2) #77.41
21.16+77.41 #98.57%



#. Comparing similarity (when shared elements) ----

## same genera (how distant are the epithets when names share the genus?)
dG <- table(genus)
posDupG <- which(dG>1) # which genera have more than one species...

length(dG) #33,926 genera
length(posDupG) #20,846 non-exclusive
sum(dG==1) #13,080

# Get all the names with non unique genera
tmpSpecies <- speciesList[which(genus %in% names(dG[posDupG]))]
tmpSpecies <- as.data.frame(do.call(rbind, tmpSpecies), stringsAsFactors = FALSE)
dim(tmpSpecies) #216,428
round(nrow(tmpSpecies)/length(speciesList)*100, 2) #94.3%

# For each genus
dupGen <- names(dG[posDupG])
minDLep <- numeric()
for(i in 1:length(dupGen)) #20,846
{
  # Get the minimum distance among its epithets
  pos <- which(tmpSpecies$V1 == dupGen[i])
  epi <- tmpSpecies$V2[pos]
  
  dlEp <- stringdistmatrix(epi, epi, method = "dl")
  dlEp[upper.tri(dlEp, diag = TRUE)] <- NA
  dlEp <- na.omit(as.vector(dlEp))
  minDLep[i] <- min(dlEp)
}


## same epithet (how distant are the genera when names share epithet?)
dE <- table(epithet)
posDupE <- which(dE>1) # which epithet is used in more than one species...

length(dE) #85,863 epithets
length(posDupE) #24,655 non-exclusive
sum(dE==1) #61,208

# Most common epithets in WoRMS
sort(dE, decreasing = T)[1:10]
#gracilis=634, australis=600, elegans=489, japonica=488, pacifica=450, antarctica=368, simplex=361, orientalis=349, elongata=344, minuta=333

# Get all the names with non unique epithet
tmpSpecies <- speciesList[which(epithet %in% names(dE[posDupE]))]
tmpSpecies <- as.data.frame(do.call(rbind, tmpSpecies), stringsAsFactors = FALSE)
dim(tmpSpecies) #168,300
round(nrow(tmpSpecies)/length(speciesList)*100, 2) #73.33%

# For each epithet
dupEps <- names(dE[posDupE])
minDLge <- numeric()
for(i in 1:length(dupEps)) #24,655
{
  # Get the minimum distance among its genera
  pos <- which(tmpSpecies$V2 == dupEps[i])
  gen <- tmpSpecies$V1[pos]
  
  dlGe <- stringdistmatrix(gen, gen, method = "dl")
  dlGe[upper.tri(dlGe, diag = TRUE)] <- NA
  dlGe <- na.omit(as.vector(dlGe))
  minDLge[i] <- min(dlGe)
}


## Summary
sum(minDLge<2)/length(minDLge)*100 #67   | 0.27
sum(minDLge<3)/length(minDLge)*100 #433  | 1.76
sum(minDLge<4)/length(minDLge)*100 #1491 | 6.05

sum(minDLep<2)/length(minDLep)*100 #303  | 1.45
sum(minDLep<3)/length(minDLep)*100 #1581 | 7.58
sum(minDLep<4)/length(minDLep)*100 #4031 | 19.34


## Histogram (Fig. 1B)
barplot(table(minDLge)/length(minDLge), space=0, col="grey80", ylim=c(0,.2), xlim=c(0,21), las=2, axisnames=F, cex.axis=1.25)
barplot(table(minDLep)/length(minDLep), space=0, col=adjustcolor("darkred",alpha.f=.1), axes=F, axisnames=F, add=T)
axis(side=1, at=seq(1,21,2)-.5, labels=NA, tck=-.035)
axis(side=1, at=seq(2,20,2)-.5, labels=NA, tck=-.02, lwd=0, lwd.ticks = 1)
axis(side=1, at=seq(1,21,2)-.5, labels=seq(1,21,2), tick=F, cex.axis=1.25)



# . Supplementary figures ----

## S1A How many genera each shared epithet has?
x <- cut(dE[posDupE], breaks=c(0,2,5,10,100,1000))
round(table(x)/length(dE[posDupE])*100, 2)
#10352	7985	3218	2984	116
#41.99% 32.39%	13.05%	12.10%	0.47%

a <- barplot(table(x)/length(dE[posDupE]), ylim=c(0,.5), las=1, axisnames=F)
axis(side=1, at=a, labels=c("[2]","[3-5]","[6-10]","[11-100]","[>100]"))
text(x=a, y=table(x)/length(dE[posDupE])+.05, table(x), cex=.8, col="darkred")


## S1B Which are the most common epithets?
library(ggplot2)
library(ggwordcloud)

p <- which(dE[posDupE]>100)
pp <- data.frame(nome=names(dE[posDupE][p]), freq=as.numeric(dE[posDupE][p]))

pp$angle <- sample(
  c(0, 90),
  size = nrow(pp),
  replace = TRUE,
  prob = c(0.8, 0.2)
)

ggplot(pp, aes(label=nome, size=freq, color=freq, angle = angle)) +
  geom_text_wordcloud(
    area_corr = TRUE,
    eccentricity = 1,
    rm_outside = F,
    shape="square"
  ) +
  scale_size_area(max_size=14) +
  scale_color_gradient(low = "grey40", high = "steelblue") +
  theme_minimal() +
  theme_void()


## S1C-D How many names repeat within and among phyla?
p <- which(epithet %in% names(dE[posDupE]))
tmpSpecies <- data.frame(species[p,], epithet=epithet[p])

df_counts <- tmpSpecies %>%
  group_by(epithet, phylum) %>%
  summarise(n = n(), .groups = "drop")

# Within phyla (C)
x <- table(df_counts$epithet)
x1 <- names(x[which(x==1)])

df_counts1 <- df_counts %>% filter(epithet %in% x1)
dim(df_counts1) #7561 names

a <- sort(table(df_counts1$phylum), decreasing = T)
length(a); max(a) #34 | 2531
barplot(rev(a[1:10]), horiz = T, las=1, xlim=c(0,3000))

# Among phyla (D)
xM <- names(x[which(x>1)])

df_countsM <- df_counts %>% filter(epithet %in% xM)
dim(df_countsM) #58930 names

a <- sort(table(df_countsM$phylum), decreasing = T)
length(a); max(a) #83 | 10478
barplot(rev(a[1:10]), horiz = T, las=1, xlim=c(0,12000))
 
#### END ####