########################################
##Accessing RaMP on MySQL through R
##Alissa Castleberry, Guy Brock
##Last updated: May 7, 2018
########################################


##set wd
setwd("~/Box Sync/OSU Brock Lab/Files for R")


########################################
##To access RaMP through Ewy's GUI
########################################

#install.packages("devtools")
library(devtools)
#install_github("mathelab/RAMP-DB")
library(RaMP)
RaMP::runRaMPapp(conpass="password")

###############################################################
##Interact directly with database file through R
##Uses RMySQL and dplyr packages
###############################################################

##From: https://www.r-bloggers.com/accessing-mysql-through-r/
library("RMySQL")
rampdb = dbConnect(MySQL(), user='root', password='password', dbname='ramp', host='localhost')
##list of all tables in RaMP
dbListTables(rampdb)
##list of all fields in a certain table
dbListFields(rampdb, 'pathway')
dbListFields(rampdb, 'analytesynonym')
dbListFields(rampdb, 'analytehaspathway')
dbListFields(rampdb, 'analytehasontology')

##disconnect after every connect
dbDisconnect(rampdb)

##From: https://datascienceplus.com/working-with-databases-in-r/
library("dplyr")
##view an individual table -- creates a dataframe
##tbl(databaseName, tableName)
pathway <- tbl(rampdb, "pathway")
##turn into data frame via tibble -- more efficient methods for data frames
pathway <- as_tibble(pathway)
#head(pathway)
#length(unique(pathway$type))


###############################################################
##Get relational database in matrix form
##From KEGG
##Compounds only - not genes
###############################################################

table(pathway$type)
kegg.names <- pathway$pathwayName[pathway$type == "kegg"]
kegg.rampID <- pathway$pathwayRampId[pathway$type == "kegg"]


analyteHasPathway <- tbl(rampdb, 'analytehaspathway')
analyteHasPathway <- as_tibble(analyteHasPathway)
#dim(analyteHasPathway)
#head(analyteHasPathway)
analyteHasPathway[analyteHasPathway$pathwayRampId == kegg.rampID[1],]

##only want to pull out compounds
compound.idx <- grep("C", analyteHasPathway$rampId)
compoundHasPathway <- analyteHasPathway[compound.idx, ]

compoundHasPathway[compoundHasPathway$pathwayRampId == kegg.rampID[1],]
compoundHasKegg <- compoundHasPathway[compoundHasPathway$pathwayRampId %in% kegg.rampID, ]
#dim(compoundHasKegg)
nCompound <- length(unique(compoundHasKegg$rampId))
nKegg <- length(unique(compoundHasKegg$pathwayRampId))
#nKegg
#nCompound


## How do we get the metabolites in each pathway and include as a matrix (0/1)
## call this 'keggCompound'
keggCompound <- matrix(0, nrow = nKegg, ncol = nCompound)
rownames(keggCompound) <- unique(compoundHasKegg$pathwayRampId)
colnames(keggCompound) <- unique(compoundHasKegg$rampId)

for (i in 1:nrow(keggCompound)) {
  pname <- rownames(keggCompound)[i]
  cnames <- compoundHasKegg$rampId[compoundHasKegg$pathwayRampId == pname]
  keggCompound[i,cnames] <- 1
}

#apply(keggCompound, 1, sum)
## first one is 'metabolic pathways' and not super helpful so probably remove it 
keggCompound <- keggCompound[-which(rownames(keggCompound) == 'RAMP_P_000051210'),]
hist(apply(keggCompound, 1, sum))

##saves to working directory
save(keggCompound, file = 'keggCompound.rda')
## to load use: load('keggCompound.rda')

###############################################################
##Get relational database in matrix form
##From HMBD
##Compounds only - not genes
###############################################################

table(pathway$type)
hmdb.names <- pathway$pathwayName[pathway$type == "hmdb"]
hmdb.rampID <- pathway$pathwayRampId[pathway$type == "hmdb"]

analyteHasPathway <- tbl(rampdb, 'analytehaspathway')
analyteHasPathway <- as_tibble(analyteHasPathway)
#dim(analyteHasPathway)
#head(analyteHasPathway)
analyteHasPathway[analyteHasPathway$pathwayRampId == hmdb.rampID[1],]

##only want to pull out compounds
compound.idx <- grep("C", analyteHasPathway$rampId)
compoundHasPathway <- analyteHasPathway[compound.idx, ]

compoundHasPathway[compoundHasPathway$pathwayRampId == hmdb.rampID[1],]
compoundHasHMDB <- compoundHasPathway[compoundHasPathway$pathwayRampId %in% hmdb.rampID, ]
#dim(compoundHasHMDB)
nCompound <- length(unique(compoundHasHMDB$rampId))
nHMDB <- length(unique(compoundHasHMDB$pathwayRampId))
#nHMDB
#nCompound


## How do we get the metabolites in each pathway and include as a matrix (0/1)
## call this 'hmdbCompound'
hmdbCompound <- matrix(0, nrow = nHMDB, ncol = nCompound)
rownames(hmdbCompound) <- unique(compoundHasHMDB$pathwayRampId)
colnames(hmdbCompound) <- unique(compoundHasHMDB$rampId)

for (i in 1:nrow(hmdbCompound)) {
  pname <- rownames(hmdbCompound)[i]
  cnames <- compoundHasHMDB$rampId[compoundHasHMDB$pathwayRampId == pname]
  hmdbCompound[i,cnames] <- 1
}

hApply <- apply(hmdbCompound, 1, sum)
## might have to remove some? 
hmdbCompound <- hmdbCompound[-which(rownames(hmdbCompound) == 'RAMP_P_000051210'),]
hist(hApply)

#filter(pathway, pathwayRampId == 'RAMP_P_000050796')

##saves to working directory
save(hmdbCompound, file = 'hmdbCompound.rda')
## to load use: load('hmdbCompound.rda')


###############################################################
##Get relational database in matrix form
##From Reactome
##Compounds only - not genes
###############################################################

table(pathway$type)
reactome.names <- pathway$pathwayName[pathway$type == "reactome"]
reactome.rampID <- pathway$pathwayRampId[pathway$type == "reactome"]


analyteHasPathway <- tbl(rampdb, 'analytehaspathway')
analyteHasPathway <- as_tibble(analyteHasPathway)
#dim(analyteHasPathway)
#head(analyteHasPathway)
analyteHasPathway[analyteHasPathway$pathwayRampId == reactome.rampID[1],]

##only want to pull out compounds
compound.idx <- grep("C", analyteHasPathway$rampId)
compoundHasPathway <- analyteHasPathway[compound.idx, ]

compoundHasPathway[compoundHasPathway$pathwayRampId == reactome.rampID[1],]
compoundHasReactome <- compoundHasPathway[compoundHasPathway$pathwayRampId %in% reactome.rampID, ]
#dim(compoundHasKegg)
nCompound <- length(unique(compoundHasReactome$rampId))
nReactome <- length(unique(compoundHasReactome$pathwayRampId))
#nReactome
#nCompound


## How do we get the metabolites in each pathway and include as a matrix (0/1)
## call this 'keggCompound'
reactomeCompound <- matrix(0, nrow = nReactome, ncol = nCompound)
rownames(reactomeCompound) <- unique(compoundHasReactome$pathwayRampId)
colnames(reactomeCompound) <- unique(compoundHasReactome$rampId)

for (i in 1:nrow(reactomeCompound)) {
  pname <- rownames(reactomeCompound)[i]
  cnames <- compoundHasReactome$rampId[compoundHasReactome$pathwayRampId == pname]
  reactomeCompound[i,cnames] <- 1
}

rApply <- apply(reactomeCompound, 1, sum)
##third entry is "metabolism" so remove
reactomeCompound <- reactomeCompound[-which(rownames(reactomeCompound) == 'RAMP_P_000050796'),]
hist(rApply)

##saves to working directory
save(reactomeCompound, file = 'reactomeCompound.rda')
## to load use: load('keggCompound.rda')

###############################################################
##Get relational database in matrix form
##From Wiki
##Compounds only - not genes
###############################################################

table(pathway$type)
wiki.names <- pathway$pathwayName[pathway$type == "wiki"]
wiki.rampID <- pathway$pathwayRampId[pathway$type == "wiki"]


analyteHasPathway <- tbl(rampdb, 'analytehaspathway')
analyteHasPathway <- as_tibble(analyteHasPathway)
#dim(analyteHasPathway)
#head(analyteHasPathway)
analyteHasPathway[analyteHasPathway$pathwayRampId == wiki.rampID[1],]

##only want to pull out compounds
compound.idx <- grep("C", analyteHasPathway$rampId)
compoundHasPathway <- analyteHasPathway[compound.idx, ]

compoundHasPathway[compoundHasPathway$pathwayRampId == wiki.rampID[1],]
compoundHasWiki <- compoundHasPathway[compoundHasPathway$pathwayRampId %in% wiki.rampID, ]
#dim(compoundHasKegg)
nCompound <- length(unique(compoundHasWiki$rampId))
nWiki <- length(unique(compoundHasWiki$pathwayRampId))
#nWiki
#nCompound


## How do we get the metabolites in each pathway and include as a matrix (0/1)
## call this 'keggCompound'
wikiCompound <- matrix(0, nrow = nWiki, ncol = nCompound)
rownames(wikiCompound) <- unique(compoundHasWiki$pathwayRampId)
colnames(wikiCompound) <- unique(compoundHasWiki$rampId)

for (i in 1:nrow(wikiCompound)) {
  pname <- rownames(wikiCompound)[i]
  cnames <- compoundHasWiki$rampId[compoundHasWiki$pathwayRampId == pname]
  wikiCompound[i,cnames] <- 1
}

wApply <- apply(wikiCompound, 1, sum)
filter(pathway, pathwayRampId == 'RAMP_P_000048981')
##second one is "Biochemical Pathways Part I so I removed it
filter(pathway, pathwayRampId == 'RAMP_P_000048918')
##amino acid metabolism is another outlier... remove? 
##did not remove
wikiCompound <- wikiCompound[-which(rownames(wikiCompound) == 'RAMP_P_000048981'),]
hist(wApply)


##saves to working directory
save(wikiCompound, file = 'wikiCompound.rda')
## to load use: load('keggCompound.rda')

