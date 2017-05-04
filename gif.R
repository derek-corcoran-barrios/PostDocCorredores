library(raster)
library(rasterVis)

###Primero armar los interpolado lineales de cada bio

e <- readRDS("e.rds")
now <- getData('worldclim', var='bio', res=10)
now <- crop(now, e)
ac50 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=50)
ac50 <- crop(ac50, e)
ac70 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=70)
ac70 <- crop(ac70, e)


a <- list()

for(i in c(2,4,5,14,15)){
b <- setValues(ac50[[1]], NA)
bio.stack <- stack( mget( rep( "b" , 8 ) ) )
bio.stack[[1]] <- now[[i]]
bio.stack[[6]] <- ac50[[i]]
bio.stack[[8]] <- ac70[[i]]

a[[length(a)+1]] <- approxNA(bio.stack)
}
Bio2000 <- list()
#levelplot(a[[1]])
for(i in 1:5){
Bio2000[[i]] <- a[[i]][[1]]
}
Bio2000<- stack(unlist(Bio2000))
###Luego agregarlos en stacks por decada
bios<- names(Bio2000)
decade <- list()
bio <- list()
for(j in 1:8){
  for(i in 1:5){
  bio[[i]] <- a[[i]][[j]]
  decade[[j]]<- stack(unlist(bio))
  }
  names(decade[[j]]) <- bios 
}

years <- as.character(seq(2000, 2070, by = 10))
names(decade) <- years

####Luego extraer y modelar datos de especies

library(rgbif)
library(ggmap)
library(dismo)
library(dplyr)
Huemul <- occ_search(scientificName = "Hippocamelus bisulcus")
Huemul <- data.frame(lon =Huemul$data$decimalLongitude, lat = Huemul$data$decimalLatitude)
Huemul <- Huemul[complete.cases(Huemul),]
Huemul <- dplyr::filter(Huemul, lon != 0 & lat != 0)
saveRDS(Huemul, "Huemul.rds")

Huemul <- readRDS("Huemul.rds")

#Set Background points

set.seed(1963)
bg <- randomPoints(decade[[1]][[1]], 500 )
colnames(bg) <- c("lon", "lat")
##Extract 
now_values <- decade[[1]]
presvals <- raster::extract(x = decade[[1]], y = Huemul)
absvals <- raster::extract(decade[[1]], bg)

pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

##Fit model

library(biomod2)

myRespName <- 'Huemul'
myResp <- sdmdata$pb
myRespXY <- rbind(Huemul, bg)
myExpl <- decade[[1]]

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('MAXENT.Phillips'),
  models.options = myBiomodOption,
  NbRunEval=1,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))

HuemulProj <- list()

for(i in 1:8){
  myBiomodProj<- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = decade[[i]],
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd')
  HuemulProj[[i]] <- myBiomodProj@proj@val
}

HuemulProj <- stack(unlist(HuemulProj))

saveGIF(for(i in 1:8){plot(BinaryTransformation(HuemulProj[[i]], myBiomodModelOut@models.evaluation@val[2,2,1,1,1]), main = paste("Huemul", years[i]))}, movie.name = "HuemulBin.gif", img.name = "Rplot", convert = "convert", clean = TRUE)


HuemulProj <- HuemulProj/(max(cellStats(HuemulProj, "max")))
saveRDS(HuemulProj, "HuemulProj.rds")
brks <- seq(0, 1, length.out = 11)
nb <- length(brks)-1 
colors <- rev(heat.colors(nb))
years <- as.character(seq(2000, 2070, by = 10))

library(animation)
saveGIF(for(i in 1:8){plot(HuemulProj[[i]], col = colors, breaks = brks, main = paste("Huemul", years[i]))}, movie.name = "Huemul.gif", img.name = "Rplot", convert = "convert", clean = TRUE)

#######Puma


library(rgbif)
library(ggmap)
library(dismo)
library(dplyr)
Puma <- occ_search(scientificName = "Puma concolor")
Puma <- data.frame(lon =Puma$data$decimalLongitude, lat = Puma$data$decimalLatitude)
Puma <- Puma[complete.cases(Puma),]
Puma <- dplyr::filter(Puma, lon != 0 & lat != 0)
Puma <- dplyr::filter(Puma, lon >= -88 & lon <= -33.83333 & lat >= -60 & lat <= 10)

saveRDS(Puma, "Puma.rds")

Puma <- readRDS("Puma.rds")

#Set Background points

set.seed(1963)
bg <- randomPoints(decade[[1]][[1]], 500 )
colnames(bg) <- c("lon", "lat")
##Extract 
now_values <- decade[[1]]
presvals <- raster::extract(x = decade[[1]], y = Puma)
absvals <- raster::extract(decade[[1]], bg)

pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

##Fit model

library(biomod2)

myRespName <- 'Puma'
myResp <- sdmdata$pb
myRespXY <- rbind(Puma, bg)
myExpl <- decade[[1]]

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('MAXENT.Phillips'),
  models.options = myBiomodOption,
  NbRunEval=1,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))

PumaProj <- list()

for(i in 1:8){
  myBiomodProj<- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = decade[[i]],
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd')
  PumaProj[[i]] <- myBiomodProj@proj@val
}

PumaProj <- stack(unlist(PumaProj))

saveGIF(for(i in 1:8){plot(BinaryTransformation(PumaProj[[i]], myBiomodModelOut@models.evaluation@val[2,2,1,1,1]), main = paste("Puma", years[i]))}, movie.name = "PumaBin.gif", img.name = "Rplot", convert = "convert", clean = TRUE)


PumaProj <- PumaProj/(max(cellStats(PumaProj, "max")))
saveRDS(PumaProj, "PumaProj.rds")
brks <- seq(0, 1, length.out = 11)
nb <- length(brks)-1 
colors <- rev(heat.colors(nb))
years <- as.character(seq(2000, 2070, by = 10))

library(animation)
saveGIF(for(i in 1:8){plot(PumaProj[[i]], col = colors, breaks = brks, main = paste("Puma", years[i]))}, movie.name = "Puma.gif", img.name = "Rplot", convert = "convert", clean = TRUE)

###############Gorrion Zonotrichia capensis

library(rgbif)
library(ggmap)
library(dismo)
library(dplyr)
Chincol <- occ_search(scientificName = "Zonotrichia capensis")
Chincol <- data.frame(lon =Chincol$data$decimalLongitude, lat = Chincol$data$decimalLatitude)
Chincol <- Chincol[complete.cases(Chincol),]
Chincol <- dplyr::filter(Chincol, lon != 0 & lat != 0)
Chincol <- dplyr::filter(Chincol, lon >= -88 & lon <= -33.83333 & lat >= -60 & lat <= 10)

saveRDS(Chincol, "Chincol.rds")

Chincol <- readRDS("Chincol.rds")

#Set Background points

set.seed(1963)
bg <- randomPoints(decade[[1]][[1]], 500 )
colnames(bg) <- c("lon", "lat")
##Extract 
now_values <- decade[[1]]
presvals <- raster::extract(x = decade[[1]], y = Chincol)
absvals <- raster::extract(decade[[1]], bg)

pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

##Fit model

library(biomod2)

myRespName <- 'Chincol'
myResp <- sdmdata$pb
myRespXY <- rbind(Chincol, bg)
myExpl <- decade[[1]]

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('MAXENT.Phillips'),
  models.options = myBiomodOption,
  NbRunEval=1,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))

ChincolProj <- list()

for(i in 1:8){
  myBiomodProj<- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = decade[[i]],
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd')
  ChincolProj[[i]] <- myBiomodProj@proj@val
}

ChincolProj <- stack(unlist(ChincolProj))


saveGIF(for(i in 1:8){plot(BinaryTransformation(ChincolProj[[i]], myBiomodModelOut@models.evaluation@val[2,2,1,1,1]), main = paste("Chincol", years[i]))}, movie.name = "ChincolBin.gif", img.name = "Rplot", convert = "convert", clean = TRUE)

ChincolProj <- ChincolProj/(max(cellStats(ChincolProj, "max")))
saveRDS(ChincolProj, "ChincolProj.rds")
brks <- seq(0, 1, length.out = 11)
nb <- length(brks)-1 
colors <- rev(heat.colors(nb))
years <- as.character(seq(2000, 2070, by = 10))

library(animation)
saveGIF(for(i in 1:8){plot(ChincolProj[[i]], col = colors, breaks = brks, main = paste("Chincol", years[i]))}, movie.name = "Chincol.gif", img.name = "Rplot", convert = "convert", clean = TRUE)


