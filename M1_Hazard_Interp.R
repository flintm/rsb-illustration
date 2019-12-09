# Seismic hazard curve input and interpolation for
# RSB Overview Illustration
# Madeleine Flint, 2019-07 to 2019-12
library(jsonlite) # rJSON package gives error
library(fields)
library(reshape2)

# load and set up log-log interpolation of USGS seismic hazard json file -------
haz.File   <- 'Charleston.BCboundary.2014DynamicConterm.json'
hazard     <- fromJSON(file.path("Data",haz.File))
Sa          <- hazard$response$metadata$xvalues[[1]]
MAFE.PGA   <- hazard$response$data[[1]]$yvalues[[1]] # Total
MAFE.SA0P2 <- hazard$response$data[[2]]$yvalues[[1]] # Total
MAFE.SA1P0 <- hazard$response$data[[3]]$yvalues[[1]] # Total
MAFE.SA2P0 <- hazard$response$data[[4]]$yvalues[[1]] # Total
T.orig <-  hazard$response$metadata$imt$value #[[4]]$yvalues[[1]] # Total
T.orig[T.orig=="PGA"] <- as.character(1e-12)
T.orig <- sub("P",".",T.orig)
T.orig <- sub("SA","",T.orig)
T.orig <- format(as.numeric(T.orig), scientific = FALSE, digits = 2, nsmall = 1)
n.T.orig <- length(T.orig)
MAFE <- matrix(c(MAFE.PGA, MAFE.SA0P2,MAFE.SA1P0,MAFE.SA2P0),nrow = 20, ncol = 4)

Salog <- log(Sa)
Tlog <- log(T.orig)
MAFElog <- log(MAFE)
MAFElog[MAFElog< -27] <- -27

haz.obj <- list(x=Salog,y=Tlog,z=MAFElog)

# Load dataframe to determine periods of applicable structures and interpolate -----
load(file.path("Data","M1.lateral.systems.RData"))
T <- unique(c(df.lat$T,1.16,0.66,0.46,0.5,0.94,1.32,0.69,0.36,T))
names(T) <- paste0('Sa',T)
T <- sort(T)
T.interp <- T
T.interp.log <- log(T.interp)
Sa.interp <- exp(seq(from=log(0.0025),to=log(3),length.out=100))
Sa.interp.log <- log(Sa.interp)
grid <- list(x=Sa.interp, y=T.interp)
grid.log <- list(x=Sa.interp.log, y = T.interp.log)
MAFE.interp.grid.log <- interp.surface.grid(haz.obj, grid.log)
MAFE.interp.grid <- exp(MAFE.interp.grid.log$z)
  
haz.interp <- as.data.frame(MAFE.interp.grid)
colnames(haz.interp) <- names(T)
haz.interp$Sa <- Sa.interp
haz.interp.melt <- melt(haz.interp, id.vars = "Sa", measure.vars = colnames(haz.interp)[grepl("Sa",colnames(haz.interp))], variable.name = "Period", value.name = "MAFE")
haz.interp.melt$Period <- as.character(haz.interp.melt$Period)
haz.interp.melt$Source <- "Interp"

haz.orig <- data.frame(Sa=Sa,MAFE=MAFE)
colnames(haz.orig)[2:(1+n.T.orig)] <- paste0('Sa',substr(as.character(T.orig),1,3))
colnames(haz.orig)[2] <- "PGA"
haz.orig.melt <- melt(haz.orig, id.vars = "Sa", variable.name = "Period", value.name = "MAFE")
haz.orig.melt$Period <- as.character(haz.orig.melt$Period)
haz.orig.melt$Source <- "Orig"

haz.all.melt <- rbind(haz.interp.melt,haz.orig.melt)
haz.all.melt$Period <- factor(haz.all.melt$Period)

haz.orig.melt$T <- as.numeric(gsub("[[:alpha:]]","",haz.orig.melt$Period))
haz.orig.melt[haz.orig.melt$Period=="PGA","T"] <- 0
haz.orig.melt$MAFE <- log(haz.orig.melt$MAFE)

save(haz.all.melt, Sa.interp, file = file.path("Data", sub(".json",".interp.RData",haz.File)))
write.table(haz.all.melt, file = file.path("Data", sub(".json",".interp.RData",haz.File)), 
            sep = "\t", row.names = FALSE)

# make hazard curve plot in log10 space ---------
Sa.interp.log10 <- log10(Sa.interp)
MAFE.interp.log10 <- log10(MAFE.interp.grid)
zmin <- -6.5
xmax <- 3
x.plot <- Sa.interp.log10[Sa.interp.log10<log10(xmax)]
y.plot <- T.interp
z.plot <- MAFE.interp.log10[Sa.interp.log10<log10(xmax),]
z.plot[z.plot < zmin] <- NA

p.melt <- haz.orig.melt[haz.orig.melt$Sa < xmax,c(1,5,3,4)]
colnames(p.melt)[1:3] <- c("x", "y", "z")
p.melt$z10 <- log10(exp(p.melt$z))

op <- par(bg = "white")
persp(x.plot,y.plot,z.plot, phi = 25, theta = 30, xlab = "Sa", ylab = "Period", zlab = "MAFE",
      xlim = c(-3,log10(xmax)), ylim = c(0,2), zlim = c(zmin,-1), ticktype = "detailed", d = 1.5)-> res
points(trans3d(x=log10(p.melt$x), y = p.melt$y, z = p.melt$z10, pmat = res))
par(op)
