suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(scales)))
suppressMessages(suppressWarnings(require(raster)))
suppressMessages(suppressWarnings(require(sp)))
suppressMessages(suppressWarnings(require(MASS)))
suppressMessages(suppressWarnings(require(rgeos)))
suppressMessages(suppressWarnings(require(plyr)))
suppressMessages(suppressWarnings(require(progress)))
suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(tidyverse)))

get_err <- function(predlocs) {
    params <- str_extract_all(predlocs,"\\(?[0-9,.]+\\)?")[[1]]
    sigma <- 1
    rep <- params[2]
    if (grepl('uniform', predlocs)) { skew <- 0 }
    else if (grepl('run_', predlocs)) { skew <- params[5] }
    else { skew <- params[4] }
    pd <- fread(predlocs)
    locs <- fread(paste0('slimulation/out/metadata/',
                         'sigma_',sigma,'_bias_0_run_',rep,'_metadata.txt'))
    names(pd) <- c('xpred','ypred','sampleID')
    pd <- merge(locs, pd, on='sampleID')
    
    plocs=as.matrix(pd[,c("xpred","ypred")])
    tlocs=as.matrix(pd[,c("x","y")])
    dists=sapply(1:nrow(plocs),function(e) 
              spDistsN1(t(as.matrix(plocs[e,])),
              t(as.matrix(tlocs[e,])),longlat=F))
    xerr=pd$xpred - pd$x
    yerr=pd$ypred - pd$y
    
    data <- data.table(err=mean(dists),
                 xerr=mean(xerr),
                 yerr=mean(yerr),
                 rep=rep,
                 sigma=sigma,
                 skew=skew)
                 
    return(data)
}
                 

directory = "slimulation/downsampled/out/predlocs/"
unf <- Sys.glob(paste0(directory, '*_predlocs.txt'))
df <- get_err(unf[1])
for (f in unf[2:length(unf)]){
    pd <- get_err(f)
    df <- rbind(df, pd)
}
      
write.table(df, 'downsample_skew_predlocs.txt', sep='\t', row.names=FALSE)
