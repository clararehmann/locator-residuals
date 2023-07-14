.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
options(warn = - 1)
library('ggplot2')
library('patchwork')
library('data.table')
library('dplyr')
library('rgeos')
library('mapproj')
library('RColorBrewer')
library('ggvoronoi')
library('rjson')
library('cowplot')
library('gridExtra')
library('grid')
library('sets')
library('sp')
library('tidyr')
library('ggnewscale')
library('scales')
library('stringr')

pal = c('#235862','#e740b4','#ffccba','#54d705','#fcfcfc')
palgen = colorRampPalette(pal)
colgen = colorRamp(pal)
figdir = '/home/crehmann/public_html/locator-residuals/' # place to save figures to

get_err <- function(predlocs) {
    params <- str_extract_all(predlocs,"\\(?[0-9,.]+\\)?")[[1]]
    sigma <- 1
    skew <- params[1]
    rep <- params[2]
    if (grepl('Locator', predlocs)) { skew <- params[2]
                                       rep <- params[3]
    #else if (grepl('run_', predlocs)) { skew <- params[5] }
    #else { skew <- params[4] }
    pd <- fread(predlocs)
    locs <- fread(paste0('data/biased_training/',
                         'sigma_1_bias_0_run_',rep,'_metadata.txt'))
    names(pd) <- c('xpred','ypred','sampleID')
    pd <- merge(locs, pd, on='sampleID')
    
    plocs=as.matrix(pd[,c("xpred","ypred")])
    tlocs=as.matrix(pd[,c("x","y")])
    dists=sapply(1:nrow(plocs),function(e) 
              spDistsN1(t(as.matrix(plocs[e,])),
              t(as.matrix(tlocs[e,])),longlat=F))
    xerr=pd$xpred - pd$x
    yerr=pd$ypred - pd$y
    }
                 
    else {
        pd <- fread(predlocs)
        colnames(pd) <- c('','','xpred','ypred','x','y')
        plocs=as.matrix(pd[,c("xpred","ypred")])
        tlocs=as.matrix(pd[,c("x","y")])
        dists=sapply(1:nrow(plocs),function(e) 
                  spDistsN1(t(as.matrix(plocs[e,])),
                  t(as.matrix(tlocs[e,])),longlat=F))
        xerr=pd$xpred - pd$x
        yerr=pd$ypred - pd$y
    }
                 
    data <- data.table(err=mean(dists),
                 xerr=mean(xerr),
                 yerr=mean(yerr),
                 rep=rep,
                 sigma=sigma,
                 skew=skew)
                 
    return(data)
}
                 
directory = "data/biased_training/SPASIBA-comp/Locator/sigma_1_*"
unf <- Sys.glob(paste0(directory, '*_predlocs.txt'))
df <- get_err(unf[1])
for (f in unf[2:length(unf)]){
    pd <- get_err(f)
    df <- rbind(df, pd)
}
                 
write.table(df, 'data/biased_training/SPASIBA-comp/Locator-predlocs.txt', sep='\t', row.names=FALSE)
locator <- df
locator$method <- 'locator'

directory = "data/biased_training/SPASIBA-comp/SPASIBA/*"
unf <- Sys.glob(paste0(directory, '*_predlocs.txt'))
df <- get_err(unf[1])
for (f in unf[2:length(unf)]){
    pd <- get_err(f)
    df <- rbind(df, pd)
}
                 
write.table(df, 'data/biased_training/SPASIBA-comp/SPASIBA-predlocs.txt', sep='\t', row.names=FALSE)
spasiba <- df
spasiba$method <- 'spasiba'
                     
df <- rbind(locator, spasiba)
df$skew <- as.numeric(df$skew)
df <- df[df$skew > 0.5,]

xerr <- ggplot(df, aes(x=factor(skew), y=xerr, fill=method)) + theme_classic() +
    geom_boxplot(lwd=0.1, outlier.size=0.001) + 
    scale_fill_manual(name="Method", values=palgen(4), labels=c('Locator','SPASIBA'))+
    ylab('Mean x-axis error')+
    xlab('Training set bias')

err <- ggplot(df, aes(x=factor(skew), y=err, fill=method)) + theme_classic() +
    geom_boxplot(lwd=0.1, outlier.size=0.001) + 
    scale_fill_manual(name="Method", values=palgen(4), labels=c('Locator','SPASIBA'))+
    ylab('Mean prediction error')+
    xlab('Training set bias') + scale_y_continuous(trans='log')

sp0 <- fread('data/biased_training//SPASIBA-comp/SPASIBA/skew_0.6_4_test50_train450_predlocs.txt')
colnames(sp0) <- c('','','xpred','ypred','x','y')
sp0 <- ggplot(sp0) + theme_bw() + 
    geom_point(mapping=aes(x=x, y=y), color='black', size=0.5) + 
    geom_segment(mapping=aes(x=x, y=y, xend=xpred, yend=ypred), color='black') +
    geom_point(mapping=aes(x=xpred, y=ypred), shape=21, color='black', fill=palgen(4)[2], stroke=0.5, size=2) +
    scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    coord_fixed()+xlab('')+ylab('')

sp1 <- fread('data/biased_training//SPASIBA-comp/SPASIBA/skew_0.9_4_test50_train450_predlocs.txt')
colnames(sp1) <- c('','','xpred','ypred','x','y')
sp1 <- ggplot(sp1) + theme_bw() + 
    geom_point(mapping=aes(x=x, y=y), color='black', size=0.5) + 
    geom_segment(mapping=aes(x=x, y=y, xend=xpred, yend=ypred), color='black') +
    geom_point(mapping=aes(x=xpred, y=ypred), shape=21, color='black', fill=palgen(4)[2], stroke=0.5, size=2) +
    scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    coord_fixed()+xlab('')+ylab('')

lc0 <- fread('data/biased_training/SPASIBA-comp/Locator/sigma_1_skew_0.6_run_4_predlocs.txt')
colnames(lc0) <- c('xpred','ypred','sampleID')
md <- fread('data/biased_training/sigma_1_bias_0_run_4_metadata.txt')
lc0 <- merge(lc0, md, by='sampleID')
lc0 <- ggplot(lc0) + theme_bw() + 
    geom_point(mapping=aes(x=x, y=y), color='black', size=0.5) + 
    geom_segment(mapping=aes(x=x, y=y, xend=xpred, yend=ypred), color='black') +
    geom_point(mapping=aes(x=xpred, y=ypred), shape=21, color='black', fill=palgen(4)[1], stroke=0.5, size=2) +
    scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    coord_fixed()+xlab('')+ylab('')

lc1 <- fread('data/biased_training/SPASIBA-comp/Locator/sigma_1_skew_0.9_run_6_predlocs.txt')
colnames(lc1) <- c('xpred','ypred','sampleID')
md <- fread('data/biased_training/sigma_1_bias_0_run_4_metadata.txt')
lc1 <- merge(lc1, md, by='sampleID')
lc1 <- ggplot(lc1) + theme_bw() + 
    geom_point(mapping=aes(x=x, y=y), color='black', size=0.5) + 
    geom_segment(mapping=aes(x=x, y=y, xend=xpred, yend=ypred), color='black') +
    geom_point(mapping=aes(x=xpred, y=ypred), shape=21, color='black', fill=palgen(4)[1], stroke=0.5, size=2) +
    scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    coord_fixed()+xlab('')+ylab('')

(xerr / err + plot_layout(guides='collect')) / (sp0 | lc0) / (sp1 | lc1) + 
    plot_annotation(tag_level='A') +
    plot_layout(heights=c(1, 0.5, 2, 2)) & 
    theme(text=element_text(size=8))
ggsave(file.path(figdir,'spasiba_comparison.pdf'), width = 6.5, height = 8, units='in')