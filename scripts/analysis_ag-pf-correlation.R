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
load("~/locator/data/cntrymap.Rdata")

pal = c('#235862','#e740b4','#ffccba','#54d705','#fcfcfc')
palgen = colorRampPalette(pal)
colgen = colorRamp(pal)
figdir = '/home/crehmann/public_html/locator-residuals/' # place to save figures to

pfdark=pal[3] 
pflight=pal[2]

agdark=pal[1]
aglight=pal[4]

cat('UNWEIGHTED:\n')
cat('--------------------------\n')

# calculate error in KM for all samples
ag <- fread('data/anopheles_plasmodium/ag1000g_centroids.txt')
md <- fread('data/anopheles_plasmodium/ag1000g_v3_gambiae.txt')
ag <- merge(ag, md, on='sampleID')
ag <- ag[,c('x','y','kd_x','kd_y')]
ag[] <- lapply(ag, as.numeric)
ag$organism <- 'anopheles'
ag <- ag[ag$y < 25,]

pf <- fread('data/anopheles_plasmodium/pf7_centroids.txt')
md <- fread('data/anopheles_plasmodium/pf7_africa_QC.txt')
pf <- merge(pf, md, on='sampleID')
pf <- pf[,c('x','y','kd_x','kd_y')]
pf[] <- lapply(pf, as.numeric)
pf$organism <- 'plasmodium'

pd <- rbind(ag, pf)
tlocs <- as.matrix(pd[,c('x','y')])
plocs <- as.matrix(pd[,c('kd_x','kd_y')])
pd$error <- sapply(1:nrow(plocs),function(e) spDistsN1(t(as.matrix(plocs[e,])),
                   t(as.matrix(tlocs[e,])),longlat = TRUE))

# summary dataframe for plotting
summary <- pd %>% group_by_at(c('x', 'y', 'organism')) %>%
                  summarize(kd_x_mean = mean(kd_x),
                  kd_y_mean = mean(kd_y),
                  error = mean(error))
                   
# correlation tests
corr = fread('data/anopheles_plasmodium/unweighted_anopheles_plasmodium_correlation.txt')
corr$type = rep('Spatial pairs', nrow(corr))
perm = fread('data/anopheles_plasmodium/unweighted_anopheles_plasmodium_permutation.txt')
perm$type = rep('Random pairs', nrow(perm))
corr_perm = rbind(corr, perm)
corr_perm$Mag <- corr_perm$mag
                   
cat('Vector correlation test:\n')
print(paste('R = ', mean(corr_perm[corr_perm$type == 'Spatial pairs',]$vector), '+/-', sd(corr_perm[corr_perm$type == 'Spatial pairs',]$vector)))
print(wilcox.test(corr_perm[corr_perm$type == 'Spatial pairs',]$vector, corr_perm[corr_perm$type == 'Random pairs',]$vector))
                   
cat('Angle correlation:\n')
print(paste('R = ', mean(corr_perm[corr_perm$type == 'Spatial pairs',]$angle), '+/-', sd(corr_perm[corr_perm$type == 'Spatial pairs',]$angle)))
print(wilcox.test(corr_perm[corr_perm$type == 'Spatial pairs',]$angle, corr_perm[corr_perm$type == 'Random pairs',]$angle))

cat('Magnitude correlation:\n')
print(paste('R = ', mean(corr_perm[corr_perm$type == 'Spatial pairs',]$Mag), '+/-', sd(corr_perm[corr_perm$type == 'Spatial pairs',]$Mag)))
print(wilcox.test(corr_perm[corr_perm$type == 'Spatial pairs',]$Mag, corr_perm[corr_perm$type == 'Random pairs',]$Mag))

# pairwise correlation between great circle distance and angle vector between Anopheles and Plasmodium predictions

anopheles <- pd[pd$organism == 'anopheles',]
anopheles$xvector <- anopheles$kd_x - anopheles$x
anopheles$yvector <- anopheles$kd_y - anopheles$y
anopheles$sampleID <- seq(nrow(anopheles))
plasmodium <- pd[pd$organism == 'plasmodium',]
plasmodium$xvector <- plasmodium$kd_x - plasmodium$x
plasmodium$yvector <- plasmodium$kd_y - plasmodium$y
plasmodium$sampleID <- seq(nrow(plasmodium))

combs <- expand.grid(anopheles$sampleID, plasmodium$sampleID)

ag <- anopheles[combs$Var1,]
pf <- plasmodium[combs$Var2,]

aglocs <- as.matrix(ag[,c('x','y')])
pflocs <- as.matrix(pf[,c('x','y')])
dists <- sapply(1:nrow(aglocs), function(e) spDistsN1(t(as.matrix(aglocs[e,])),
                                                      t(as.matrix(pflocs[e,])), longlat=TRUE))

hypot <- function(x, y) { return(sqrt(x**2 + y**2)) }

xvector <- ag[,'xvector'] + pf[,'xvector']
yvector <- ag[,'yvector'] + pf[,'yvector']
summags <- hypot(xvector, yvector)

df <- data.frame(dists = dists, mags = summags$xvector)

N <- sample(length(dists), 10000)
unweighted_vsumcorr <- ggplot(df[N,], aes(x=dists, y=mags)) + theme_classic()+
     geom_jitter(color='gray', alpha=1/(df[N,]$dists/500), width=0.1, height=0.1) + 
     geom_smooth(method='lm', color=agdark, fill=agdark) +scale_y_continuous(trans='log10') +scale_x_continuous(trans='log10') +
     ylab('Magnitude of vector sum (km)') + xlab('Distance between samples (km)')

T <- cor.test(dists, summags$xvector, method='spearman')
cat('Correlation between vector sums and spatial distance:\n')
print(T)
                   
                                      
cat('WEIGHTED:\n')
cat('--------------------------\n')

# calculate error in KM for all samples
ag <- fread('data/anopheles_plasmodium/ag1000g_weighted_centroids.txt')
md <- fread('data/anopheles_plasmodium/ag1000g_v3_gambiae.txt')
ag <- merge(ag, md, on='sampleID')
ag <- ag[,c('x','y','kd_x','kd_y')]
ag[] <- lapply(ag, as.numeric)
ag$organism <- 'anopheles'
ag <- ag[ag$y < 25,]

pf <- fread('data/anopheles_plasmodium/pf7_africa_weighted_centroids.txt')
md <- fread('data/anopheles_plasmodium/pf7_africa_QC.txt')
pf <- merge(pf, md, on='sampleID')
pf <- pf[,c('x','y','kd_x','kd_y')]
pf[] <- lapply(pf, as.numeric)
pf$organism <- 'plasmodium'

pd <- rbind(ag, pf)
tlocs <- as.matrix(pd[,c('x','y')])
plocs <- as.matrix(pd[,c('kd_x','kd_y')])
pd$error <- sapply(1:nrow(plocs),function(e) spDistsN1(t(as.matrix(plocs[e,])),
                   t(as.matrix(tlocs[e,])),longlat = TRUE))

# summary dataframe for plotting
summary <- pd %>% group_by_at(c('x', 'y', 'organism')) %>%
                  summarize(kd_x_mean = mean(kd_x),
                  kd_y_mean = mean(kd_y),
                  error = mean(error))
                   


print(c('mean Anopheles prediction error:', mean(pd[pd$organism=='anopheles',]$error)))
print(c('median Anopheles prediction error:',median(pd[pd$organism=='anopheles',]$error)))
print(c('mean Plasmodium prediction error:',mean(pd[pd$organism=='plasmodium',]$error)))
print(c('median Plasmodium prediction error:',median(pd[pd$organism=='plasmodium',]$error)))
                   
# correlation tests
corr = fread('data/anopheles_plasmodium/weighted_anopheles_plasmodium_correlation.txt')
corr$type = rep('Spatial pairs', nrow(corr))
perm = fread('data/anopheles_plasmodium/weighted_anopheles_plasmodium_permutation.txt')
perm$type = rep('Random pairs', nrow(perm))
corr_perm = rbind(corr, perm)
corr_perm$Mag <- corr_perm$mag
                   
cat('Vector correlation test:\n')
print(paste('R = ', mean(corr_perm[corr_perm$type == 'Spatial pairs',]$vector), '+/-', sd(corr_perm[corr_perm$type == 'Spatial pairs',]$vector)))
print(wilcox.test(corr_perm[corr_perm$type == 'Spatial pairs',]$vector, corr_perm[corr_perm$type == 'Random pairs',]$vector))
                   
cat('Angle correlation:\n')
print(paste('R = ', mean(corr_perm[corr_perm$type == 'Spatial pairs',]$angle), '+/-', sd(corr_perm[corr_perm$type == 'Spatial pairs',]$angle)))
print(wilcox.test(corr_perm[corr_perm$type == 'Spatial pairs',]$angle, corr_perm[corr_perm$type == 'Random pairs',]$angle))

cat('Magnitude correlation:\n')
print(paste('R = ', mean(corr_perm[corr_perm$type == 'Spatial pairs',]$Mag), '+/-', sd(corr_perm[corr_perm$type == 'Spatial pairs',]$Mag)))
print(wilcox.test(corr_perm[corr_perm$type == 'Spatial pairs',]$Mag, corr_perm[corr_perm$type == 'Random pairs',]$Mag))

# pairwise correlation between great circle distance and angle vector between Anopheles and Plasmodium predictions

anopheles <- pd[pd$organism == 'anopheles',]
anopheles$xvector <- anopheles$kd_x - anopheles$x
anopheles$yvector <- anopheles$kd_y - anopheles$y
anopheles$sampleID <- seq(nrow(anopheles))
plasmodium <- pd[pd$organism == 'plasmodium',]
plasmodium$xvector <- plasmodium$kd_x - plasmodium$x
plasmodium$yvector <- plasmodium$kd_y - plasmodium$y
plasmodium$sampleID <- seq(nrow(plasmodium))

combs <- expand.grid(anopheles$sampleID, plasmodium$sampleID)

ag <- anopheles[combs$Var1,]
pf <- plasmodium[combs$Var2,]

aglocs <- as.matrix(ag[,c('x','y')])
pflocs <- as.matrix(pf[,c('x','y')])
dists <- sapply(1:nrow(aglocs), function(e) spDistsN1(t(as.matrix(aglocs[e,])),
                                                      t(as.matrix(pflocs[e,])), longlat=TRUE))

hypot <- function(x, y) { return(sqrt(x**2 + y**2)) }

xvector <- ag[,'xvector'] + pf[,'xvector']
yvector <- ag[,'yvector'] + pf[,'yvector']
summags <- hypot(xvector, yvector)

df <- data.frame(dists = dists, mags = summags$xvector)

N <- sample(length(dists), 10000)
weighted_vsumcorr <- ggplot(df[N,], aes(x=dists, y=mags)) + theme_classic()+
     geom_jitter(color='gray', alpha=1/(df[N,]$dists/500), width=0.1, height=0.1) + 
     geom_smooth(method='lm', color=agdark, fill=agdark) +scale_y_continuous(trans='log10') +scale_x_continuous(trans='log10') +
     ylab('Magnitude of vector sum (km)') + xlab('Distance between samples (km)')

T <- cor.test(dists, summags$xvector, method='spearman')
cat('Correlation between vector sums and spatial distance:\n')
print(T)                   
                   
unweighted_vsumcorr + weighted_vsumcorr + plot_annotation(tag_levels='A') 
ggsave(file.path(figdir,'anopheles_plasmodium_distance_vectorsum_correlation.pdf'), width = 6.5, height = 5, units='in')
