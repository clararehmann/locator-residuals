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

# read in SLiM centroids
predlocs <- data.frame()

for (run in c('0','1','2','3','4','5','6','8')){
    md <- fread(paste0('data/simulated_anopheles/',run,'_0.05scale_metadata.txt'))
    pd <- fread(paste0('data/simulated_anopheles/anopheles_slim_',run,'_centroids.txt'))
    df <- merge(md, pd, on='sampleID')
    predlocs <- rbind(predlocs, df)
    }

tlocs <- as.matrix(predlocs[,c('x','y')])
plocs <- as.matrix(predlocs[,c('kd_x','kd_y')])
predlocs$error <- sapply(1:nrow(plocs),function(e) spDistsN1(t(as.matrix(plocs[e,])),
                   t(as.matrix(tlocs[e,])),longlat = TRUE))                      
summary <- predlocs %>% group_by_at(c('x', 'y')) %>%
                  summarize(kd_x_mean = mean(kd_x),
                  kd_y_mean = mean(kd_y),
                  error = mean(error))
summary$type <- 'simulation'
write.table(predlocs, 'data/simulated_anopheles/slim_centroids.txt',quote=FALSE, col.names=NA, sep='\t')
                         
# read in Anopheles centroids
ag <- fread('data/anopheles_plasmodium/ag1000g_weighted_centroids.txt')
md <- fread('data/anopheles_plasmodium/ag1000g_v3_gambiae.txt')
ag <- merge(ag, md, on='sampleID')
ag <- ag[,c('x','y','kd_x','kd_y')]
ag[] <- lapply(ag, as.numeric)
ag <- ag[ag$y < 25,]
tlocs <- as.matrix(ag[,c('x','y')])
plocs <- as.matrix(ag[,c('kd_x','kd_y')])
ag$error <- sapply(1:nrow(plocs),function(e) spDistsN1(t(as.matrix(plocs[e,])),
                   t(as.matrix(tlocs[e,])),longlat = TRUE))

agsummary <- ag %>% group_by_at(c('x', 'y')) %>%
                  summarize(kd_x_mean = mean(kd_x),
                  kd_y_mean = mean(kd_y),
                  error = mean(error))
agsummary$type <- 'empirical'  

# combine dataframes and plot
summary <- rbind(summary, agsummary)
                   
mapplot <- ggplot()+theme_classic()+ labs(xlab='Longitude',ylab='Latitude',error='Mean prediction error') + 
            # plot map
            coord_map(projection = "mollweide",
                   xlim=c(min(na.omit(c(summary$kd_x_mean,summary$x)))-2.5,
                          max(na.omit(c(summary$kd_x_mean,summary$x)))+2.5),
                   ylim=c(min(na.omit(c(summary$kd_y_mean,summary$y)))-2.5,
                          max(na.omit(c(summary$kd_y_mean,summary$y)))+2.5))+
            geom_polygon(data=fortify(map),aes(x=long,y=lat,group=group),fill="#ebebeb",color="white",lwd=0.2)+
                            theme(axis.title = element_blank())+
            geom_point(summary, mapping=aes(x=x, y=y), color=aglight, size=0.1) +

            geom_point(summary, mapping=aes(x=kd_x_mean, y=kd_y_mean, size=error), alpha=0.5, color=agdark) +
            geom_segment(summary, mapping=aes(x=x, y=y, xend=kd_x_mean, yend=kd_y_mean), alpha=0.6, color=aglight) +
            guides(size = guide_legend(title = "Mean prediction\nerror (km)")) +
            facet_grid(~type) + theme( strip.background = element_blank())

ag$type <- 'empirical'
ag$alpha <- 0.9
predlocs$type <- 'simulated'
predlocs$alpha <- 1
df <- rbind(predlocs[,c('x','y','kd_x','kd_y','error','type', 'alpha')], ag)

errdist <- ggplot(df) + theme_classic() + geom_density(mapping=aes(x=error, group=type, fill=type, color=type)) +
    scale_fill_manual(values=c(aglight, agdark)) +
    scale_color_manual(values=c(aglight, agdark)) + 
    scale_x_continuous(trans='log10') +
    geom_density(ag, mapping=aes(x=error, group=type, fill=type, color=type), alpha=0.4) 

(mapplot+theme(plot.margin=unit(c(0,0,0,0), 'pt'))) / 
(plot_spacer() + errdist+ plot_spacer() + plot_layout(widths=c(0.25, 1, 0.25))) +
plot_layout(heights=c(2, 0.75), widths=c(1, 0.5)) + 
plot_annotation(tag_level='A') & theme(text=element_text(size=8))
ggsave(file.path(figdir,'simulated_anopheles_predlocs.pdf'), width = 6.5, height = 4, units='in')

# error magnitudes                   
print(c('mean Anopheles prediction error:', mean(ag$error)))
print(c('median Anopheles prediction error:',median(ag$error)))
print(c('mean SLiM prediction error:',mean(predlocs$error)))
print(c('median SLiM prediction error:',median(predlocs$error)))                   
                   
# correlation tests?
corr <- fread('data/simulated_anopheles/simulated_empirical_residual_correlation.txt')
perm <- fread('data/simulated_anopheles/simulated_empirical_residual_permutation.txt')
corr_perm = rbind(corr, perm)
corr_perm$Mag <- corr_perm$mag
                   
cat('Vector correlation test:\n')
print(paste('R = ', mean(corr$vector), '+/-', sd(corr$vector)))
print(wilcox.test(corr$vector, perm$vector))
                   
cat('Angle correlation:\n')
print(paste('R = ', mean(corr$angle), '+/-', sd(corr$angle)))
print(wilcox.test(corr$angle, perm$angle))

cat('Magnitude correlation:\n')
print(paste('R = ', mean(corr$mag), '+/-', sd(corr$mag)))
print(wilcox.test(corr$mag, perm$mag))
                   
