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
                   
# plot map and predictions
pred_plot <- ggplot()+theme_classic()+ labs(xlab='Longitude',ylab='Latitude',error='Mean prediction error') + 
            coord_map(projection = "mollweide",
                   xlim=c(min(na.omit(c(summary$kd_x_mean,summary$x)))-2.5,
                          max(na.omit(c(summary$kd_x_mean,summary$x)))+2.5),
                   ylim=c(min(na.omit(c(summary$kd_y_mean,summary$y)))-2.5,
                          max(na.omit(c(summary$kd_y_mean,summary$y)))+2.5))+
            geom_polygon(data=fortify(map),aes(x=long,y=lat,group=group),fill="#ebebeb",color="white",lwd=0.2)+
                            theme(axis.title = element_blank())+
            geom_point(summary, mapping=aes(x=kd_x_mean, y=kd_y_mean, group=organism, color=organism, size=error), alpha=0.5) +
            guides(size = guide_legend(title = "Mean prediction error (km)")) +
            scale_color_manual(name="Mean predicted location", values=c("anopheles"=agdark, "plasmodium"=pflight), labels=c("Anopheles","Plasmodium")) +

            new_scale_colour() +
            geom_point(summary, mapping=aes(x=x, y=y, group=organism, color=organism), size=0.5) +
            scale_color_manual(name="True sampling location", values=c("anopheles"=aglight, "plasmodium"=pfdark), labels=c("Anopheles","Plasmodium")) +   
            geom_segment(summary, mapping=aes(x=x, y=y, xend=kd_x_mean, yend=kd_y_mean, group=organism, color=organism), alpha=0.6) +
            facet_wrap(~organism, labeller = labeller(organism = c('anopheles'='Anopheles','plasmodium'='Plasmodium'))) + theme( strip.background = element_blank() )
                   
# distribution of error magnitudes
agemperr <- pd[pd$organism=='anopheles',]$error
agd <- ggplot(mapping = aes(pd[pd$organism=='anopheles',]$error))+theme_classic()+
        geom_density(fill=aglight, color=agdark, alpha=0.5)+
        ylab('')+xlab('Anopheles prediction error (km)')+
        theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +scale_x_continuous(trans='log10')

pfemperr <- pd[pd$organism=='plasmodium',]$error
pfd <- ggplot(mapping = aes(pd[pd$organism=='plasmodium',]$error))+theme_classic()+
        geom_density(fill=pflight, color=pfdark, alpha=0.5)+
        ylab('')+xlab('Plasmodium prediction error (km)')+
        theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
        expand_limits(x=0, y=0)+scale_x_continuous(trans='log10')

print(c('mean Anopheles prediction error:', mean(pd[pd$organism=='anopheles',]$error)))
print(c('median Anopheles prediction error:',median(pd[pd$organism=='anopheles',]$error)))
print(c('mean Plasmodium prediction error:',mean(pd[pd$organism=='plasmodium',]$error)))
print(c('median Plasmodium prediction error:',median(pd[pd$organism=='plasmodium',]$error)))
                   
# correlation tests
corr = fread('data/anopheles_plasmodium/unweighted_anopheles_plasmodium_correlation.txt')
corr$type = rep('Spatial pairs', nrow(corr))
perm = fread('data/anopheles_plasmodium/unweighted_anopheles_plasmodium_permutation.txt')
perm$type = rep('Random pairs', nrow(perm))
corr_perm = rbind(corr, perm)
corr_perm$Mag <- corr_perm$mag
                   
# correlation plots
vector <- ggplot()+theme_classic()+
    geom_boxplot(corr_perm, lwd=0.1, outlier.size=0.1, mapping=aes(x=vector, fill=reorder(type, -vector)), show.legend = FALSE)+coord_flip()+
    scale_fill_manual(name="Test", values=c("Spatial pairs"=aglight, "Random pairs"=pfdark))+
    ylab('Vector correlation')+xlab('')+theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())

angle <- ggplot()+theme_classic()+
    geom_boxplot(corr_perm, lwd=0.1, outlier.size=0.1, mapping=aes(x=angle, fill=reorder(type, -vector)), show.legend = FALSE)+coord_flip()+
    scale_fill_manual(name="Test", values=c("Spatial pairs"=aglight, "Random pairs"=pfdark))+
    ylab('Angle correlation')+xlab('')+theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())

magn <- ggplot()+theme_classic()+
    geom_boxplot(corr_perm, lwd=0.1, outlier.size=0.25, mapping=aes(x=Mag, fill=reorder(type, -vector)))+coord_flip()+
    scale_fill_manual(name="Test", values=c("Spatial pairs"=aglight, "Random pairs"=pfdark))+
    ylab('Magnitude correlation')+xlab('')+theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())

corrs <- vector + angle + magn + plot_layout(ncol=3)
(pred_plot + theme(legend.spacing.y = unit(-0.01, "cm"),legend.spacing.x = unit(.01, 'cm')))/(agd + pfd) / corrs + plot_layout(heights=c(10, 2, 3)) + plot_annotation(tag_levels='A') & theme(text=element_text(size=8))
ggsave(file.path(figdir,'empirical_anopheles_plasmodium.pdf'), width = 6.5, height = 5, units='in')
                   
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
                   
# plot map and predictions
pred_plot <- ggplot()+theme_classic()+ labs(xlab='Longitude',ylab='Latitude',error='Mean prediction error') + 
            coord_map(projection = "mollweide",
                   xlim=c(min(na.omit(c(summary$kd_x_mean,summary$x)))-2.5,
                          max(na.omit(c(summary$kd_x_mean,summary$x)))+2.5),
                   ylim=c(min(na.omit(c(summary$kd_y_mean,summary$y)))-2.5,
                          max(na.omit(c(summary$kd_y_mean,summary$y)))+2.5))+
            geom_polygon(data=fortify(map),aes(x=long,y=lat,group=group),fill="#ebebeb",color="white",lwd=0.2)+
                            theme(axis.title = element_blank())+
            geom_point(summary, mapping=aes(x=kd_x_mean, y=kd_y_mean, group=organism, color=organism, size=error), alpha=0.5) +
            guides(size = guide_legend(title = "Mean prediction error (km)")) +
            scale_color_manual(name="Mean predicted location", values=c("anopheles"=agdark, "plasmodium"=pflight), labels=c("Anopheles","Plasmodium")) +

            new_scale_colour() +
            geom_point(summary, mapping=aes(x=x, y=y, group=organism, color=organism), size=0.5) +
            scale_color_manual(name="True sampling location", values=c("anopheles"=aglight, "plasmodium"=pfdark), labels=c("Anopheles","Plasmodium")) +   
            geom_segment(summary, mapping=aes(x=x, y=y, xend=kd_x_mean, yend=kd_y_mean, group=organism, color=organism), alpha=0.6) +
            facet_wrap(~organism, labeller = labeller(organism = c('anopheles'='Anopheles','plasmodium'='Plasmodium'))) + theme( strip.background = element_blank() )

agerror <- pd[pd$organism=='anopheles',]$error
agerror <- data.frame('Weighted'=agerror,
             'Unweighted'=agemperr) %>% pivot_longer(cols=c('Unweighted','Weighted'), names_to='Training', values_to='error')

pferror <- pd[pd$organism=='plasmodium',]$error
pferror <- data.frame('Weighted'=pferror,
             'Unweighted'=pfemperr) %>% pivot_longer(cols=c('Unweighted', 'Weighted'), names_to='Training', values_to='error')                                     
# distribution of error magnitudes
agd <- ggplot(agerror, mapping=aes(x=error, fill=Training, color=Training, linetype=Training, alpha=Training, size=Training)) + theme_classic() + 
    geom_density() + scale_fill_manual(values=c(agdark, aglight),guide = "none") + scale_alpha_manual(values=c(0, 0.5),guide = "none") +
    scale_linetype_manual(values=c('dashed','solid')) +
    scale_color_manual(values=c(agdark, aglight),guide = "none") + scale_size_manual(values=c(0.5, 1),guide = "none") +
        ylab('')+xlab('Anopheles prediction error (km)')+
        theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +scale_x_continuous(trans='log10')
pfd <- ggplot(pferror, mapping=aes(x=error, fill=Training, color=Training, linetype=Training, alpha=Training, size=Training)) + theme_classic() + 
    geom_density() + scale_fill_manual(values=c(pflight, pfdark),guide = "none") + scale_alpha_manual(values=c(0, 0.5),guide = "none") +
    scale_linetype_manual(values=c('dashed','solid')) +
    scale_color_manual(values=c(pflight, pfdark),guide = "none") + scale_size_manual(values=c(0.5, 1),guide = "none") +
        ylab('')+xlab('Anopheles prediction error (km)')+
        theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +scale_x_continuous(trans='log10')

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


# correlation plots
vector <- ggplot()+theme_classic()+
    geom_boxplot(corr_perm, lwd=0.1, outlier.size=0.1, mapping=aes(x=vector, fill=reorder(type, -vector)), show.legend = FALSE)+coord_flip()+
    scale_fill_manual(name="Test", values=c("Spatial pairs"=aglight, "Random pairs"=pfdark))+
    ylab('Vector correlation')+xlab('')+theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())

angle <- ggplot()+theme_classic()+
    geom_boxplot(corr_perm, lwd=0.1, outlier.size=0.1, mapping=aes(x=angle, fill=reorder(type, -vector)), show.legend = FALSE)+coord_flip()+
    scale_fill_manual(name="Test", values=c("Spatial pairs"=aglight, "Random pairs"=pfdark))+
    ylab('Angle correlation')+xlab('')+theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())

magn <- ggplot()+theme_classic()+
    geom_boxplot(corr_perm, lwd=0.1, outlier.size=0.25, mapping=aes(x=Mag, fill=reorder(type, -vector)))+coord_flip()+
    scale_fill_manual(name="Test", values=c("Spatial pairs"=aglight, "Random pairs"=pfdark))+
    ylab('Magnitude correlation')+xlab('')+theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())

corrs <- vector + angle + magn + plot_layout(ncol=3)
(pred_plot + theme(legend.spacing.y = unit(-0.01, "cm"),legend.spacing.x = unit(.01, 'cm')))/(agd + pfd + plot_layout(guides='collect')) / corrs + plot_layout(heights=c(10, 2, 3)) + plot_annotation(tag_levels='A') & theme(text=element_text(size=8))
ggsave(file.path(figdir,'weighted_empirical_anopheles_plasmodium.pdf'), width = 6.5, height = 5, units='in')                          