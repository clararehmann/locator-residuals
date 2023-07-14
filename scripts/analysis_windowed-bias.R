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

distances = fromJSON(file = 'data/biased_training/distance.json')
distances = as.data.frame(distances)
density = fromJSON(file = 'data/biased_training/density.json')
density = as.data.frame(density)

tsa = fread('data/biased_training/sigma_1.0_train_450_meanresidual.txt')
tsa$tssize <- 450
tsb = fread('data/biased_training/sigma_1.0_train_900_meanresidual.txt')
tsb$tssize <- 900
tsc = fread('data/biased_training/sigma_1.0_train_2250_meanresidual.txt')
tsc$tssize <- 2250
tsd = fread('data/biased_training/sigma_1.0_train_4500_meanresidual.txt')
tsd$tssize <- 4500
tsres <- rbind(tsa, tsb, tsc, tsd)

tsres$dx2 <- tsres$dx ** 2
tsres$dy2 <- tsres$dy ** 2

s05 <-  ggplot() + theme_bw() + labs(title='0.5') +
        geom_point(mapping=(aes(x=density$'X0.5', y=distances$'X0.5')), size=0.001, color='gray')+
        geom_smooth(method='lm', mapping=(aes(x=density$'X0.5', y=distances$'X0.5')), color='black', lwd=0.25)+
        ylab('Prediction uncertainty')+xlab('Local sample density')+
        scale_y_continuous(limits=c(0, 60))
s06 <-  ggplot() + theme_bw() + labs(title='0.6') +
        geom_point(mapping=(aes(x=density$'X0.6', y=distances$'X0.6')), size=0.001, color='gray')+
        geom_smooth(method='lm', mapping=(aes(x=density$'X0.6', y=distances$'X0.6')), color='black', lwd=0.25)+
        ylab('')+xlab('Local sample density')+
        scale_y_continuous(limits=c(0, 60))
s07 <-  ggplot() + theme_bw() + labs(title='0.7') +
        geom_point(mapping=(aes(x=density$'X0.7', y=distances$'X0.7')), size=0.001, color='gray')+
        geom_smooth(method='lm', mapping=(aes(x=density$'X0.7', y=distances$'X0.7')), color='black', lwd=0.25)+
        ylab('')+xlab('Local sample density')+
        scale_y_continuous(limits=c(0, 60))
s08 <-  ggplot() + theme_bw() + labs(title='0.8') +
        geom_point(mapping=(aes(x=density$'X0.8', y=distances$'X0.8')), size=0.001, color='gray')+
        geom_smooth(method='lm', mapping=(aes(x=density$'X0.8', y=distances$'X0.8')), color='black', lwd=0.25)+
        ylab('')+xlab('Local sample density')+
        scale_y_continuous(limits=c(0, 60))
s09 <-  ggplot() + theme_bw() + labs(title='0.9') +
        geom_point(mapping=(aes(x=density$'X0.9', y=distances$'X0.9')), size=0.001, color='gray')+
        geom_smooth(method='lm', mapping=(aes(x=density$'X0.9', y=distances$'X0.9')), color='black', lwd=0.25)+
        ylab('')+xlab('Local sample density')+
        scale_y_continuous(limits=c(0, 60))

# error
xax <- ggplot()+theme_classic()+
        geom_boxplot(tsres,lwd=0.1, outlier.size=0.001, mapping=aes(x=factor(skew), y=dx2, fill=factor(tssize)))+
        scale_fill_manual(name='Training set size', values=spal)+ylab('Mean squared x-axis error')+xlab('Training set skew')
yax <- ggplot()+theme_classic()+
        geom_boxplot(tsres, lwd=0.1, outlier.size=0.001, mapping=aes(x=factor(skew), y=dy2, fill=factor(tssize)))+
        scale_fill_manual(name='Training set size', values=spal)+ylab('Mean squared y-axis error')+xlab('Training set skew')
mgn <- ggplot()+theme_classic()+
        geom_boxplot(tsres, lwd=0.1, outlier.size=0.001,  mapping=aes(x=factor(skew), y=mg, fill=factor(tssize)))+
        scale_fill_manual(name='Training set size', values=spal)+ylab('Mean prediction error')+xlab('Training set skew')

corr <- (s05 | s06 | s07 | s08 | s09) & 
        theme(text=element_text(size=8), axis.title.x=element_text(size=5), axis.title.y=element_text(size=5))
fg <- xax / (yax | mgn) + plot_layout(guides='collect', heights=c(7.5, 4)) & theme(text=element_text(size=8))

for (n in c(2, 3, 4, 5)){
    corr[[n]] <- corr[[n]]+plot_layout(tag_level = 'new')
}

fig <- corr/fg + plot_layout(heights = c(0.75, 5))
fig + plot_annotation(tag_level='A')
ggsave(file.path(figdir,'sigma_1_skewed_training_set_size_corr.pdf'), width = 6.5, height = 6, units='in')