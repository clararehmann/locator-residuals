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

palgen <- colorRampPalette(pal[1:4])
figdir = '/home/crehmann/public_html/locator-residuals/'

# weight values as lambda, bandwidth change

df <- fread('data/biased_training/weighted/lambda_sample_weights_bw2.0.txt')
lplot <- ggplot(df[df$lambda %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)],
       aes(x=sample_weight, color=factor(lambda), fill=factor(lambda))) + theme_classic()+
  geom_density(alpha=0.25) + scale_x_continuous(trans='log10')+ 
    scale_color_manual(values = palgen(9), name='Lambda')+
    scale_fill_manual(values = palgen(9), name='Lambda')+
    xlab('log(sample weight)') + ylab('Density') +
    theme(legend.key.size = unit(0.5,"line"))
    #guides(fill=guide_legend(override.aes = list(size=0.1)))


df <- fread('data/biased_training/weighted/bandwidth_sample_weights_lambda5.0.txt')
bplot <- ggplot(df[df$bandwidth %in% c(1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4)], 
       aes(x=sample_weight, color=factor(bandwidth), fill=factor(bandwidth))) + theme_classic()+
  geom_density(alpha=0.25) + scale_x_continuous(trans='log10')+ 
    scale_color_manual(values = palgen(8), name='Bandwidth')+
    scale_fill_manual(values = palgen(8), name='Bandwidth')+
    xlab('log(sample weight)')+ ylab('Density') + 
     theme(legend.key.size = unit(0.5,"line"))
    #guides(fill=guide_legend(override.aes = list(size=0.1)))

(lplot | bplot) + plot_annotation(tag_level='A') & theme(text=element_text(size=8))
#ggsave(file.path(figdir,'lambda_bandwidth_sampleweightdensity.pdf'), width =6.5, height = 2.5, units='in')

ls <- fread('data/biased_training/weighted/skew_0.9_bandwidth_1.0_lambda_0.001_run_0_sample_weights.txt')
lm <- fread('data/biased_training/weighted/skew_0.9_bandwidth_1.0_lambda_2.511886431509581_run_0_sample_weights.txt')
lg <- fread('data/biased_training/weighted/skew_0.9_bandwidth_1.0_lambda_10.0_run_0_sample_weights.txt')

ls <- ggplot(ls, aes(x = x, y = y, color = sample_weight)) + theme_bw() + geom_point(size=0.5) +
            scale_color_gradientn(colors = palgen(100),
                                  name = "log(sample weight)",
                          trans='log10',limits=c(min(df$sample_weight),max(df$sample_weight)))+
            ylab('')+xlab('')+
            scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
            scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
            coord_fixed()+ 
                    theme(legend.key.width = unit(dev.size()[1] / 15, "line"),
                         legend.key.size = unit(0.75,"line"))
lm <- ggplot(lm, aes(x = x, y = y, color = sample_weight)) + theme_bw() + geom_point(size=0.5) +
            scale_color_gradientn(colors = palgen(100),
                                  name = "log(sample weight)",
                          trans='log10',limits=c(min(df$sample_weight),max(df$sample_weight)))+
            ylab('')+xlab('')+
            scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
            scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
            coord_fixed()+ 
                    theme(legend.key.width = unit(dev.size()[1] / 15, "line"),
                         legend.key.size = unit(0.75,"line"))
lg <- ggplot(lg, aes(x = x, y = y, color = sample_weight)) + theme_bw() + geom_point(size=0.5) +
            scale_color_gradientn(colors = palgen(100),
                                  name = "log(sample weight)",
                          trans='log10',limits=c(min(df$sample_weight),max(df$sample_weight)))+
            ylab('')+xlab('')+
            scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
            scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
            coord_fixed()+ 
                    theme(legend.key.width = unit(dev.size()[1] / 15, "line"),
                         legend.key.size = unit(0.75,"line"))


(lplot | bplot) / ((ls | lm | lg) + plot_layout(guides='collect')) + 
        plot_annotation(tag_level='A') + plot_layout(heights=c(1, 1))  & theme(text=element_text(size=8))
ggsave(file.path(figdir,'lambda_bandwidth_sampleweights.pdf'), width =6.5, height = 4, units='in')