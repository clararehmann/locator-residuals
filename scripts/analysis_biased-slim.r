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

pal = c('#235862','#e740b4','#ffccba','#54d705','#fcfcfc')
palgen = colorRampPalette(pal)
colgen = colorRamp(pal)
figdir = '/home/crehmann/public_html/locator-residuals/' # place to save figures to

cat('UNWEIGHTED:\n')
cat('--------------------------\n')

df <- fread('data/biased_training/all_skew_predlocs.txt')
df$skew_scale <- rescale(df$skew)
df$sigma_scale <- rescale(df$sigma)

# linear relationship between bias and error?
cat('Linear relationship between error and training set bias at sigma=1.0\n')
xfit <- glm(xerr ~ 0 + skew_scale, data=df[df$sigma==1,])
summary(xfit)

yfit <- glm(yerr ~ 0 + skew_scale, data=df[df$sigma==1,])
summary(yfit)

unweighted_sumdf <- df %>% group_by_at(c('sigma','skew','rep')) %>%
            summarize(xerr_mean = mean(xerr),
                      yerr_mean = mean(yerr),
                      err_mean = mean(err))

unweighted_error <- mean(unweighted_sumdf[unweighted_sumdf$sigma == 1 & unweighted_sumdf$skew == 0.9,]$err_mean)
unweighted_xerr <- mean(unweighted_sumdf[unweighted_sumdf$sigma == 1 & unweighted_sumdf$skew == 0.9,]$xerr_mean)
unweighted_yerr <- mean(unweighted_sumdf[unweighted_sumdf$sigma == 1 & unweighted_sumdf$skew == 0.9,]$yerr_mean)

print(paste('Mean x-axis error at bias=0.9, sigma=1.0:', unweighted_xerr))
print(paste('Mean y-axis error at bias=0.9, sigma=1.0:', unweighted_yerr))
print(paste('Mean error magnitude at bias=0.9, sigma=1.0:', unweighted_error))

# example training sets
ts5 = fread('data/biased_training/sigma_1_bias_0_run_0_uniform_ts.txt')
ts6 = fread('data/biased_training/sigma_1_bias_0_run_0_skew_0.6_ts.txt')
ts7 = fread('data/biased_training/sigma_1_bias_0_run_0_skew_0.7_ts.txt')
ts8 = fread('data/biased_training/sigma_1_bias_0_run_0_skew_0.8_ts.txt')
ts9 = fread('data/biased_training/sigma_1_bias_0_run_0_skew_0.9_ts.txt')

tsplot <- function(ts, titl) {
    plot = ggplot()+theme_bw()+
            labs(title=titl)+
            geom_point(ts, mapping=aes(x=x, y=y), size=0.001)+
            ylab('')+xlab('')+
            scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
            scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
            coord_fixed()
    return(plot)
    }

# skewed training sets
p5 <- tsplot(ts5, '0.5')
p6 <- tsplot(ts6, '0.6')
p7 <- tsplot(ts7, '0.7')
p8 <- tsplot(ts8, '0.8')
p9 <- tsplot(ts9, '0.9')

spal <- palgen(4)

# error boxplots
xax <- ggplot()+theme_classic()+
        geom_boxplot(unweighted_sumdf,lwd=0.1, outlier.size=0.001,  mapping=aes(x=factor(skew), y=xerr_mean, fill=factor(sigma)))+
        scale_fill_manual(name="Dispersal rate", values=rev(spal))+ylab('Mean x-axis error')+xlab('Training set bias')
yax <- ggplot()+theme_classic()+
        geom_boxplot(unweighted_sumdf, lwd=0.1, outlier.size=0.001,  mapping=aes(x=factor(skew), y=yerr_mean, fill=factor(sigma)))+
        scale_fill_manual(name="Dispersal rate", values=rev(spal))+ylab('Mean y-axis error')+xlab('Training set bias')
mgn <- ggplot()+theme_classic()+
        geom_boxplot(unweighted_sumdf, lwd=0.1, outlier.size=0.001,  mapping=aes(x=factor(skew), y=err_mean, fill=factor(sigma)))+
        scale_fill_manual(name="Dispersal rate", values=rev(spal))+ylab('Mean prediction error')+xlab('Training set bias')

ts <- p5 | p6 | p7 | p8 | p9
fig <- ts / xax / (yax | mgn) +
  plot_layout(guides = 'collect', heights=c(2, 7.5, 4)) & theme(text=element_text(size=8))
for (n in c(2, 3, 4, 5)){
    fig[[1]][[n]] <- fig[[1]][[n]]+plot_layout(tag_level = 'new')
    }
fig + plot_annotation(tag_level='A')
ggsave(file.path(figdir,'sigma_skewed_training_set.pdf'), width = 6.5, height = 6, units='in')

