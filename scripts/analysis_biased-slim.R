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

cat('BIASED TRAINING SETS:\n')
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
    # function to make a training set scatterplot
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

cat('SAMPLE WEIGHTING:\n')
cat('--------------------------\n')

df <- fread('data/biased_training/weighted_all_skew_predlocs.txt')
# rescale to same 0-1 linear as unweighted (0.5 not included)
df$skew_scaled <- rescale(df$skew, to=c(1/5, 1))
df$lm <- round(df$lm, 6)

lmCoefficient <- function(skew_scaled, xerr) {
    # function to return lm coefficient between bias & error
    if (length(skew_scaled) > 0) {
        sumres <- summary(lm(xerr ~ 0 + skew_scaled))
        return(sumres$coefficient[1])
    }
    else { return(NA) }
}
    
lmSignificance <- function(skew_scaled, xerr) {
    # function to return lm significance between bias & error
    if (length(skew_scaled) > 0) {
        sumres <- summary(lm(xerr ~ 0 + skew_scaled))
        if (is.numeric(sumres$fstatistic[1L])) {
            return(pf(sumres$fstatistic[1L], sumres$fstatistic[2L], sumres$fstatistic[3L], lower.tail = FALSE))
        }
         else { return(NA) }
    }
    else { return(NA) }
}

# summarize error and lm

lmdf <- df %>% group_by_at(c('bw','lm')) %>%
         summarize(R = lmCoefficient(skew_scaled, xerr),
                   p = lmSignificance(skew_scaled, xerr))

ylmdf <- df %>% group_by_at(c('bw','lm')) %>%
         summarize(R = lmCoefficient(skew_scaled, yerr),
                   p = lmSignificance(skew_scaled, yerr))

err_df <- df %>% group_by_at(c('bw','lm','skew')) %>%
         summarize(skew = mean(skew),
                   xerr_mean = mean(xerr),
                   xerr_median = median(xerr),
                   err_mean = mean(err),
                   err_median = median(err),
                   yerr_mean = mean(yerr))

sum_df <- df %>% group_by_at(c('bw','lm','skew','rep')) %>%
         summarize(skew = mean(skew),
                   xerr_mean = mean(xerr),
                   xerr_median = median(xerr),
                   err_mean = mean(err),
                   err_median = median(err),
                   yerr_mean = mean(yerr))


# tile plots of grid search error values 
err <- ggplot(err_df[err_df$skew == 0.9,], aes(x=bw, y=factor(lm), fill=err_mean, width=0.5, height=0.1)) + theme_classic() + 
            geom_tile(height=1) + scale_fill_viridis_c(name="Mean\nerror") + ylab('')+xlab('Bandwidth')+
            theme(legend.key.width = unit(dev.size()[1] / 15, "line"),
                         legend.key.size = unit(0.75,"line"))

xerr <- ggplot(err_df[err_df$skew == 0.9,], aes(x=bw, y=factor(lm), fill=xerr_mean, width=0.5, height=0.1)) + theme_classic() + 
            geom_tile(height=1) + scale_fill_viridis_c(name="Mean\nx-axis\nerror")+ylab('Lambda')+xlab('Bandwidth')+
            theme(legend.key.width = unit(dev.size()[1] / 15, "line"),
                         legend.key.size = unit(0.75,"line"))

yerr <- ggplot(err_df[err_df$skew == 0.9,], aes(x=bw, y=factor(lm), fill=yerr_mean, width=0.5, height=0.1)) + theme_classic() + 
            geom_tile(height=1) + scale_fill_viridis_c(name="Mean\ny-axis\nerror")+ylab('')+xlab('Bandwidth')+
            theme(legend.key.width = unit(dev.size()[1] / 15, "line"),
                         legend.key.size = unit(0.75,"line"))

(xerr + yerr + err) + plot_annotation(tag_level='A') & theme(text=element_text(size=8))
ggsave(file.path(figdir,'weight_grid_search.pdf'), width = 6.5, height = 2, units='in')

# which weights work best?
df <- merge(err_df, lmdf, on=c('bw','lm'))
xdf <- df[df$skew == 0.9,]
xdf <- xdf[xdf$err_mean <= unweighted_error,]
print(head(xdf[order(xdf$err_mean), ]))

bw = 0.25
lam = 1e-5

print(lmdf[(lmdf$bw == bw) & (lmdf$lm == lam),])

cat('DOWNSAMPLING:\n')
cat('--------------------------\n')

downsampled <- fread('data/biased_training/downsample_skew_predlocs.txt')
downsampled$skew_scale <- rescale(downsampled$skew, to=c(1/5, 1))
downsampled$sigma_scale <- rescale(downsampled$sigma)
downsampled$weight <- 1
# rename columns so dataframes play nice
downsampled$xerr_mean <- downsampled$xerr
downsampled$yerr_mean <- downsampled$yerr
downsampled$err_mean <- downsampled$err
xfit <- glm(xerr ~ 0 + skew_scale, data=downsampled)
summary(xfit)

# make comparison plot!
rescaled_unwsumdf <- unweighted_sumdf[unweighted_sumdf$skew > 0.5 & unweighted_sumdf$sigma == 1,]
rescaled_unwsumdf <- rescaled_unwsumdf[,c('skew','rep','xerr_mean','err_mean','yerr_mean')]
rescaled_unwsumdf$weight <- 0    
weighted_df <- sum_df[sum_df$bw == bw & sum_df$lm == lam,]
weighted_df <- weighted_df[,c('skew','rep','xerr_mean','err_mean','yerr_mean')]
weighted_df$weight <- 2
comp_df <- rbind(rescaled_unwsumdf, weighted_df, downsampled)

xplot <- ggplot()+theme_classic()+
            geom_boxplot(comp_df, lwd=0.1, mapping=aes(x=factor(skew), y=xerr_mean, fill=factor(weight))) + 
            scale_fill_manual(name='', values=palgen(4), labels=c('Unweighted','Downsampled','Weighted'))+ylab('Mean x-axis error')+xlab('Training set bias')
yplot <- ggplot()+theme_classic()+
            geom_boxplot(comp_df, lwd=0.1, mapping=aes(x=factor(skew), y=yerr_mean, fill=factor(weight))) + 
            scale_fill_manual(name='', values=palgen(4), labels=c('Unweighted','Downsampled','Weighted'))+ylab('Mean y-axis error')+xlab('Training set bias')
eplot <- ggplot()+theme_classic()+
            geom_boxplot(comp_df, lwd=0.1, mapping=aes(x=factor(skew), y=err_mean, fill=factor(weight))) + 
            scale_fill_manual(name='', values=palgen(4), labels=c('Unweighted','Downsampled','Weighted'))+ylab('Mean error')+xlab('Training set bias')

xR <- lmdf[lmdf$bw == bw & lmdf$lm == lam,]$R
xp<- lmdf[lmdf$bw == bw & lmdf$lm == lam,]$p

yR <- ylmdf[ylmdf$bw == bw & ylmdf$lm == lam,]$R
yp<- ylmdf[ylmdf$bw == bw & ylmdf$lm == lam,]$p

df <- fread('data/biased_training/weighted/skew_0.9_bandwidth_0.25_lambda_0.0001_run_0_predlocs.txt')
df$xpred <- df$x
df$ypred <- df$y
df <- df[,c('xpred','ypred','sampleID')]
pd <- fread('data/biased_training/sigma_1_bias_0_run_0_metadata.txt')
weighted <- merge(df, pd)
weighted$mode <- 2

wpred <- ggplot(weighted) + theme_bw() + 
    geom_point(mapping=aes(x=x, y=y), color='black', size=0.5) + 
    geom_segment(mapping=aes(x=x, y=y, xend=xpred, yend=ypred), color='black') +
    geom_point(mapping=aes(x=xpred, y=ypred), shape=21, color='black', fill=palgen(4)[3], stroke=0.5, size=2) +
    scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    coord_fixed()+xlab('')+ylab('')

df <- fread('data/biased_training/sigma_1_bias_0_run_0_skew_0.9_predlocs.txt')
df$xpred <- df$x
df$ypred <- df$y
df <- df[,c('xpred','ypred','sampleID')]
pd <- fread('data/biased_training/sigma_1_bias_0_run_0_metadata.txt')
df <- merge(df, pd)

unwpred <- ggplot(df) + theme_bw() + 
    geom_point(mapping=aes(x=x, y=y), color='black', size=0.5) + 
    geom_segment(mapping=aes(x=x, y=y, xend=xpred, yend=ypred), color='black') +
    geom_point(mapping=aes(x=xpred, y=ypred), shape=21, color='black', fill=palgen(4)[1], stroke=0.5, size=2) +
    scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    coord_fixed()+xlab('')+ylab('')

df <- fread('data/biased_training/downsampled/skew_0.9_run_0_predlocs.txt')
df$xpred <- df$x
df$ypred <- df$y
df <- df[,c('xpred','ypred','sampleID')]
pd <- fread('data/biased_training/sigma_1_bias_0_run_0_metadata.txt')
df <- merge(df, pd)

dspred <- ggplot(df) + theme_bw() + 
    geom_point(mapping=aes(x=x, y=y), color='black', size=0.5) + 
    geom_segment(mapping=aes(x=x, y=y, xend=xpred, yend=ypred), color='black') +
    geom_point(mapping=aes(x=xpred, y=ypred), shape=21, color='black', fill=palgen(4)[2], stroke=0.5, size=2) +
    scale_x_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    scale_y_continuous(limits=c(0, 50), breaks=c(0, 25, 50))+
    coord_fixed()+xlab('')+ylab('')

(xplot / (yplot | eplot) + plot_layout(guides='collect')) / (unwpred | dspred | wpred) + plot_annotation(tag_level='A') + plot_layout(heights=c(7.5, 4, 5)) & theme(text=element_text(size=8))
ggsave(file.path(figdir,'weight_error.pdf'), width = 6.5, height = 8, units='in')