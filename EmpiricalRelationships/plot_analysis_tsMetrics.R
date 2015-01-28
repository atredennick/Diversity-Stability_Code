##########################################################################
#### Script to plot time series metrics and do quick analysis.
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-28-2015
####

#Clear the workspace
rm(list=ls(all=TRUE))

####
#### Load necessary libraries ----------------------------------------
####
library(ggplot2); library(plyr); library(reshape2); library(gridExtra); library(lme4)

####
#### Read in data (a list) and reformat as dataframe -----------------
####
quad_ts_metrics <- readRDS("quad_ts_metrics_AllSites.rds")
ts_metrics_df <- data.frame(quad=NA, stability=NA, comm_synchrony=NA, site=NA)
for(i in 1:length(quad_ts_metrics)){
  tmp_df <- as.data.frame(quad_ts_metrics[[i]])
  tmp_df$site <- names(quad_ts_metrics[i])
  ts_metrics_df <- rbind(ts_metrics_df, tmp_df)
}
ts_metrics_df <- ts_metrics_df[2:nrow(ts_metrics_df),]

####
#### Plot stability as a function of synchrony ---------------------
####
png("synchrony_vs_stability.png", width = 500, height = 1000)
g1 <- ggplot(ts_metrics_df, aes(x=comm_synchrony, y=log(stability)))+
  geom_point(size=4, aes(color=site))+
  geom_point(size=4, shape=1)+
  stat_smooth(method="lm", aes(color=site),size=1, se=FALSE)+
  stat_smooth(method="lm", color="black", size=2, linetype=1, se=FALSE)+
  xlab(expression(phi[C]))+
  ylab("ln(TS)")+
  scale_x_continuous(limits=c(0,1))+
  ggtitle("A")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1))
g2 <- ggplot(ts_metrics_df)+
  geom_boxplot(aes(x=site, y=comm_synchrony, fill=site), width=0.5)+
  xlab("Site")+
  ylab(expression(phi[C]))+
  ggtitle("B")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1))
g3 <- ggplot(ts_metrics_df)+
  geom_boxplot(aes(x=site, y=log(stability), fill=site), width=0.5)+
  xlab("Site")+
  ylab("ln(TS)")+
  ggtitle("C")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1))
three_plots <- grid.arrange(g1,g2,g3,ncol=1,nrow=3)
print(three_plots)
dev.off()

####
#### Linear model analysis with site random effect ------------------
####
mixed_mod <- lmer(log(stability)~comm_synchrony+(1|site), data=ts_metrics_df)
summary(mixed_mod)

#By site models
mod_az <- lm(log(stability)~comm_synchrony, 
             data=subset(ts_metrics_df, subset = site=="Arizona"))
summary(mod_az)

mod_id <- lm(log(stability)~comm_synchrony, 
             data=subset(ts_metrics_df, subset = site=="Idaho"))
summary(mod_id)

mod_ks <- lm(log(stability)~comm_synchrony, 
             data=subset(ts_metrics_df, subset = site=="Kansas"))
summary(mod_ks)

mod_mt <- lm(log(stability)~comm_synchrony, 
             data=subset(ts_metrics_df, subset = site=="Montana"))
summary(mod_mt)

mod_nm <- lm(log(stability)~comm_synchrony, 
             data=subset(ts_metrics_df, subset = site=="NewMexico"))
summary(mod_nm)




