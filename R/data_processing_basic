setwd("C:/Users/awray/Desktop/quicR")
library(reshape2)
library(tidyverse)
library(data.table)
library(wesanderson)
library(knitr)

df<-data.frame(read.csv(file="RTQ_simulated_data.csv"))
meta<-data.frame(read.csv(file="RTQ_simulated_metadata.csv"))
input information about your run (either how many cycles you ran or how many cycles you want to cutoff to):

n_cycles<-100
n_neg_ctrl_reps<-6
max_fl<-150000
meta_vars<-ncol(meta)
df2<-merge(meta[,1:5], df, by="well")

names(df2)[1:5]<-c("well", "sample_ID", "type", "concentration", "run_ID")
names(df2)[6:(n_cycles+5)]<-paste("t", c(1:n_cycles), sep="")

df_melt<-df2%>%
 pivot_longer(cols= -c(well, sample_ID, type, concentration, run_ID), names_to="timepoint",
              values_to="raw_fluor")
df_melt$percent_fluor<-df_melt$raw_fluor/max_fl
df_melt$timepoint<-str_remove(df_melt$timepoint, "t")
df_melt$timepoint<-as.integer(df_melt$timepoint)


raw_fluor_plot1<-ggplot(aes(x=timepoint, y=as.numeric(raw_fluor), color=as.factor(type), lty=as.factor(concentration)), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(well)))+ theme_classic()+
  labs(x="cycle",y="raw fluorescence", color="type",lty="concentration")
raw_fluor_plot


percent_fluor_plot<-ggplot(aes(x=as.integer(timepoint), y=percent_fluor, color=as.factor(type), lty=as.factor(concentration)), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(well)))+ theme_classic()+
  labs(x="cycle",y="percent fluorescence", color="type",lty="concentration")
percent_fluor_plot


neg_means<-df2%>% 
  filter(type == "neg_ctrl")%>%
  select(t2,t3,t4,t5,t6,t7,t8)%>%
  gather()%>%
  summarize(min=min(value),
            max=max(value),
            median=median(value),
            mean=mean(value),
            sd=sd(value))

print(neg_means)

threshold_neg5<-neg_means$mean+5*neg_means$sd
threshold_neg10<-neg_means$mean+10*neg_means$sd


all_means<-df2%>% 
  select(t2,t3,t4,t5,t6,t7,t8)%>%
  gather()%>%
  summarize(min=min(value),
            max=max(value),
            median=median(value),
            mean=mean(value),
            sd=sd(value))

threshold_all5<-all_means$mean+5*all_means$sd
threshold_all10<-all_means$mean+10*all_means$sd

print(c(threshold_neg5, threshold_neg10, threshold_all5, threshold_all10))
max_fl<-max(df[,2:n_cycles])


threshold_t1<-(mean(df2$t1)+(10*sd(df2$t1)))
threshold_t1


raw_fluor_plot<-ggplot(aes(x=as.integer(timepoint), y=raw_fluor, color=as.factor(type), lty=as.factor(concentration)), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(well)))+ theme_classic()+
  labs(x="cycle",y="raw fluorescence", color="type",lty="concentration")+
  geom_hline(yintercept=threshold_all10, lty=2, col="gray70")+
  annotate("text", label="threshold", x=4, y=threshold_all10+(threshold_all10/20))
raw_fluor_plot


percent_fluor_plot<-ggplot(aes(x=as.integer(timepoint), y=percent_fluor, color=as.factor(type), lty=as.factor(concentration)), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(well)))+ theme_classic()+
  labs(x="cycle",y="percent fluorescence", color="type",lty="concentration")+
  geom_hline(yintercept=threshold_all10/max_fl, lty=2, col="gray70")+
  annotate("text", label="threshold", x=4, y=threshold_all10/max_fl+(threshold_all10/(max_fl*20)))
percent_fluor_plot

df_samples_melt<-subset(df_melt, df_melt$type == "sample")

n_samples<-n_distinct(df_samples_melt$sample_ID)
pal2 <- wesanderson::wes_palette(n_samples, name="Zissou1",type = "continuous")
threshold_percent<-threshold_all10/max_fl

p<-ggplot(aes(x=timepoint, y=percent_fluor, group=well, color=sample_ID), data=df_samples_melt)+geom_point()+
  geom_line(alpha=0.8)+theme_classic()+scale_color_manual(values=pal2)+
  geom_hline(yintercept=threshold_percent, lty=2, col="gray70")+
  annotate("text", label="threshold", x=4, y=threshold_percent+0.02)+
  labs(color="Sample ID")+ylab("percent total fluorescence")+xlab("cycle")
p

df_over<-df_samples_melt%>%
  filter(raw_fluor > threshold_all10)

df_tt<-df_over%>%
  group_by(well)%>%
  slice(which.min(timepoint))

df_tt$timepoint_avg<-df_tt$timepoint-((df_tt$raw_fluor-threshold_all10)/threshold_all10)

df_tt$minutes<-(df_tt$timepoint_avg-1)*15
df_tt$hours<-df_tt$minutes/60
df_tt$seconds<-df_tt$minutes*60

df_tt$amyloid_rate<-(1/df_tt$hours)

avg_tt<-df_tt%>%
  group_by(sample_ID)%>%
  summarize(avg_tt_hours=mean(hours),
            sd_tt_hours=sd(hours),
            avg_amyloid_rate=mean(amyloid_rate),
            sd_amyloid_rate=sd(amyloid_rate))


avg_tt<-as.data.frame(avg_tt)

q<-ggplot(aes(x=as.factor(sample_ID), y=avg_tt_hours, color=as.factor(sample_ID)), data=avg_tt)+geom_point()+ geom_pointrange(mapping=aes(ymin=avg_tt_hours-sd_tt_hours, ymax=avg_tt_hours+sd_tt_hours), data=avg_tt, size=1)+ labs(color="Sample ID") + ylab("hours to threshold (mean \u00b1 SD) ")+xlab("sample")+theme_classic()
q
