####FROM THE MACHINE######
#after the run, export the raw data using the "excel report" button.
#do not include the run info. do not transpose.
#to use this script, make sure
#that data are in the format where the COLUMN A is wells, COLUMN B is "content"
#(eg., whatever you named the samples in the machine's plate layout), with ROW 1
#showing labels (e.g., 1C should say something like "raw data" then the filter info).
#ROW 2 should show the time period -- it doesn't really matter how it is formatted. 
#you will also need a metadata file that corresponds to your Sample ID's somehow.

##working on a simulated dataset so this code can run as a vignette###
#note that some info is missing in this script 

setwd("~/Desktop/RT_QuIC/data")

library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyverse)
library(xlsx)
library(DataCombine)

######clean up that data######


#check out data structure:
dim(df1)

#now, for this run, we actually did more than 24 hours but I want to cut it off
#at the 24 hr mark so I will select only the cycles that we want:

df1<-df1[,1:99]
df2<-df2[,1:99]
df3<-df3[,1:99]

#add run date info for merging later:

#combine, there's probably a more elegant way to do this but this works:
df12<-rbind(df1, df2)
df123<-rbind(df12, df3)

#rename columns
names(df123)[1]<-"Well"
names(df123)[2]<-"SampleID"

#rename dataset:
df<-df123

#I personally do not like numbers in the labels for a column, so we are 
#pasting a "t" for timepoint, even though it will be removed later
names(df)[3:99]<-paste("t", c(1:97), sep="")

#get name formats the same for merging with metadata file
#you have to run this twice because of the negative control name.
df$SampleID <- sub(" ", "_", df$SampleID)
df$SampleID <- sub(" ", "_", df$SampleID)

#I also find single digit sample IDs to be extremely distressing. Use
#find & replace to fix it so they can be in order:

replaces<-data.frame(from=c("Sample_X1", "Sample_X2", "Sample_X3","Sample_X4",
                            "Sample_X5","Sample_X6","Sample_X7","Sample_X8",
                            "Sample_X9"), 
                     to=c("Sample_X01","Sample_X02","Sample_X03","Sample_X04",
                          "Sample_X05","Sample_X06","Sample_X07","Sample_X08",
                          "Sample_X09"))
df_renamed<-FindReplace(data=df, Var="SampleID", replaceData=replaces, 
                 from="from", to="to", exact=TRUE)

#change back to df for simplicity:
df2<-df_renamed

#for multiple runs, name your run with the date (MMDDYYYY) so you can eventually
#merge with the corresponding metadata

#get rid of blank wells in the plate, if you have them
df3<- 
  df2%>%
  filter(SampleID != "Blank_B")

#per Haley: 
#####mean & sd fluorescence for negative control wells#####
#cycles 2 to 8, as well as standard deviation across cycles (calc 2 ways):

neg<-df3%>% 
  filter(SampleID == "Negative_control_N")

neg_means<-df3%>% 
  filter(SampleID == "Negative_control_N")%>%
  select(t2,t3,t4,t5,t6,t7,t8)%>%
  gather()%>%
  summarize(min=min(value),
            max=max(value),
            median=median(value),
            mean=mean(value),
            sd=sd(value))

neg_mean<-mean(as.matrix(neg[,c(4:10)]))
neg_sd<-sd(as.matrix(neg[,c(4:10)]))

print(c(neg_mean, neg_means$mean))
print(c(neg_sd, neg_means$sd))

all_means<-df2%>% 
  select(t2,t3,t4,t5,t6,t7,t8)%>%
  gather()%>%
  summarize(min=min(value),
            max=max(value),
            median=median(value),
            mean=mean(value),
            sd=sd(value))

threshold_all<-all_means$mean+10*all_means$sd
threshold_all

max_fl<-max(df3[,4:99])

#so these are the same...gooood. but still I'm not totally sure if we want
#the overall mean, or a mean of means & a mean of the SDs by well.
#I *think* because it says "across cycles", the the std dev and mean
#of all cells cumulatively from t2-t8 is what we want though.

####thresholds#####
#the "threshold" is the mean +10 standard dev, so:
threshold10<-neg_mean+(10*neg_sd)
threshold10
threshold10/max_fl

#another option: take first time point + 10 SD
threshold_t1<-(mean(df3$t1)+(10*sd(df3$t1)))
threshold_t1
threshold_t1/max_fl

#I'm not really sure whta the best way to calculate threshold is.
#I think for an R package, giving people the option to set which type of 
#threshold they want to use would be a nice option. they are all fairly close
#at least. I think I prefer using the average of all samples for t2-t8 +10 SD
#personally so maybe that will be the default. 

####calculating time to threshold for samples#####

df_samples<-df3[,1:99]%>%
  filter(SampleID !="Negative_control_N")

df_samples_melt<-reshape2::melt(df_samples, id.vars=c("SampleID","Well"))
df_samples_melt$timepoint<-str_remove(df_samples_melt$variable, "t")
df_samples_melt$timepoint<-as.integer(df_samples_melt$timepoint)

df_samples_melt$percent_total_fluorescence<-df_samples_melt$value/max_fl

df_samples_melt_avg<-
  df_samples_melt%>%
  group_by(SampleID, timepoint)%>%
  summarize(mean_fl=mean(percent_total_fluorescence),
            sd_fl=sd(percent_total_fluorescence)) 
  
threshold_percent<-threshold_all/max_fl

pal2 <- wes_palette(15, name="Zissou1",type = "continuous")



p0b<-ggplot(aes(x=timepoint, y=mean_fl, group=dfab_meta, color=as.factor(mean_elisa),
                shape=dfab_meta2), data=df_sums2)+geom_point()+
  geom_line(alpha=0.8)+theme_classic()+scale_color_manual(values=pal2)+
  geom_hline(yintercept=threshold_percent, lty=2, col="gray70")+
  annotate("text", label="threshold", x=4, y=0.22)+
  labs(color="ELISA test value")+ylab("percent total fluorescence")+xlab("cycle")
p0b

p0c<-ggplot(df_vals, aes(x=elisa))+geom_histogram(color="black",fill="gray70")+theme_classic()+
  xlab("ELISA test value")
p0c


my_pal<-rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4'))
pal <- wes_palette("Zissou1", 9, type = "continuous")

p0<-ggplot(aes(x=timepoint, y=mean_fl, group=SampleID, color=as.factor(no_pos), 
               shape=as.factor(no_pos)), data=df_samples_melt_avg)+geom_point()+
  geom_line(alpha=0.8)+theme_classic()+scale_color_manual(values=pal)+
  #geom_ribbon(aes(ymin=mean_fl-sd_fl, ymax=mean_fl+sd_fl), alpha=0.1)+
  scale_shape_manual(values=c(1, rep(19,8)), guide="none")+
  geom_hline(yintercept=threshold_percent, lty=2, col="gray70")+
  annotate("text", label="threshold", x=4, y=0.22)+
  labs(color="no. positives in pool")+ylab("percent total fluorescence")+xlab("cycle")
p0

p0+p0b+ plot_layout(widths=c(3,4))+plot_annotation(tag_levels="A")








df_over<-df_samples_melt%>%
  filter(value > threshold_all)

df_tt<-df_over%>%
  group_by(Well, RunDate)%>%
  slice(which.min(timepoint))

df_tt$timepoint_avg<-df_tt$timepoint-((df_tt$value-threshold_all)/threshold_all)

##convert "timepoint" to minutes, then to hours (subtract 1 because t1 = 0)
df_tt$minutes<-(df_tt$timepoint_avg-1)*15
df_tt$hours<-df_tt$minutes/60
df_tt$seconds<-df_tt$minutes*60

df_tt$amyloid_rate<-(1/df_tt$hours)

pos_ctrl_tts<-subset(df_tt, df_tt$SampleID=="Sample_X01")
pos_avg_ttt<-mean(pos_ctrl_tts$seconds)

df_tt$rel_rate<-pos_avg_ttt/df_tt$seconds

avg_tt<-df_tt%>%
  group_by(SampleID, RunDate)%>%
  summarize(avg_tt_hours=mean(hours))

meta$RunDate<-as.character(meta$RunDate)
meta$RunDate<-paste("0",meta$RunDate, sep="")
tt_merged<-merge(avg_tt, meta, by=c("RunDate","SampleID"), drop=TRUE)

##select wells that never reach threshold, set amyloid rate at zero:


#mean & std dev for plotting if you want:
#####PLOTZ#####
meta3<-merge(meta, df2[,c(1:2)], by="SampleID")
df_tt_merged<-merge(meta3[,c(1,3,6,7)], df_tt, by=c("SampleID","Well"), all=TRUE)
df_tt_merged[is.na(df_tt_merged)]<-0

df_tt_summary<-df_tt_merged%>%
  group_by(SampleID)%>%
  summarize(mean=mean(hours),
            sd=sd(hours))
df_tt_summary$concentration<-as.factor(df_tt_summary$concentration)
p1<-ggplot(aes(x=SampleID,color=type, y=mean, shape=concentration), data=df_tt_summary)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2)+
  geom_hline(yintercept=0.1, lty=2)+
  geom_hline(yintercept=0.05, lty=2)+
  annotate("text", label="10 hours", x="Sample_X30", y=0.105)+
  annotate("text", label="20 hours", x="Sample_X30", y=0.055)+
  ylab("amyloid formation rate (1/h)")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45, hjust=1))
  
p1


p1b<-ggplot(aes(x=SampleID, y=mean, color=type), data=df_tt_summary)+
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd))+
  coord_flip()+
  facet_wrap(~type, scales="free_y",ncol=1)+theme_classic()
p1b

#or change however you want and incorporate metadata. 

###you can also plot the regular output type curves:
df_ag<-aggregate(.~SampleID+RunDate, data=df2[,c(2:98,100)], mean)
df_melt<-reshape2::melt(df_ag, id.vars=list("SampleID","RunDate"))
df_melt$percent_total<-as.numeric(df_melt$value)/max_fl

threshold10/max
#now you have the average for all samples at all time points
#you can recreate the heinous figures from MARS like this:

df_merged<-merge(meta, df_melt, by=c("SampleID","RunDate"), drop=TRUE)
df_merged2<-subset(df_merged, df_merged$variable !="t1" & df_merged$type!="NA")
df_merged2$cycle<-str_remove(df_merged2$variable, "t")
df_merged2$cycle<-as.integer(df_merged2$cycle)

df_subset<-subset(df_merged2, df_merged2$LabID !="pos_02")
#df_subset<-subset(df_merged2, df_merged2$concentration==1.0 | df_merged2$concentration==0.1)

df_subset$IHC_status<-factor(df_subset$IHC_status, levels=c("ND","ISF","POS"))

df_ctrls<-df_merged2%>%
  filter(LabID == "pos_02")%>%
  group_by(SampleID, RunDate, concentration, type, IHC_status, variable, cycle)%>%
  summarize(percent_total=mean(percent_total),
            value=mean(value))
df_ctrls$LabID<-"X01"

df_ctrls2<-df_ctrls[,c(1,2,10,3:6,9,8,7)]
df_ctrls2$IHC_status<-"control"

df_new<-rbind(df_subset, df_ctrls2)
df_new$IHC_status<-factor(df_new$IHC_status, levels=c("control","ND","ISF","POS"))

df_new$concentration<-as.factor(df_new$concentration)
df_new$ID_conc_date<-paste(df_new$LabID,df_new$concentration,df_new$RunDate, sep="_")
df_new$hour<-as.numeric(df_new$cycle*15)/60


p0<-ggplot(aes(x=hour, y=percent_total, group=ID_conc_date, color=IHC_status, alpha=concentration), 
              data=df_new)+
  geom_point(size=0.5)+
  theme_classic()+geom_line()+geom_hline(yintercept=threshold_all/max_fl, lty=2, col="gray60", alpha=0.8)+
  ylab("proportion total fluorescence")+
 # annotate("text", label="threshold", x=2, y=(threshold_all/max_fl)+0.06)+
  facet_wrap(~IHC_status, ncol=2)+
  scale_x_continuous(breaks=c(0,6, 12, 18, 24))+
  scale_color_manual(values=c("gray40","#3B9AB2", "#EBCC2A","#F21A00"))+
  guides(alpha = guide_legend(reverse=T))
p0


dd<-data.frame(cbind(IHC_status, animal_ID, avg_tt, follicles_pos))
dd$follicles_pos<-as.numeric(dd$follicles_pos)
dd$avg_tt<-as.numeric(dd$avg_tt)

dd_mean<-dd%>%
  group_by(IHC_status)%>%
  summarize(avg_tt3=mean(avg_tt, na.rm=TRUE),
            sd_tt=sd(avg_tt, na.rm=TRUE))

p0b<-ggplot(aes(x=animal_ID, y=avg_tt3, group=animal_ID, color=IHC_status, shape=IHC_status), data=dd_mean)+
  geom_pointrange(aes(ymax=avg_tt3+sd_tt, ymin=avg_tt3-sd_tt), size=0.6)+
  scale_y_reverse()+
  scale_color_manual(values=c("#EBCC2A","#F21A00"))+
  ylab(paste("hours to threshold","\n","(mean ± SD)"))+
  xlab("original animal ID")+
  theme_classic()
p0b

p0/p0b+plot_layout(heights=c(3,1))+plot_annotation(tag_levels="A")

equal_breaks <- function(n = 3){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}

install.packages('scales')
library(scales)
p0+facet_wrap(~type, scales="free", ncol=1)+theme(legend.position="none")+
  scale_y_continuous(breaks=scales::pretty_breaks(n=5))

##eventually, merge with metadata & combine different runs etc. but you know.



