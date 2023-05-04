# quicR

## Getting Started
### Exporting data from BMG FLUOstar omega microplate reader & similar machines: 
After the RT-QuIC run, export the raw data using the "excel report" button. Do not include the run info. Do not transpose. To use this script, make sure that data are in the format where the COLUMN A is wells, COLUMN B is "content" (for example, whatever you named the samples in the machine's plate layout), with ROW 1 showing labels (e.g., 1C should say something like ###"raw data" then the filter info). ROW 2 should show the time period -- it doesn't really matter how it is formatted. You will also need a metadata file that corresponds to your Sample IDs somehow.


```{r echo=FALSE, message=FALSE, warning=FALSE}
setwd("C:/Users/awray/Desktop/quicR")
library(reshape2)
library(tidyverse)
library(data.table)
library(wesanderson)
library(knitr)
```

### Cleaning up your dataset:
To start off, you will need 2 data frames: 1 with the raw output from the plate reader, for which the first column should identify each well. You will also need a metadata file that relates your metadata information to each well. this MUST contain the exact format in the exact order where column names must start with "well", "names" (sample names, any unique identifier at the SAMPLE level [replicates are ok]), "type" must be one of the following: pos_ctrl, sample, neg_ctrl, or blank. "concentration" must include the dilution factor from the original tissue (so if the ORIGINAL tissue aliquot is diluted 1:1000, the dilution factor would be -3, then if you diluted that 1:10 when adding to wells it would be -4). Finally, you MUST include "run_ID", which can be any unique identifier that you use to differentiate a unique plate or run. Other metadata (species, tissue type, animal ID, etc) may be included in later rows. 

```{r echo=TRUE}
df<-data.frame(read.csv(file="RTQ_simulated_data.csv"))
meta<-data.frame(read.csv(file="RTQ_simulated_metadata.csv"))
```

```{r echo=FALSE}
dim(df)
```
input information about your run (either how many cycles you ran or how many cycles you want to cutoff to):
```{r echo=TRUE}
n_cycles<-100
n_neg_ctrl_reps<-6
max_fl<-150000
meta_vars<-ncol(meta)
```

merge with your metadata:
```{r echo=TRUE}
df2<-merge(meta[,1:5], df, by="well")
```

rename columns if they are not named EXACTLY as "well", "sample_ID","type","concentration","run_ID"
```{r echo=TRUE}
names(df2)[1:5]<-c("well", "sample_ID", "type", "concentration", "run_ID")
```

I personally do not like numbers in the labels for a column, so we are pasting a "t" for timepoint, even though it will be removed later:

```{r echo=TRUE}
names(df2)[6:(n_cycles+5)]<-paste("t", c(1:n_cycles), sep="")
```

### Creating plots
Plotting a typical fluorescence vs. time plot:

```{r echo=TRUE}
df_melt<-df2%>%
 pivot_longer(cols= -c(well, sample_ID, type, concentration, run_ID), names_to="timepoint",
              values_to="raw_fluor")
df_melt$percent_fluor<-df_melt$raw_fluor/max_fl
df_melt$timepoint<-str_remove(df_melt$timepoint, "t")
df_melt$timepoint<-as.integer(df_melt$timepoint)
```
Plot raw fluorescence:
```{r echo=TRUE}
raw_fluor_plot<-ggplot(aes(x=timepoint, y=as.numeric(raw_fluor), color=as.factor(type), lty=as.factor(concentration)), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(well)))+ theme_classic()+
  labs(x="cycle",y="raw fluorescence", color="type",lty="concentration")
```

  
```{r echo=FALSE}
raw_fluor_plot
```
![image](https://user-images.githubusercontent.com/70864453/236318197-e67fefc8-01ac-4733-ab59-06ca5965c6b5.png)

Plot percent fluorescence:
```{r echo=TRUE}
percent_fluor_plot<-ggplot(aes(x=as.integer(timepoint), y=percent_fluor, color=as.factor(type), lty=as.factor(concentration)), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(well)))+ theme_classic()+
  labs(x="cycle",y="percent fluorescence", color="type",lty="concentration")
percent_fluor_plot
```

### Determining thresholds
Calculate mean & sd fluorescence for negative control wells from cycles 2 to 8, as well as standard deviation across cycles:

```{r echo=TRUE}
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
```

Another option is to calculate this same value for ALL samples from cycle 2-8:
(the rationale for ignoring t1 is that sometimes t1 readings are artificially high, then drops in the following cycles)

```{r echo=TRUE}
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
```

Yet another option that other publications have used: take first time point + 10 SD:
```{r echo=TRUE}
threshold_t1<-(mean(df2$t1)+(10*sd(df2$t1)))
threshold_t1
```

### Plot fluorescence data with the threshold of your choice:

```{r echo=TRUE}
raw_fluor_plot<-ggplot(aes(x=as.integer(timepoint), y=raw_fluor, color=as.factor(type), lty=as.factor(concentration)), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(well)))+ theme_classic()+
  labs(x="cycle",y="raw fluorescence", color="type",lty="concentration")+
  geom_hline(yintercept=threshold_all10, lty=2, col="gray70")+
  annotate("text", label="threshold", x=4, y=threshold_all10+(threshold_all10/20))
raw_fluor_plot
```
![image](https://user-images.githubusercontent.com/70864453/236318265-51ecf31c-4afa-4342-b7e5-3f16448f4604.png)


### Plot percent fluorescence:

```{r echo=TRUE}
percent_fluor_plot<-ggplot(aes(x=as.integer(timepoint), y=percent_fluor, color=as.factor(type), lty=as.factor(concentration)), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(well)))+ theme_classic()+
  labs(x="cycle",y="percent fluorescence", color="type",lty="concentration")+
  geom_hline(yintercept=threshold_all10/max_fl, lty=2, col="gray70")+
  annotate("text", label="threshold", x=4, y=threshold_all10/max_fl+(threshold_all10/(max_fl*20)))
percent_fluor_plot
```

I'm not really sure what the best way to calculate threshold is. I think for the R package, giving people the option to set which type of threshold they want to use would be a nice option. I think I prefer using the average of all samples for t2-t8 +10 SD personally so maybe that will be the default. 

## Calculating time to threshold for samples

```{r echo=TRUE}
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
```
![image](https://user-images.githubusercontent.com/70864453/236318005-ba773232-5975-4958-98c6-542e5cc89f3e.png)

### Determine ranges for samples that cross threshold
Note that here, we are selecting replicates that cross the threshold, so samples with inconsistent replicates are not accounted for. I'm not sure if there is currently a preferred way to calculate the average time to threshold for samples where some replicates never cross the threshold, but keep in mind for data interpretation to use caution with inconsistent within-sample replicates.
```{r echo=TRUE}
df_over<-df_samples_melt%>%
  filter(raw_fluor > threshold_all10)

df_tt<-df_over%>%
  group_by(well)%>%
  slice(which.min(timepoint))

df_tt$timepoint_avg<-df_tt$timepoint-((df_tt$raw_fluor-threshold_all10)/threshold_all10)
```

### Convert "timepoint" to minutes, then to hours. For this, you will need to know the intervals at which fluorescence readings were taken & convert that into minutes. 
Subtract 1 because t1 = 0
```{r echo=TRUE}
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
kable(avg_tt)
```

### Plot ranges of samples that cross threshold:
```{r echo=FALSE}
avg_tt<-as.data.frame(avg_tt)
```

```{r echo=TRUE}
q<-ggplot(aes(x=as.factor(sample_ID), y=avg_tt_hours, color=as.factor(sample_ID)), data=avg_tt)+geom_point()+ geom_pointrange(mapping=aes(ymin=avg_tt_hours-sd_tt_hours, ymax=avg_tt_hours+sd_tt_hours), data=avg_tt, size=1)+ labs(color="Sample ID") + ylab("hours to threshold (mean \u00b1 SD) ")+xlab("sample")+theme_classic()
q
```
![image](https://user-images.githubusercontent.com/70864453/236318122-5dcd5adf-d127-47f3-9020-09c285cf2d4d.png)
