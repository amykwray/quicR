###simulated dataset: there's probably a slicker way to do this but I'm not that good at for loops!
###creates simulation of RT-QuIC data for half of a plate, where all controls & samples have 3 replicates
###by AK Wray 21 April 2023

setwd("C:/Users/awray/Desktop/quicR")


# Load required packages
library(tidyverse)
library(reshape2)

# Set parameters
n_cycles <- 100  # Number of RT-QuIC cycles
n_replicates_pos<-3
n_replicates_pos_samples<-8
n_replicates_neg<-28
max_fluor <- 150000 # Maximum fluorescence value
min_fluor<-6000
plateau <- 150000  # Plateau fluorescence value
slope <- 0.10    # Slope of sigmoidal function
noise_sd <- 6000   # Standard deviation of noise

##### create data for positive controls with dilution series ######
####the goal here is to create simulated data with a fair amount of noise
####also, as pos ctrls are more "diluted", replicates should have more variability

data1 <- matrix(NA, nrow = n_replicates_pos, ncol = n_cycles)
for (j in 1:n_replicates_pos) {
  cycle <- seq_len(n_cycles)
  fluor <- min_fluor + (3*j/3.003) * max_fluor * (1 / (1 + exp(-slope * (cycle - (n_cycles / 2.8))))) # Sigmoidal function
  fluor <- fluor + rnorm(n_cycles, 0, noise_sd)                           # Add noise
  fluor[fluor > plateau] <- plateau      # Apply plateau effect
  data1[j,]<-fluor
}

data2<-matrix(NA, nrow = n_replicates_pos, ncol = n_cycles)
for (j in 1:n_replicates_pos) {
  cycle <- seq_len(n_cycles)
  fluor <- min_fluor +  (3*j/3.05) * max_fluor * (1 / (1 + exp(-slope * (cycle - (n_cycles / 1.8))))) # Sigmoidal function
  fluor <- fluor + rnorm(n_cycles, 0, noise_sd)                             # Add noise
  fluor[fluor > plateau] <- plateau      # Apply plateau effect
  data2[j,]<-fluor
}

data3 <- matrix(NA, nrow = n_replicates_pos, ncol = n_cycles)
for (j in 1:n_replicates_pos) {
  cycle <- seq_len(n_cycles)
  fluor <- min_fluor +  (3*j/3.1) * max_fluor * (1 / (1 + exp(-slope * (cycle - (n_cycles / 1.2))))) # Sigmoidal function
  fluor <- fluor + rnorm(n_cycles, 0, noise_sd*1.1)                             # Add noise
  fluor[fluor > plateau] <- plateau      # Apply plateau effect
  data3[j,]<-fluor
}

data4<-matrix(NA, nrow = n_replicates_pos, ncol = n_cycles)
for (j in 1:n_replicates_pos) {
  cycle <- seq_len(n_cycles)
  fluor <- min_fluor +  (3*j/3.4) * max_fluor * (1 / (1 + exp(-slope * (cycle - (n_cycles / 1.001))))) # Sigmoidal function
  fluor <- fluor + rnorm(n_cycles, 0, noise_sd*1.2)                             # Add noise
  fluor[fluor > plateau] <- plateau      # Apply plateau effect
  data4[j,]<-fluor
}

pos_ctrls<-rbind(data1,data2,data3,data4)
pos_ctrls[pos_ctrls<min_fluor]<-min_fluor
wells<-c("A01","A02","A03","A04","A05","A06","A07","A08","A09","A10","A11","A12")
names<-c(rep("pos_ctrl_1",12))
type<-c(rep("pos_ctrl",12))

concentration<-c(rep("-4",3), rep("-5",3), rep("-6",3),rep("-7",3))
pos_ctrls_names<-as.data.frame(cbind(wells,names,type, concentration, pos_ctrls))

####create data for negative controls & presumed nonamplifying samples #####
#let's say you have 13 samples and 2 negative controls, each with 3 replicates
#of those samples, 3 are positive (2 for which all 3 replicates are positive, and 1 for which 2/3 replicates are positive)

data5 <- matrix(NA, nrow = n_replicates_neg, ncol = n_cycles)
for (j in 1:n_replicates_neg) {
  cycle <- seq_len(n_cycles)
  fluor <- min_fluor  # Sigmoidal function
  fluor <- fluor + rnorm(n_cycles, 0, noise_sd/2)   # Add noise
  fluor[fluor > plateau] <- plateau      # Apply plateau effect
  data5[j,]<-fluor
}

##now, assign well numbers. lets say, B10-12, C04-C06, and D01-D02 are the positives (leave out)
wells2<-c("B01","B02","B03","B04","B05","B06","B07","B08","B09",
          "C01","C02","C03","C07","C08","C09","C10","C11","C12",
          "D03","D04","D05","D06","D07","D08","D09","D10","D11","D12",
          "B10","B11","B12","C04","C05","C06","D01","D02")

###create data for presumed positive samples #####
#2 positive triplets + 2 positive/3 replicates:

data6<-matrix(NA, nrow = n_replicates_pos_samples, ncol = n_cycles)
for (j in 1:n_replicates_pos_samples) {
  cycle <- seq_len(n_cycles)
  fluor <- min_fluor +  (8*j/8.3) * max_fluor * (1 / (1 + exp(-slope * (cycle - (n_cycles / 1.001))))) # Sigmoidal function
  fluor <- fluor + rnorm(n_cycles, 0, noise_sd*1.2)                             # Add noise
  fluor[fluor > plateau] <- plateau      # Apply plateau effect
  data6[j,]<-fluor
}

all_samples<-as.data.frame(rbind(data5,data6))
all_samples_named<-as.data.frame(cbind(wells2, all_samples))
all_samples_ordered<-with(all_samples_named,  all_samples_named[order(wells2) , ])

names2<-c(rep("sample_01",3), rep("sample_02",3), rep("sample_03",3),rep("sample_04",3), 
          rep("sample_05",3), rep("sample_06",3),rep("sample_07",3), rep("sample_08",3), 
          rep("sample_09",3), rep("sample_10",3), rep("neg_ctrl_1",3), rep("neg_ctrl_2",3))

type2<-c(rep("sample",30),rep("neg_ctrl",6))
concentration2<-paste(rep("-4", 36))

all_samples_ordered[all_samples_ordered < min_fluor] <- min_fluor
all_samples_combo<-cbind(names2, concentration2, type2, all_samples_ordered)
all_samples_final<-as.data.frame(all_samples_combo[,c(4,1,3,2,5:104)])

######combine all simulated data#####
colnames(all_samples_final)<-c("well","names","type", "concentration",paste(c(1:100)))
colnames(pos_ctrls_names)<-c("well","names","type", "concentration",paste(c(1:100)))

fake_data<-rbind(pos_ctrls_names, all_samples_final)
df<-fake_data

df_melt<-reshape2::melt(df, id.vars=c("well", "names","type", "concentration"))
df_melt$variable<-as.numeric(df_melt$variable)
df_melt$value<-as.numeric(df_melt$value)

p1<-ggplot(aes(x=variable, y=value, color=as.factor(type), lty=concentration), data=df_melt)+
  geom_point()+geom_line(aes(group=as.factor(wells)))
p1+theme_classic()

#separate data "output" from metadata
df2<-df[,-c(2:4)]
meta<-df[,c(1:4)]
meta$run_ID<-paste("run1")
meta$tissue_type<-paste("RLN")
meta$species<-paste("elk")

#write.csv(df2, file="RTQ_simulated_data.csv", row.names=FALSE)
#write.csv(meta, file="RTQ_simulated_metadata.csv", row.names=FALSE)
