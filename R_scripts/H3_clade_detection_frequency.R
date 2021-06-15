library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(reshape2)
library(scales)
library(patchwork)

axis_size <- 17
axis_title_size <- 20
x_axis_angle <- 45
coeff <-5 

setwd("//iastate/lss/research/pcgauger-lab/Megan/IV_a Paper/c-iva/data_collection/HA/")

#read in table and format dates
df <- read.csv("classification_for_detection_frequency.csv", header = T)
#df <- df[df$most.similar.blast.clade != "gamma",]
df$format_date <- as.Date(df$date, "%m/%d/%Y")

#deal with year only dates
for(i in 1:length(df$format_date)) {
  if(is.na(df$format_date[i])) {
    df$format_date[i] <- as.Date(df$date[i], "%Y")
    print(as.Date(df$date[i], "%Y"))
  }
}

#deal with month-year only dates
for(i in 1:length(df$format_date)) {
  if(is.na(df$format_date[i])) {
    df$format_date[i] <- as.Date(parse_date_time(df$date[i], orders = "my"))
    #print(as.Date(df$date[i], "%b-%y"))
  }
}

#convert date to decimal
df$decimal_date <- decimal_date(df$format_date)
df$year <- floor(df$decimal_date)
df$most.similar.blast.clade <- as.factor(df$most.similar.blast.clade)

levels(df$most.similar.blast.clade) <-  c("2010.1", "2010.2", "Cluster IV","Cluster IVA","Cluster IVB","Other","Other","Other","Other","Other","Cluster I") 
df$most.similar.blast.clade <- factor(df$most.similar.blast.clade, levels = c("2010.1", "2010.2","Cluster I","Cluster IV","Cluster IVA","Cluster IVB","Other") )

#get relative frequency of clades by year
df_agg <- df %>%
  group_by(year,most.similar.blast.clade) %>%
  summarize(N = n()) %>%
  mutate(freq = N / sum(N)) %>%
  ungroup() %>%
  complete(year, most.similar.blast.clade, fill=list(N=0,freq=0))

df_agg$year <- as.Date(ISOdate(df_agg$year,1,1))

#relative frequency line graph
#(df_agg, aes(x=year,y=freq, color=most.similar.blast.clade)) +
#  geom_line(size=1.1) +
#  labs(x= "Year", y = "Percent")

#relative frequency area graph
ggplot(df_agg, aes(x=year,y=freq * 100, fill=most.similar.blast.clade)) +
  geom_area(color="white") +
  scale_x_date(labels = date_format("%Y"), breaks=date_breaks("1 year")) +
  labs( x = "Time", y ="H3 Clade Frequency of Detection - By Proportion", fill = "Clade")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "gray90"),
    panel.grid.minor = element_line(colour = "gray90"),
    axis.text = element_text(size=axis_size, face="bold"),
    axis.text.x = element_text(angle=x_axis_angle,vjust=0.5),
    axis.title = element_text(size=axis_title_size, vjust = -0.35, face="bold"),
    legend.title = element_text(size=20, face="bold"),
    legend.text = element_text(size=20, face="bold")
  ) 


#### EPS Graph ####
setwd("//iastate/lss/research/pcgauger-lab/Megan/civa/civa_diversity_result/")

#get C-IVA detection frequency from dataframe above
df_civa <- df[df$most.similar.blast.clade == "Cluster IVA",]
dataDetFreq <- df_civa %>%
  group_by(year) %>%
  summarize(Count = n())
colnames(dataDetFreq) <- c("Time","Count")

#get VDL detection frequency
#vdlDetFreq <- read.csv("cIVa_VDL_detection_frequency_updated.csv", header=TRUE)
#vdlDetFreq$Time <- decimal_date(as.Date(vdlDetFreq$Time))
#vdlDetFreq$Time <- floor(vdlDetFreq$Time)
#vdlDetFreqYear <- aggregate(vdlDetFreq$Count, by = list(Category=vdlDetFreq$Time), FUN=sum)
#colnames(vdlDetFreqYear) <- c("Time", "Count")
#vdlDetFreq <- vdlDetFreqYear

#get IRD EPS data from csv
dataIRD <- read.csv("run1_ird_results/civa.csv", header=TRUE, sep="\t")

#Graph with overlay
ggplot(dataIRD, aes(x = Time, y = Median)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.5) +
  geom_line(aes(y = Median, color = "Median EPS"), size = 1.5) + 
  geom_line(data = dataDetFreq , aes(y= Count/coeff, color = 'Public Detection'), size = 2) +
  #geom_line(data = vdlDetFreq, aes(y=Count/coeff,color='VDL Detection'), size = 2) +
  scale_y_continuous(
    name = "Effective Population Size (EPS)",
    sec.axis = sec_axis(~.*coeff, name = "Detection Frequency")
  ) +
  scale_x_continuous(
    name = "Time",
    breaks = c(2011, 2012, 2013,2014,2015,2016,2017,2018,2019,2020), 
    limits = c(2010.565,2020.213)
  ) +
  scale_color_manual(values = c(
    'Median EPS' = 'darkblue',
    'Public Detection' = 'darkorange',
    'VDL Detection' = 'darkred')) + 
  labs(y = "Effective Population Size (EPS)", x = "Time", color = "") +
  theme(legend.position = "bottom", 
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"),
        title = element_text(size=25,face="bold"),
        axis.title = element_text(size=axis_title_size,face="bold"),
        axis.text = element_text(size=axis_size,face="bold"),
        axis.text.x = element_text(angle=x_axis_angle, vjust = 0.5),
        axis.title.x = element_text(vjust=-1.0),
        axis.title.y.right = element_text(vjust=1.5,hjust=0.5),
        legend.title = element_text(size=25,),
        legend.text = element_text(size=20, face="bold"),
  )
