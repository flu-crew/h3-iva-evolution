library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(reshape2)
library(scales)
library(patchwork)
library(viridis)

axis_size <- 17
axis_title_size <- 18
x_axis_angle <- 90
coeff <- 0.5 

setwd("//iastate/lss/research/pcgauger-lab/Megan/manuscripts/IV_a Paper/h3-iva-evolution/data_collection/HA/")

#read in table and format dates
df <- read.csv("classification_for_detection_frequency.csv", header = T)
#df <- df[df$most.similar.blast.clade != "gamma",]
df$format_date <- as.Date(df$date, "%m/%d/%Y")

#deal with year only dates
for(i in 1:length(df$format_date)) {
  if(is.na(df$format_date[i])) {
    df$format_date[i] <- as.Date(df$date[i], "%Y")
  }
}

#deal with month-year only dates
for(i in 1:length(df$format_date)) {
  if(is.na(df$format_date[i])) {
    df$format_date[i] <- as.Date(parse_date_time(df$date[i], orders = "my"))
  }
}

#convert date to decimal
df$decimal_date <- decimal_date(df$format_date)
df$year <- floor(df$decimal_date)
df$most.similar.blast.clade <- as.factor(df$most.similar.blast.clade)

levels(df$most.similar.blast.clade) <-  c("2010.1", "2010.2", "Cluster IV","Cluster IVA","Cluster IVB","Other","Other","Other","Other","Other","Other") 
df$most.similar.blast.clade <- factor(df$most.similar.blast.clade, levels = c("2010.1", "2010.2","Cluster IV","Cluster IVA","Cluster IVB","Other") )

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
plot3 <- ggplot(df_agg, aes(x=year,y=freq * 100, fill=most.similar.blast.clade)) +
  geom_area(color="white") +
  scale_x_date(labels = date_format("%Y"), breaks=date_breaks("1 year"),
               limits = as.Date(c('2011-01-01','2021-03-01'))) +
  labs( x = "Time", y ="Relative Freq. of Detection", fill = "Clade")+
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "gray90"),
    panel.grid.minor = element_line(colour = "gray90"),
    axis.text = element_text(size=axis_size, face="bold"),
    axis.text.x = element_text(angle=x_axis_angle,vjust=0.5),
    axis.title.x = element_blank(),
    #axis.title.x = element_text(vjust=-1.0),
    axis.title = element_text(size=axis_title_size, vjust = -0.35, face="bold"),
    legend.title = element_blank(),
    #legend.title = element_text(size=axis_title_size, face="bold"),
    legend.text = element_text(size=axis_size, face="bold")
  ) 

#### HA-NA Detection Frequency Graph ####
clades <- read.csv("//iastate/lss/research/pcgauger-lab/Megan/manuscripts/IV_a Paper/h3-iva-evolution/data_collection/Figure1/ha-na_by_year_for_R.csv", na.strings = 0)
clades[is.na(clades)] <- 0
longClades <- melt(clades, id = c("year"))

normClades <- clades/rowSums(clades[2:7])
normClades[1] <- clades[1]
normLongClades <- melt(normClades, id = c("year"))

normLongClades$variable <- factor(normLongClades$variable, levels = c("H3_CIVA_2002A","H3_CIVA_2002B","H3_2010.1_2002A",
                                                                      "H3_2010.1_2002B","H3_CIVA_Other","H3_2010.1_Other"))

plot1 <- ggplot(data = normLongClades, aes(x=year, y=value, color=variable, fill=variable))+
  geom_area(alpha=1, size=.5, colour="white") +
  #pretty good colors
  #scale_fill_viridis(discrete = T, option="plasma", name="") +
  scale_fill_viridis(discrete = T, option="inferno", name="") +
  labs(y = "Proportion of Detection", color="test") +
  scale_x_continuous(
    name = "Time",
    breaks = c(2011, 2012, 2013,2014,2015,2016,2017,2018,2019,2020,2021),
    limits = c(2011,2021.230)
  ) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(legend.position = "bottom", 
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"),
        axis.title = element_text(size=axis_title_size,face="bold"),
        axis.text = element_text(size=axis_size,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        legend.title = element_text(size=axis_title_size,),
        legend.text = element_text(size=axis_size,face="bold"),
  )


#### EPS Graph ####
setwd("//iastate/lss/research/pcgauger-lab/Megan/manuscripts/IV_a Paper/h3-iva-evolution/beast_analysis/")

#get C-IVA detection frequency from dataframe above
df_civa <- df[df$most.similar.blast.clade == "Cluster IVA",]
df_civa$month <- month(df_civa$format_date)

df_civa$my <- paste(df_civa$monthly, df_civa$year)

dataDetFreq <- df_civa %>%
  group_by(month, year) %>%
  #group_by(year) %>%
  summarize(Count = n())
colnames(dataDetFreq) <- c("Month","Year", "Count")

dataDetFreq <- dataDetFreq %>%
  group_by(Year) %>%
  summarize(Count = mean(Count))
dataDetFreq$Time <- dataDetFreq$Year


#get VDL detection frequency
#vdlDetFreq <- read.csv("cIVa_VDL_detection_frequency_updated.csv", header=TRUE)
#vdlDetFreq$Time <- decimal_date(as.Date(vdlDetFreq$Time))
#vdlDetFreq$Time <- floor(vdlDetFreq$Time)
#vdlDetFreqYear <- aggregate(vdlDetFreq$Count, by = list(Category=vdlDetFreq$Time), FUN=sum)
#colnames(vdlDetFreqYear) <- c("Time", "Count")
#vdlDetFreq <- vdlDetFreqYear

#get IRD EPS data from csv
dataIRD <- read.csv("Results/c-iva_sequences_subset_clean_final_alignment_gmrf_skyride_reconstruction.log.tsv", header=TRUE, sep="\t")

#Graph with overlay
plot2 <- ggplot(dataIRD, aes(x = Time, y = Median)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.5) +
  geom_line(aes(y = Median, color = "Median EPS"), size = 1.5) + 
  geom_line(data = dataDetFreq , aes(y= Count/coeff, color = 'Detection'), size = 2) +
  scale_y_continuous(
    name = "Effective Population Size (EPS)",
    sec.axis = sec_axis(~.*coeff, name = "Avg. Detection Freq. per Month")
  ) +
  scale_x_continuous(
    name = "Time",
    breaks = c(2011, 2012, 2013,2014,2015,2016,2017,2018,2019,2020,2021),
    limits = c(2011,2021.230)
  ) +
  scale_color_manual(values = c(
    'Median EPS' = 'darkblue',
    'Detection' = 'darkorange',
    'VDL Detection' = 'darkred')) + 
  labs(y = "Effective Population Size (EPS)", x = "Time", color = "") +
  theme(legend.position = "bottom", 
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"),
        title = element_text(size=25,face="bold"),
        axis.title = element_text(size=axis_title_size,face="bold"),
        axis.text = element_text(size=axis_size,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=x_axis_angle, vjust = 0.5),
        #axis.title.x = element_text(vjust=-1.0),
        axis.title.y.right = element_text(vjust=1.5,hjust=0.5),
        legend.title = element_text(size=axis_title_size,),
        legend.text = element_text(size=axis_size, face="bold"),
  )

setwd("//iastate/lss/research/pcgauger-lab/Megan/manuscripts/IV_a Paper/h3-iva-evolution/figures/")

tiff("Figure_1_HA-NA_detection_and_demography.tiff", units = "in", width = 8.5, height = 11, res = 300, compression = "lzw")

plot1+plot_spacer()+plot2 + plot_layout(width=c(1.0,0.01,1.0)) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"))

plot1 / plot2 + plot_layout(width=c(1.0,0.01,1.0)) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"))


dev.off()
