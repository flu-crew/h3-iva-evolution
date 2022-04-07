library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(reshape2)
library(scales)
library(patchwork)
library(viridis)
library(png)
library(stringr)

axis_size <- 9
axis_title_size <- 12
x_axis_angle <- 90
coeff <- 0.5 

#read in table and format dates
setwd("//iastate/lss/research/pcgauger-lab/Megan/manuscripts/IV_a Paper/h3-iva-evolution/data_collection")
df <- read.csv("swine-surveillance-data.txt", sep="\t", header = T)

#convert date to decimal
df$format_date <- as.Date(df$Date, "%Y-%m-%d")
df$decimal_date <- decimal_date(df$format_date)
df$year <- floor(df$decimal_date)

df <- df %>%
  mutate(HA_type = case_when(!is.na(H1) ~ "H1",!is.na(H3) ~ "H3")) %>%
  mutate(HA = coalesce(H1,H3,HA_type)) %>%
  filter(!is.na(HA)) %>%
  unite('HA_clade',HA_type,HA,sep="-",remove=FALSE, na.rm=TRUE) 
  
#df$HA <- as.factor(df$HA)
#levels(df$most.similar.blast.clade) <-  c("2010.1", "2010.2", "Cluster IV","Cluster IVA","Cluster IVB","Other","Other","Other","Other","Other","Other") 
#df$most.similar.blast.clade <- factor(df$most.similar.blast.clade, levels = c("2010.1", "2010.2","Cluster IV","Cluster IVA","Cluster IVB","Other") )

#get relative frequency of clades by year
df_agg <- df %>%
  group_by(year,HA_clade) %>%
  summarize(N = n()) %>%
  mutate(freq = N / sum(N)) %>%
  ungroup() %>%
  complete(year, HA_clade, fill=list(N=0,freq=0))

df_agg$year <- as.Date(ISOdate(df_agg$year,1,1))

#relative frequency area graph
plot1 <- ggplot(df_agg, aes(x=year,y=freq * 100, fill=HA_clade)) +
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
    legend.text = element_text(size=axis_size)
  ) 



#### HA-NA Heatmap ####
df <- df %>%
  mutate(N2 = replace(N2, N2 == "Human-like", "Hu-like")) %>%
  mutate(N1 = replace(N1, N1 == "Classical", "Classic")) %>%
  mutate(N1 = replace(N1, N1 == "Pandemic", "Pdm")) %>%
  mutate(NA_type = case_when(!is.na(N1) ~ "N1",!is.na(N2) ~ "N2")) %>%
  mutate(NA_clade = coalesce(N1,N2)) %>%
  filter(!is.na(NA_clade)) %>%
  unite('NA_clade',NA_type,NA_clade,sep="-",remove=FALSE, na.rm=TRUE) 

df_agg_heatmap <- df %>%
  group_by(HA_clade,NA_clade) %>%
  summarize(N = n()) %>%
  mutate(freq = N / sum(n())) %>%
  ungroup() %>%
  complete(HA_clade, NA_clade, fill=list(N=0,freq=0))

plot2 <- ggplot(df_agg_heatmap, aes(HA_clade,NA_clade,fill=N)) +
  geom_tile() +
  geom_text(aes(label = N),size=2) +
  scale_fill_gradient(low="white", high="blue") +
  labs(y="NA Clade", x="HA Clade", fill="Detections") +
  theme( panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"),
        axis.title = element_text(size=axis_title_size,face="bold"),
        axis.text = element_text(size=axis_size,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        legend.title = element_text(size=axis_title_size,face="bold"),
        legend.text = element_text(size=axis_size),
  )


### HA-NA-Constellation ###

df <- df %>%
  mutate(Constellation = str_replace_all(Constellation, "-", "~")) %>%
  mutate(Constellation = replace(Constellation, N2 == "Human-like", "Hu-like")) %>%
  unite('HA_NA_Constellation',HA_clade,NA_clade,Constellation,sep="/",remove=FALSE, na.rm=TRUE)

#df <- df %>%
#  mutate(HA_NA_Constellation = replace(HA_NA_Constellation, grepl("~", HA_NA_Constellation), "Other"))

df_r <- df %>%
  filter(WGS == 'TRUE')

df_agg_2 <- df_r %>%
  group_by(year,HA_NA_Constellation) %>%
  summarize(N = n()) %>%
  mutate(freq = N / sum(N)) %>%
  ungroup() %>%
  #filter(freq > 0.05) %>%
  mutate(HA_NA_Constellation = replace(HA_NA_Constellation, freq <= 0.05, "Other")) %>%
  filter(year > 2010)
  #filter(HA_NA_Constellation != "Other")

df_agg_2$HA_NA_Constellation <- as.factor(df_agg_2$HA_NA_Constellation)
df_agg_2$HA_NA_Constellation <- relevel(df_agg_2$HA_NA_Constellation, "Other")
df_agg_2$year <- as.factor(df_agg_2$year)

plot3 <- ggplot(df_agg_2, aes(fill=HA_NA_Constellation, y=freq, x=year)) +
  geom_bar(position="stack", stat = "identity") +
  labs(y="Relative Freq. of Detection", x="Year", fill="HA Clade/NA Clade/Constellation") +
  guides(fill=guide_legend(ncol=2)) +
  theme( 
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "gray90"),
    panel.grid.minor = element_line(colour = "gray90"),
    axis.title = element_text(size=axis_title_size,face="bold"),
    axis.text = element_text(size=axis_size,face="bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size=axis_size),
    legend.position="bottom"
  )

#export plots as .tiff, combine in ppt
setwd("//iastate/lss/research/pcgauger-lab/Megan/manuscripts/IV_a Paper/h3-iva-evolution/figures")

tiff("Supp_Figure_2c_btm.tiff", units = "in", width = 5, height = 8.5, res = 300, compression = "lzw")
plot3
dev.off()

tiff("Supp_Figure_2ab_btm.tiff", units = "in", width = 6, height = 8.5, res = 300, compression = "lzw")
plot1 / plot2
dev.off()

#try to format plots with patchwork, couldn't get to work well
tiff("Supp_Figure_2_try.tiff", units = "in", width = 20, height =15, res = 300, compression = "lzw")

#plot1+plot_spacer()+plot2 + plot_layout(width=c(1.0,0.01,1.0)) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"))
#plot1 / plot2 + plot_layout(width=c(1.0,0.01,1.0)) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"))
#(plot1 / plot2) | plot3 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"))
plot1 + plot2 + plot3 + plot_layout(ncol = 2) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face="bold"))

dev.off()
