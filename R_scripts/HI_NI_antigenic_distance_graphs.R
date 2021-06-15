library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(patchwork)
library(colorspace)

min_text_size = 15
min_point_size = 5

### HI Graph ###

data <- readxl:: read_excel("//iastate/lss/research/pcgauger-lab/Megan/civa/civa_antigenic_study/41C/ACMACS/Output/60187e1c330274d831bf6e04-all-distances.xlsx")

data <- data %>%
  rename("Antigen1" = "Name1") %>%
  rename("Antigen2" = "Name2") %>%
  mutate(Antigen1 = factor(Antigen1)) %>%
  filter(Antigen1 == "A/SWINE/OKLAHOMA/A01770191/2015 MDCK2" | 
           Antigen1 == "A/SWINE/NORTH CAROLINA/A02245294/2019 MDCK2" |
           Antigen1 == "A/SWINE/MINNESOTA/A02266068/2018 MDCK2"
           ) %>%
  mutate(Antigen2 = factor(Antigen2)) %>%
  filter(Antigen2 == "A/SWINE/NEW YORK/A01104005/2011 MDCK4" |
           Antigen2 == "A/SWINE/IOWA/A01480656/2014 MDCK3" |
           Antigen2 == "A/SWINE/NEBRASKA/A01567651/2015 MDCK4" |
           Antigen2 == "A/SWINE/MINNESOTA/A01280592/2013 MDCK2" |
           #Antigen2 == "A/BEIJING/32/1992 MDCK2" |
           #Antigen2 == "A/WUHAN/359/1995 MDCK2" |
           Antigen2 == "A/SWINE/NORTH CAROLINA/A02245294/2019 MDCK2" |
           Antigen2 == "A/SWINE/WYOMING/A01444562/2013 MDCK2" |
           Antigen2 == "A/SWINE/INDIANA/A01202866/2011 MDCK3" |
           
           Antigen2 == "A/SWINE/OKLAHOMA/A01770191/2015 MDCK2" |
           Antigen2 == "A/SWINE/MINNESOTA/A02266068/2018 MDCK2" | 
           
           #Antigen2 == "A/SWINE/WYOMING/A01397400/2013 MDCK3" |
           #Antigen2 == "A/SWINE/OKLAHOMA/A01776961/2016 MDCK3" |
           Antigen2 == "A/SWINE/NORTH CAROLINA/A01730588/2016 MDCK3" |
           Antigen2 == "A/SWINE/IOWA/A01432826/2013 MDCK3" 
           #Antigen2 == "A/SWINE/IOWA/A01290671/2013 MDCK3" 
           #Antigen2 == "A/SWINE/ILLINOIS/A01475273/2014 MDCK3" |
           #Antigen2 == "A/SWINE/ILLINOIS/A01503899/2014 MDCK3" |
           #Antigen2 == "A/SWINE/ILLINOIS/A01555638/2015 MDCK4"
            ) %>%
  mutate(Antigen1 = factor(Antigen1)) %>%
  mutate(Antigen2 = factor(Antigen2)) %>%
  group_by(Antigen1)

data$strain <- data$No2
data$xlab <- data$No1

data <- data %>%
  mutate(strain = case_when(No2 == "5 AG" ~ "BJ/92", 
                           No2 == "6 AG" ~ "WU/95", 
                           No2 == "0 AG" ~ "NY/11", 
                           No2 == "1 AG" ~ "IA/14", 
                           No2 == "7 AG" ~ "NE/15", 
                           No2 == "10 AG" ~ "MN/13",
                           No2 == "8 AG" ~ "OK/15",
                           No2 == "2 AG" ~ "NC/19",
                           No2 == "9 AG" ~ "MN/18",
                           No2 == "3 AG" ~ "WY/13",
                           No2 == "4 AG" ~ "IN/11",
                           #merged in antigens (marcus):
                           No2 == "52 AG" ~ "IA/13",
                           No2 == "53 AG" ~ "IA/A01290671/2013",
                           No2 == "60 AG" ~ "IL/A01475273/2014",
                           No2 == "61 AG" ~ "IL/A01503899/2014",
                           No2 == "65 AG" ~ "WY/A01397400/2013",
                           No2 == "97 AG" ~ "NC/16",
                           No2 == "100 AG" ~ "OK/16",
                           No2 == "104 AG" ~ "IL/15"
  )) %>%
  mutate(strain = factor(strain)) %>%
  mutate(id = case_when(No2 == "5 AG" ~ "1", 
                        No2 == "6 AG" ~ "2", 
                        No2 == "4 AG" ~ "3", 
                        No2 == "0 AG" ~ "4", 
                        No2 == "52 AG" ~ "5", 
                        No2 == "53 AG" ~ "6",
                        No2 == "10 AG" ~ "7",
                        No2 == "3 AG" ~ "8",
                        #No2 == "65 AG" ~ "9",
                        No2 == "1 AG" ~ "9",
                        #No2 == "60 AG" ~ "11",
                        #No2 == "61 AG" ~ "12",
                        #No2 == "104 AG" ~ "13",
                        No2 == "7 AG" ~ "10",
                        No2 == "8 AG" ~ "11",
                        No2 == "97 AG" ~ "12",
                        #No2 == "100 AG" ~ "17",
                        No2 == "9 AG" ~ "13",
                        No2 == "2 AG" ~ "14"
                        
  )) %>%
  mutate(id = factor(id)) %>%
  mutate(motif = case_when(No2 == "5 AG" ~ "NHKEYR", 
                           No2 == "6 AG" ~ "KHKEYS", 
                           No2 == "0 AG" | No2 == "8 AG" | No2 == "9 AG" | No2 == "52 AG" |
                             No2 == "53 AG" | No2 == "60 AG" | No2 == "61 AG" | No2 == "65 AG" |
                             No2 == "97 AG" | No2 == "100 AG" | No2 == "104 AG" ~ "NYNNYK", 
                           No2 == "1 AG" | No2 == "3 AG" ~ "KYNNYK", 
                           No2 == "7 AG" ~ "SYKNYK", 
                           No2 == "2 AG" | No2 == "10 AG" ~ "NYHNYK",
                           No2 == "4 AG" ~ "NYHGHE",

  )) %>%
  mutate(motif = factor(motif)) %>%
  mutate(clade = case_when(No2 == "5 AG" | No2 == "6 AG" ~ "HuVac", 
                           No2 == "0 AG" | No2 == "1 AG" | No2 == "7 AG" | No2 == "2 AG" |
                             No2 == "3 AG" | No2 == "8 AG" | No2 == "9 AG"|
                             No2 == "60 AG" | No2 == "61 AG" | No2 == "65 AG" | No2 == "97 AG" |
                             No2 == "100 AG" | No2 == "104 AG" ~ "CIV-A", 
                           No2 == "10 AG" | No2 == "52 AG" | No2 == "53 AG" ~ "CIV-B",
                           No2 == "4 AG" ~ "CIV-C"
  )) %>%
  mutate(clade = factor(clade)) %>%
  mutate(xlab = case_when(No1 == "8 AG" ~ "OK/15 \n (NYNNYK)",
                          No1 == "2 AG" ~ "NC/19 \n (NYHNYK)",
                          No1 == "9 AG" ~ "MN/18 \n (NYNNYK)"
                          
                            )) %>%
  mutate(xlab = factor(xlab))

data$motif <- factor(data$motif, levels = c("NYNNYK", "NYHNYK", "KYNNYK", "SYKNYK", "NHKEYR", "NYHGHE", "KHKEYS"))
data$xlab <- factor(data$xlab, levels = c("OK/15 \n (NYNNYK)", "NC/19 \n (NYHNYK)", "MN/18 \n (NYNNYK)"))
str(data)

g1 <- ggplot(data, aes(xlab, distance)) +
  geom_hline(yintercept = 3, linetype = "longdash") +
  geom_point(aes(color = motif, shape = clade), size = min_point_size) +
  #geom_label_repel(aes(label=strain, color=motif,segment.size = 1), size = min_point_size/1.5, fontface = "bold", show.legend = F, min.segment.length =  0, box.padding = 0.5) + 
  scale_shape_manual(values=c(16,15,17,7)) +
  scale_color_manual(values=c("#ED6141", "#f37b59", "#00b81f", "#00A5FF", "#C77Cff", "#C77Cff", "#FC61D5")) + 
  labs(x = "Antigen", y = "HA Antigenic Distance", shape = "H3 Clade", color = "Antigenic Motif") +
  theme (
    axis.text = element_text(size = min_text_size - 1, face = "bold"),
    axis.text.x = element_text(color = c("#ED6141", "#f37b59", "#ED6141")),
    axis.title = element_text(size = min_text_size, face = "bold"),
    legend.title = element_text(size = min_text_size, face= "bold"),
    legend.text = element_text(size = min_text_size)
  ) +
  guides(shape = guide_legend(override.aes = list(size=7)))

#### NI graph ###
data <- readxl:: read_excel("//iastate/lss/research/pcgauger-lab/Megan/civa/na_diversity/ACMACS_results/full_NI_ACMACS-all-distances.xlsx")

data <- data %>%
  rename("Antigen1" = "Name1") %>%
  rename("Antigen2" = "Name2") %>%
  mutate(Antigen1 = factor(Antigen1)) %>%
  filter(Antigen1 == "A(H1N2)/SWINE/OKLAHOMA/A01770190/2015 MDCK2" | 
           Antigen1 == "A(H1N2)/SWINE/NORTH CAROLINA/A02245294/2019 MDCK2" |
           Antigen1 == "A(H1N2)/SWINE/MINNESOTA/A02266068/2018 MDCK2"
  ) %>%
  mutate(Antigen2 = factor(Antigen2)) %>%
  filter(Antigen2 == "A(H1N2)/SWINE/NEBRASKA/A01492366/2014 MDCK3" |
           Antigen2 == "A(H1N2)/SWINE/NEW YORK/A01104005/2011 MDCK4" |
           Antigen2 == "A(H1N2)/SWINE/MINNESOTA/A02266068/2018 MDCK2" |
           Antigen2 == "A(H3N2)/SWINE/MINNESOTA/A01678475/2016 MDCK5" |
           Antigen2 == "A(H3N2)/SWINE/IOWA/A01480656/2014 MDCK4" |
           Antigen2 == "A(H1N2)/SWINE/OKLAHOMA/A01770190/2015 MDCK2" |
           Antigen2 == "A(H1N2)/SWINE/NORTH CAROLINA/A02245294/2019 MDCK2"
  ) %>%
  mutate(Antigen1 = factor(Antigen1)) %>%
  mutate(Antigen2 = factor(Antigen2)) %>%
  group_by(Antigen1)

data$strain <- data$No2
data$xlab <- data$No1
data$distance <- as.numeric(data$distance)

data <- data %>%
  mutate(clade = case_when(No2 == "2 AG" ~ "02A.1",
                           No2 == "0 AG" | No2 == "4 AG" ~ "02A.2",
                           No2 == "3 AG" ~ "02B.1",
                           No2 == "1 AG" | No2 == "5 AG" | No2 == "6 AG" ~ "02B.2",
                           
  )) %>%
  mutate(clade = factor(clade)) %>%
  mutate(xlab = case_when(No1 == "5 AG" ~ "OK/15\n(02B.2)",
                          No1 == "6 AG" ~ "NC/19\n(02B.2)",
                          No1 == "4 AG" ~ "MN/18\n(02A.2)"
                          
  )) %>%
  mutate(xlab = factor(xlab))

data$xlab <- factor(data$xlab, levels = c("OK/15\n(02B.2)", "NC/19\n(02B.2)", "MN/18\n(02A.2)"))

g2 <- ggplot(data, aes(xlab, distance)) +
  geom_hline(yintercept = 3, linetype = "longdash") +
  #geom_point(aes(color = clade), size = min_point_size, show.legend = F) +
  geom_point(aes(color = clade), size = min_point_size, show.legend = T) +
  scale_color_hue(l=40) +
  #geom_label_repel(min.segment.length =  0, size = min_point_size, aes(label = clade, color = clade, segment.size = 1), show.legend = F, box.padding = 0.5) +
  #geom_label(aes(label=clade,color=clade), size = 5, fontface = "bold", show.legend = F) +
  labs(x="Antigen", y="NA Antigenic Distance", color = "NA Lineage") +
  theme(
    axis.text = element_text(size = min_text_size - 1, face = "bold"),
    #axis.text.x = element_text(color = c("#C77CFF", "#C77CFF", "#7CAE00")),
    axis.text.x = element_text(color = darken(c("#C77CFF", "#C77CFF", "#7CAE00"),amount=0.4)),
    axis.title = element_text(size = min_text_size, face = "bold"),
    legend.title = element_text(size = min_text_size, face="bold"),
    legend.text = element_text(size=min_text_size)
  )

setwd("//iastate/lss/research/pcgauger-lab/Megan/IV_a Paper/c-iva/figures/")


g_combo <- g1 + g2 + plot_layout(width=c(1,1)) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = min_text_size, face="bold"))
g_combo <- g_combo & ylim(0,6.5)

tiff("Figure_4_HI_NI_antigenic_distance.tiff", units = "in", width = 11, height = 8.5, res = 300, compression = "lzw")
g_combo
dev.off()