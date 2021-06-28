library(ggtree)
library(ggplot2)
library(ape)
library(tidytree)
library(lubridate)

setwd("//iastate/lss/research/pcgauger-lab/Megan/IV_a Paper/c-iva/")

#read in tree
tree <- read.tree("beast_analysis/c-iva_sequences_alignment_RAxML_tree.tre")

ggtree(tree) + geom_text(aes(label=node))

tree <- root(tree,466)
nodelabels(cex = .75, bg = "yellow")

#split tip label into desired columns
tips <- tree$tip.label
df <- as.data.frame(strsplit(tips, "[|]"))
df<- as.data.frame(t(df))
df$label <- tips
colnames(df) <- c("strain","date","state","label")

#create Date objects
df$format_date <- as.Date(df$date, "%m/%d/%Y")
#deal with month-year only dates
for(i in 1:length(df$format_date)) {
  if(is.na(df$format_date[i])) {
    df$format_date[i] <- as.Date(parse_date_time(df$date[i], orders = "my"))
  }
}
#deal with year only dates
for(i in 1:length(df$format_date)) {
  if(is.na(df$format_date[i])) {
    df$format_date[i] <- as.Date(df$date[i], "%Y")
  }
}

#add column with year information
df$year <- year(df$format_date)

#join data to phylogenetic tree stored as tibble
tib <- as_tibble(tree)
y <- full_join(tib, df, by= 'label')

#convert to treedata
z <- as.treedata(y)


tiff("supplementary files/Supplementary_Figure_1_ML_Tree.tiff", units = "in", width = 7, height = 5, res = 300, compression = "lzw")

#graph tree
ggtree(z,aes(color=year)) +
  scale_color_continuous(low="blue", high="red") +
  labs(color="Time")
  
dev.off()



