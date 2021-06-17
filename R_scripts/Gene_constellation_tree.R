library(ggtree)
library(ggplot2)
library(ape)
setwd("//iastate/lss/research/pcgauger-lab/Megan/IV_a Paper/c-iva/")

#read in genotype file
genotypes <- read.csv("nextstrain/metadata.tsv", sep="\t", stringsAsFactors = FALSE)
internal_genes <- c("strain", "PB2", "PB1", "PA", "NP", "MP", "NS")
genotypes <- genotypes[internal_genes]
colnames(genotypes) <- c("strain", "PB2", "PB1", "PA", "NP", "M", "NS")
#replace ? with NA
genotypes$PB2 <- gsub("?","",genotypes$PB2, fixed = TRUE)
genotypes$PB1 <- gsub("?","",genotypes$PB1, fixed = TRUE)
genotypes$PA <- gsub("?","",genotypes$PA, fixed = TRUE)
genotypes$NP <- gsub("?","",genotypes$NP, fixed = TRUE)
genotypes$M <- gsub("?","",genotypes$M, fixed = TRUE)
genotypes$NS <- gsub("?","",genotypes$NS, fixed = TRUE)
#reduce to only tips with full WGS
genotypes <- genotypes[!(is.na(genotypes$PB2) | genotypes$PB2==""),]
genotypes <- genotypes[!(is.na(genotypes$PB1) | genotypes$PB1==""),]
genotypes <- genotypes[!(is.na(genotypes$PA) | genotypes$PA==""),]
genotypes <- genotypes[!(is.na(genotypes$NP) | genotypes$NP==""),]
genotypes <- genotypes[!(is.na(genotypes$M) | genotypes$M==""),]
genotypes <- genotypes[!(is.na(genotypes$NS) | genotypes$NS==""),]

#set rownames to tips
rownames(genotypes) <- genotypes$strain
genotypes <- subset(genotypes, select=-(strain))
#change to correct label
genotypes[genotypes == "pdm"] <- "pdm09"

tree <- read.tree("nextstrain/nextstrain_selected-WGS-for-tree_timetree.nwk")
#tree <- drop.tip(tree, "A/swine/Ohio/12TOSU374/2012", trim.internal = FALSE)
circ <- ggtree(tree, mrsd = "2020-04-11", size = 0.8)+ theme_tree2()
circ
#df <- data.frame(matrix(NA, nrow=604,ncol=1))
#df$matrix.NA..nrow...604..ncol...1. <- "trig"
#colnames(df) <- "PB2"
#rownames(df) <- tree$tip.label

min_text_size <- 8

##need to edit tiff in GIMP to remove strain name residue 
tiff("figures/Figure_3_genotype_constellation_tree.tiff", units = "in", width = 7, height = 5, res = 300, compression = "lzw")

gheatmap(circ,genotypes,offset=0.1,width=0.4,colnames_angle = 45,colnames_offset_y = 0.25, colnames = F)+
  scale_x_ggtree() + 
  scale_fill_manual("Lineage", na.translate=F, values=c("yellow2","red3","green4")) +
 # scale_fill_discrete("Lineage", na.translate=F) +
  xlab("Time") +
  #geom_tiplab(align=T,size=4,linesize = .3) +
  theme(legend.position = c(.05, y=.945),
        legend.text = element_text(size=min_text_size, face = "bold"),
        legend.title = element_text(size=min_text_size, face="bold"),
        axis.text = element_text(size=min_text_size, face="bold"),
        axis.title = element_text(size=min_text_size+2, face="bold"),
        axis.line = element_line(size=1),
        axis.ticks = element_line(size=1))

dev.off()

##change tip labels to include whether we have WGS or not
full_tree <- read.tree("civa/test_tree.tre")
ggtree(full_tree) +geom_tiplab()
tips <- full_tree$tip.label
wgs_tips <- rownames(genotypes)
in_wgs <- ifelse(tips %in% wgs_tips, "|WGS", "|not")
df <- data.frame(tips, in_wgs)
df$new_header <- paste(df$tips,df$in_wgs, sep = "")
full_tree$tip.label <- df$new_header
ggtree(full_tree) + geom_tiplab()
write.tree(full_tree, file="civa/test_tree_WGS_appended.tree")

full_tree <- read.tree("civa/test_tree_WGS_appended.tre")
ggtree(full_tree) + geom_tiplab()
tips <- full_tree$tip.label
major <- read.table("civa/civa_major_clade.txt")
major_clade <- major$V1
minor <- read.table("civa/civa_minor_clade.txt")
minor_clade <- minor$V1
in_major_clade <- ifelse(tips %in% major_clade, "|major", "")
in_minor_clade <- ifelse(tips %in% minor_clade, "|minor", "")
df2 <- data.frame(tips,in_major_clade, in_minor_clade)
df2$new_header <- paste(df2$tips,df2$in_major_clade,df2$in_minor_clade, sep = "")
full_tree$tip.label <- df2$new_header
ggtree(full_tree) + geom_tiplab()
write.tree(full_tree, file="civa/test_tree_WGS_appended_major_minor.tre")


