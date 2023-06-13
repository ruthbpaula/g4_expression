## get_g4_at_gene_locations.sh
## Written by: Ruth De Paula and Albino Bacolla

# GENE EXPRESSION ANALYSIS - UHLEN:

# 1. Intersect coord of genomic region of interest with G4 coord
#### "get_g4_at_gene_locations.sh" - Generates lists of genes with G4 at exons, introns, 5'UTR, 3'UTR, TSS, TSS -45, any location, or no G4.

# 2. Convert genes to ensembl ids (https://biit.cs.ut.ee/gprofiler/convert) and grep with Uhlen table:
#### "preparation_for_uhlen_graph.sh":
# Converts genes to ensembl ids and greps with Uhlen table
# Runs script to make log2 of the mean of expression of each gene in Uhlen subtables
# Makes "uhlen_table_toplot_ok" with columns 3 and 4 from "uhlen_tableS18_mean" and column 3 switched with column 4

# 3. Run R script to make plots with the output from 3:

setwd("C:/Users/ruthb/Documents/UTH-MDA")

library(ggplot2)
library(ggpubr)
library("extrafont")

infile1 <- "uhlen_table_toplot_ok"
outfile1 <- "uhlen_table_toplot_ok"

myData <- read.table(file = infile1, sep="\t", header=FALSE)
colnames(myData) <- c("Group", "Gene_Expr")
head(myData)
xlabs <- paste(levels(myData$Group), "\n(n=", table(myData$Group),")", sep="")

compare_means(Gene_Expr ~ Group, data = myData)

ppi <- 300

png(file = paste(paste(outfile1),".png", sep=""), width = 5, height = 3, units = 'in', res = ppi)
plot <- ggplot(myData, aes(x=Group, y=Gene_Expr, fill = Group)) + 
  geom_boxplot(width = 0.45) +
  theme_classic() +
  #     geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.1) +
  #     ggtitle("Expression of genes that start at G4 DNA") +
  ylab("Mean expression log2(fpkm + 1)") +
  scale_x_discrete(labels = xlabs) +
  theme(plot.title = element_text(family = "Times New Roman", face = "bold.italic", color = "darkred", size = 8, hjust = 0.5),
        axis.title.x = element_text(family = "Helvetica", face = "plain", color = "blue", size = 8),
        axis.title.y = element_text(family = "Helvetica", face = "plain", color = "blue", size = 8),
        axis.text.x = element_text(family = "Helvetica", face = "italic", color = "black", size = 7, angle = 45),
        axis.text.y = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
        legend.text = element_text(family = "Helvetica", face = "bold", color = "darkgreen", size = 7),
        legend.title = element_text(family = "Helvetica", face = "italic", color = "blue", size = 7)) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, col = "black", fill = "black") #+ 
  #stat_compare_means(label.x = 1.1, label.y = 15.0, size = 2.7, color = "red")

plot

dev.off()


# Wilcoxon
list1 <- unique(myData$Group)

wilcox_matrix <- as.data.frame(matrix(nrow = length(list1), ncol = length(list1)))
colnames(wilcox_matrix) <- list1
rownames(wilcox_matrix) <- list1

for (i in 1:length(list1)) {
  region1 <- list1[i]
  for (j in 1:length(list1)) {
    region2 <- list1[j]
    set1 <- myData[myData$Group %in% region1,][,2]
    set2 <- myData[myData$Group %in% region2,][,2]
    wilcox_matrix[i,j] <- wilcox.test(set1, set2)$p.value
  }
}

write.table(wilcox_matrix, "uhlen_table_wilcox_matrix.txt")
