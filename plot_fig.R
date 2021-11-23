
library(reshape2)
library(ggplot2)
library(scales)
library("viridis")

# requires running Protocol_1 first
pltd <- P1_results

txnamF <- rownames(pltd[[3]])

write.table(pltd[[7]], file = 'P1_heatmap.csv', sep = ",")
txordF <- sub("^(ub-).+", "14\\1BA",
sub("^(ub-)(AP).+", "15\\1\\2",
sub("^(di-)(RZ).+", "07\\1\\2",
sub("^(di-)(ST).+", "09\\1\\2",
sub("^(di-)(AL).+", "08\\1\\2",
sub("^(di-)(CR).+", "10\\1\\2",
sub("^(di-)(RP).+", "13\\1\\2",
sub("^(di-)(PL).+", "11\\1\\2",
sub("^(di-)(CH).+", "12\\1\\2",
sub("^(ex-)(EG).+", "06\\1\\2",
sub("^(ex-)(JA).+", "05\\1\\2",
sub("^(ex-)(HL).+", "04\\1\\2",
sub("^(am-)(AM).+", "03\\1\\2",
sub("^(am-)(FU).+", "02\\1\\2",
sub("^(am-)(HO).+", "01\\1\\2", txnamF)))))))))))))))


taxaMads <- rowMads(pltd[[7]])
ordF <- order(txordF, -taxaMads, decreasing = T)

mgm01 <- pltd[[3]][ordF,]  #missing data full
mgm02 <- pltd[[4]][ordF,]  #missing data full + top 10%
mgm03 <- pltd[[5]][ordF,]  #missing data full + top 10% + 50%
mgm04 <- pltd[[7]][ordF,]  #outlying scores #

txcolrF <- sub("^....-AM", "cornflowerblue",
sub("^....-FU", "deepskyblue",
sub("^....-HO", "blue",
sub("^....-AL", "olivedrab",
sub("^....-CH", "green3",
sub("^....-CR", "goldenrod",
sub("^....-PL", "green4",
sub("^....-RP", "chartreuse3",
sub("^....-RZ", "seagreen",
sub("^....-ST", "yellowgreen",
sub("^....-EG", "coral",
sub("^....-HL", "red",
sub("^....-JA", "magenta",
sub("^....-AP", "grey60",
sub("^....-BA", "black", txordF[ordF])))))))))))))))

rownames(mgm01) <-
  rownames(mgm02) <- 
  rownames(mgm03) <- 
  rownames(mgm04) <- #
  sub("^..-[^_]+_", "" , rownames(mgm01))

mgm04x <- mgm04 #

antrectx <- which(rowSums(mgm03)==0)
antrecty <- which(colSums(mgm03)==0)

mgm04x[,] <- 4 #
mgm04x[mgm02==0] <- 2 #
mgm04x[mgm01==0] <- 1 #

meltmgm04 <- melt(mgm04) #
meltmgm04x <- melt(mgm04x) #

## Data missing, black markers for removed genes and taxa
ggplot(data = meltmgm04x, aes(x=Var1, y=Var2, fill=factor(value))) +
labs(x = "Taxa", y = "Genes", fill = "Value") +
  geom_tile(color = "#999999") +
  scale_fill_manual(breaks = levels(factor(meltmgm04x$value)), values = c("white", "indianred1", "skyblue2"), labels =  c("Missing Data", "Top 10%" ,"Rest 90%")) +
   theme(text = element_text(size=7),
        axis.text.x = element_text(color = txcolrF, angle = 75, hjust = 1.1, vjust = 1.1, face = "bold"),
        legend.text=element_text(size=10),
        legend.title = element_blank(),
        legend.key = element_rect(color="black"),
        axis.title=element_text(size=14,face="bold")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  annotate("rect", xmin=antrectx+0.2, xmax=antrectx-.2, ymin=-0.2, ymax=.7, alpha=1, color='black', size=0.5, fill="black") +
  annotate("rect", ymin=antrecty+0.2, ymax=antrecty-.2, xmin=-0.2, xmax=.7, alpha=1, color='black', size=0.5, fill="black") 

# selected taxa for protocol 2
tx2 <- scan("intermediate_data/P2/P2_loci/taxaUsed_P2.txt", what = "character")
tx2 <- gsub('$','_', gsub('^\\w\\w-[^_]+_', '', tx2) )
tx2in <- which(rownames(mgm03) %in% tx2)

# Heat map, blue mark for selected taxa for P2
ggplot(data = meltmgm04, aes(x=Var1, y=Var2, fill=value)) +
  labs(x = "Taxa", y = "Genes", fill = "Value") +
  geom_tile(color = "black") +
  theme(text = element_text(size=7),
        axis.text.x = element_text(color = txcolrF, angle = 75, hjust = 1.1, vjust = 1.1, face = "bold"),
        legend.text=element_text(size=10),
        legend.title = element_blank(),
        legend.key = element_rect(color="black"),
        axis.title=element_text(size=14,face="bold")) +
  scale_fill_viridis(direction = -1, option = "B", na.value = "#FFFFFF", limits=c(quantile(meltmgm04$value, 0.20), quantile(meltmgm04$value, 0.98)), oob=squish) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  annotate("rect", xmin=tx2in-0.3, xmax=tx2in+.3, ymin=-0.5, ymax=.7, alpha=1, color='grey30', size=0.3, fill="dodgerblue") 

