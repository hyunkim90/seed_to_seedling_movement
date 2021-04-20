### community dissimilarity



### OTU tables
b.otu <- otu_table(bac.clean.log)
f.otu <- otu_table(fun.clean.log)

## Customized function
get_upper_tri <- function(cormat2){
  cormat2[lower.tri(cormat2)]<- NA
  return(cormat2)
}

### Fungal distance
bray.dist.fun<-vegdist(t(f.otu), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

df.bray.dist.fun<- as.matrix(bray.dist.fun)
rownames(df.bray.dist.fun)



upper_tri <- get_upper_tri(df.bray.dist.fun)
# Melt the correlation matrix
melted_tab.fun <- melt(upper_tri, na.rm = TRUE)
head(melted_tab.fun)
melted_tab.fun <- subset(melted_tab.fun, value > 0)

melted_tab.fun$kingdom <- "fungi"


melted_tab.fun.leaf.14 <- subset(melted_tab.fun,melted_tab.fun$Var1 %in% c("MSD14.L1", "MSD14.L2","MSD14.L3") & melted_tab.fun$Var2 %in% c("MSI14.L1", "MSI14.L2", "MSI14.L3"))
melted_tab.fun.stem.14 <- subset(melted_tab.fun,melted_tab.fun$Var1 %in% c("MSD14.S1", "MSD14.S2","MSD14.S3") & melted_tab.fun$Var2 %in% c("MSI14.S1", "MSI14.S2", "MSI14.S3"))
melted_tab.fun.root.14 <- subset(melted_tab.fun,melted_tab.fun$Var1 %in% c("MSD14.R1", "MSD14.R2","MSD14.R3") & melted_tab.fun$Var2 %in% c("MSI14.R1", "MSI14.R2", "MSI14.R3"))

melted_tab.fun.shoot.7 <- subset(melted_tab.fun,melted_tab.fun$Var1 %in% c("MSD7.Sh1", "MSD7.Sh2","MSD7.Sh3") & melted_tab.fun$Var2 %in% c("MSI7.Sh1", "MSI7.Sh2"))
melted_tab.fun.root.7 <- subset(melted_tab.fun,melted_tab.fun$Var1 %in% c("MSD7.R1", "MSD7.R2","MSD7.R3") & melted_tab.fun$Var2 %in% c("MSI7.R1", "MSI7.R2"))

melted_tab.fun.grain.4 <- subset(melted_tab.fun,melted_tab.fun$Var1 %in% c("MSD4.1", "MSD4.2","MSD4.3") & melted_tab.fun$Var2 %in% c("MSI4.1", "MSI4.2","MSI4.3"))
melted_tab.fun.grain.1 <- subset(melted_tab.fun,melted_tab.fun$Var1 %in% c("MSD1.1", "MSD1.2","MSD1.3") & melted_tab.fun$Var2 %in% c("MSI1.1", "MSI1.2","MSI1.3"))
melted_tab.fun.grain.0 <- subset(melted_tab.fun,melted_tab.fun$Var1 %in% c("MSD0.1", "MSD0.2","MSD0.3") & melted_tab.fun$Var2 %in% c("MSI0.1", "MSI0.2","MSI0.3"))

melted_tab.fun.leaf.14$Compartment <- "Leaf_14"
melted_tab.fun.stem.14$Compartment <- "Stem_14"
melted_tab.fun.root.14$Compartment <- "Root_14"
melted_tab.fun.root.7$Compartment <- "Root_7"
melted_tab.fun.shoot.7$Compartment <- "Shoot_7"
melted_tab.fun.grain.4$Compartment <- "Grain_4"
melted_tab.fun.grain.1$Compartment <- "Grain_1"
melted_tab.fun.grain.0$Compartment <- "Grain_0"

### Bacterial distance
bray.dist.bac<-vegdist(t(b.otu), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

df.bray.dist.bac<- as.matrix(bray.dist.bac)
rownames(df.bray.dist.bac)



upper_tri <- get_upper_tri(df.bray.dist.bac)
# Melt the correlation matrix
melted_tab.bac <- melt(upper_tri, na.rm = TRUE)
head(melted_tab.bac)
melted_tab.bac <- subset(melted_tab.bac, value > 0)

melted_tab.bac$kingdom <- "bacteria"



melted_tab.bac.leaf.14 <- subset(melted_tab.bac,melted_tab.bac$Var1 %in% c("H14L1.17", "H14L2.17","H14L3.17") & melted_tab.bac$Var2 %in% c("U14L1.17", "U14L2.17","U14L3.17"))
melted_tab.bac.stem.14 <- subset(melted_tab.bac,melted_tab.bac$Var1 %in% c("H14S1.17", "H14S2.17","H14S3.17") & melted_tab.bac$Var2 %in% c("U14S1.17", "U14S2.17","U14S3.17"))
melted_tab.bac.root.14 <- subset(melted_tab.bac,melted_tab.bac$Var1 %in% c("H14R1.17", "H14R2.17","H14R3.17") & melted_tab.bac$Var2 %in% c("U14R1.17", "U14R2.17","U14R3.17"))

melted_tab.bac.shoot.7 <- subset(melted_tab.bac,melted_tab.bac$Var1 %in% c("H7Sh1.17", "H7Sh2.17","H7Sh3.17") & melted_tab.bac$Var2 %in% c("U7Sh1.17", "U7Sh2.17"))
melted_tab.bac.root.7 <- subset(melted_tab.bac,melted_tab.bac$Var1 %in% c("H7R1.17", "H7R2.17","H7R3.17") & melted_tab.bac$Var2 %in% c("U7R1.17", "U7R2.17"))

melted_tab.bac.grain.4 <- subset(melted_tab.bac,melted_tab.bac$Var1 %in% c("H4G1.17", "H4G2.17","H4G3.17") & melted_tab.bac$Var2 %in% c("U4G1.17", "U4G2.17","U4G3.17"))
melted_tab.bac.grain.1 <- subset(melted_tab.bac,melted_tab.bac$Var1 %in% c("H1G1.17", "H1G2.17","H1G3.17") & melted_tab.bac$Var2 %in% c("U1G1.17", "U1G2.17","U1G3.17"))
melted_tab.bac.grain.0 <- subset(melted_tab.bac,melted_tab.bac$Var1 %in% c("H0G1.17", "H0G2.17","H0G3.17") & melted_tab.bac$Var2 %in% c("U0G1.17", "U0G2.17","U0G3.17"))


melted_tab.bac.leaf.14$Compartment <- "Leaf_14"
melted_tab.bac.stem.14$Compartment <- "Stem_14"
melted_tab.bac.root.14$Compartment <- "Root_14"
melted_tab.bac.root.7$Compartment <- "Root_7"
melted_tab.bac.shoot.7$Compartment <- "Shoot_7"
melted_tab.bac.grain.4$Compartment <- "Grain_4"
melted_tab.bac.grain.1$Compartment <- "Grain_1"
melted_tab.bac.grain.0$Compartment <- "Grain_0"

### Merging tables
melted_tab.merged<-rbind(melted_tab.bac.grain.0, melted_tab.bac.grain.1,melted_tab.bac.grain.4,melted_tab.bac.shoot.7,melted_tab.bac.root.7,
      melted_tab.bac.leaf.14,melted_tab.bac.stem.14,melted_tab.bac.root.14,
      melted_tab.fun.grain.0, melted_tab.fun.grain.1,melted_tab.fun.grain.4,melted_tab.fun.shoot.7,melted_tab.fun.root.7,
      melted_tab.fun.leaf.14,melted_tab.fun.stem.14,melted_tab.fun.root.14)
      


x <- melted_tab.bac.grain.0$value
y <- melted_tab.fun.grain.0$value
wilcox.test(x, y, conf.int = TRUE)  #0.2973

x <- melted_tab.bac.grain.1$value
y <- melted_tab.fun.grain.1$value
wilcox.test(x, y, conf.int = TRUE)  #0.0004936

x <- melted_tab.bac.grain.4$value
y <- melted_tab.fun.grain.4$value
wilcox.test(x, y, conf.int = TRUE)  #4.114e-05

x <- melted_tab.bac.shoot.7$value
y <- melted_tab.fun.shoot.7$value
wilcox.test(x, y, conf.int = TRUE)  #0.2403

x <- melted_tab.bac.root.7$value
y <- melted_tab.fun.root.7$value
wilcox.test(x, y, conf.int = TRUE)  #0.002165


#### plotting
library(ggpubr)

p <- ggboxplot(data = melted_tab.merged, x="Compartment", y="value", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("bacteria" = "#6699CC", "fungi" = "#CC9900"))+
  ylab("Bray-Curtis dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = kingdom), label = "p.signif", method = "wilcox.test")+
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed")

p
