### Niche differentiation by Jaccard dissimilarity
### OTU tables
b.otu <- otu_table(bac.clean.log)
f.otu <- otu_table(fun.clean.log)

## Customized function
get_upper_tri <- function(cormat2){
  cormat2[lower.tri(cormat2)]<- NA
  return(cormat2)
}

### Fungal distance
jaccard.dist.fun<-vegdist(t(f.otu), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

df.jaccard.dist.fun<- as.matrix(jaccard.dist.fun)
rownames(df.jaccard.dist.fun)



upper_tri <- get_upper_tri(df.jaccard.dist.fun)
# Melt the correlation matrix
melted_tab.fun <- melt(upper_tri, na.rm = TRUE)
head(melted_tab.fun)
melted_tab.fun <- subset(melted_tab.fun, value > 0)

melted_tab.fun$kingdom <- "fungi"

### Comparison with seeds at the sowing stage
melted_tab.fun.leaf.14 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSD14.L1", "MSD14.L2","MSD14.L3") & melted_tab.fun$Var1 %in% c("MSD0.1", "MSD0.2","MSD0.3"))
melted_tab.fun.stem.14 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSD14.S1", "MSD14.S2","MSD14.S3") & melted_tab.fun$Var1 %in% c("MSD0.1", "MSD0.2","MSD0.3"))
melted_tab.fun.root.14 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSD14.R1", "MSD14.R2","MSD14.R3") & melted_tab.fun$Var1 %in% c("MSD0.1", "MSD0.2","MSD0.3"))

melted_tab.fun.shoot.7 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSD7.Sh1", "MSD7.Sh2","MSD7.Sh3") & melted_tab.fun$Var1 %in% c("MSD0.1", "MSD0.2","MSD0.3"))
melted_tab.fun.root.7 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSD7.R1", "MSD7.R2","MSD7.R3") & melted_tab.fun$Var1 %in% c("MSD0.1", "MSD0.2","MSD0.3"))

melted_tab.fun.grain.4 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSD4.1", "MSD4.2","MSD4.3") & melted_tab.fun$Var1 %in% c("MSD0.1", "MSD0.2","MSD0.3"))
melted_tab.fun.grain.1 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSD1.1", "MSD1.2","MSD1.3") & melted_tab.fun$Var1 %in% c("MSD0.1", "MSD0.2","MSD0.3"))



melted_tab.fun.leaf.14$Compartment <- "Leaf_14_Hulled"
melted_tab.fun.stem.14$Compartment <- "Stem_14_Hulled"
melted_tab.fun.root.14$Compartment <- "Root_14_Hulled"
melted_tab.fun.root.7$Compartment <- "Root_7_Hulled"
melted_tab.fun.shoot.7$Compartment <- "Shoot_7_Hulled"
melted_tab.fun.grain.4$Compartment <- "Grain_4_Hulled"
melted_tab.fun.grain.1$Compartment <- "Grain_1_Hulled"



##Intact seeds
intact_melted_tab.fun.leaf.14 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSI14.L1", "MSI14.L2","MSI14.L3") & melted_tab.fun$Var1 %in% c("MSI0.1", "MSI0.2","MSI0.3"))
intact_melted_tab.fun.stem.14 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSI14.S1", "MSI14.S2","MSI14.S3") & melted_tab.fun$Var1 %in% c("MSI0.1", "MSI0.2","MSI0.3"))
intact_melted_tab.fun.root.14 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSI14.R1", "MSI14.R2","MSI14.R3") & melted_tab.fun$Var1 %in% c("MSI0.1", "MSI0.2","MSI0.3"))

intact_melted_tab.fun.shoot.7 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSI7.Sh1", "MSI7.Sh2","MSI7.Sh3") & melted_tab.fun$Var1 %in% c("MSI0.1", "MSI0.2","MSI0.3"))
intact_melted_tab.fun.root.7 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSI7.R1", "MSI7.R2","MSI7.R3") & melted_tab.fun$Var1 %in% c("MSI0.1", "MSI0.2","MSI0.3"))

intact_melted_tab.fun.grain.4 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSI4.1", "MSI4.2","MSI4.3") & melted_tab.fun$Var1 %in% c("MSI0.1", "MSI0.2","MSI0.3"))
intact_melted_tab.fun.grain.1 <- subset(melted_tab.fun,melted_tab.fun$Var2 %in% c("MSI1.1", "MSI1.2","MSI1.3") & melted_tab.fun$Var1 %in% c("MSI0.1", "MSI0.2","MSI0.3"))



intact_melted_tab.fun.leaf.14$Compartment <- "Leaf_14_Intact"
intact_melted_tab.fun.stem.14$Compartment <- "Stem_14_Intact"
intact_melted_tab.fun.root.14$Compartment <- "Root_14_Intact"
intact_melted_tab.fun.root.7$Compartment <- "Root_7_Intact"
intact_melted_tab.fun.shoot.7$Compartment <- "Shoot_7_Intact"
intact_melted_tab.fun.grain.4$Compartment <- "Grain_4_Intact"
intact_melted_tab.fun.grain.1$Compartment <- "Grain_1_Intact"


### Bacterial distance
jaccard.dist.bac<-vegdist(t(b.otu), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

df.jaccard.dist.bac<- as.matrix(jaccard.dist.bac)
rownames(df.jaccard.dist.bac)



upper_tri <- get_upper_tri(df.jaccard.dist.bac)
# Melt the correlation matrix
melted_tab.bac <- melt(upper_tri, na.rm = TRUE)
head(melted_tab.bac)
melted_tab.bac <- subset(melted_tab.bac, value > 0)

melted_tab.bac$kingdom <- "bacteria"



### Comparison with seeds at the sowing stage
melted_tab.bac.leaf.14 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("H14L1.17", "H14L2.17","H14L3.17") & melted_tab.bac$Var1 %in% c("H0G1.17", "H0G2.17","H0G3.17"))
melted_tab.bac.stem.14 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("H14S1.17", "H14S2.17","H14S3.17") & melted_tab.bac$Var1 %in% c("H0G1.17", "H0G2.17","H0G3.17"))
melted_tab.bac.root.14 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("H14R1.17", "H14R2.17","H14R3.17") & melted_tab.bac$Var1 %in% c("H0G1.17", "H0G2.17","H0G3.17"))

melted_tab.bac.shoot.7 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("H7Sh1.17", "H7Sh2.17","H7Sh3.17") & melted_tab.bac$Var1 %in% c("H0G1.17", "H0G2.17","H0G3.17"))
melted_tab.bac.root.7 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("H7R1.17", "H7R2.17","H7R3.17") & melted_tab.bac$Var1 %in% c("H0G1.17", "H0G2.17","H0G3.17"))

melted_tab.bac.grain.4 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("H4G1.17", "H4G2.17","H4G3.17") & melted_tab.bac$Var1 %in% c("H0G1.17", "H0G2.17","H0G3.17"))
melted_tab.bac.grain.1 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("H1G1.17", "H1G2.17","H1G3.17") & melted_tab.bac$Var1 %in% c("H0G1.17", "H0G2.17","H0G3.17"))



melted_tab.bac.leaf.14$Compartment <- "Leaf_14_Hulled"
melted_tab.bac.stem.14$Compartment <- "Stem_14_Hulled"
melted_tab.bac.root.14$Compartment <- "Root_14_Hulled"
melted_tab.bac.root.7$Compartment <- "Root_7_Hulled"
melted_tab.bac.shoot.7$Compartment <- "Shoot_7_Hulled"
melted_tab.bac.grain.4$Compartment <- "Grain_4_Hulled"
melted_tab.bac.grain.1$Compartment <- "Grain_1_Hulled"



##Intact seeds
intact_melted_tab.bac.leaf.14 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("U14L1.17", "U14L2.17","U14L3.17") & melted_tab.bac$Var1 %in% c("U0G1.17", "U0G2.17","U0G3.17"))
intact_melted_tab.bac.stem.14 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("U14S1.17", "U14S2.17","U14S3.17") & melted_tab.bac$Var1 %in% c("U0G1.17", "U0G2.17","U0G3.17"))
intact_melted_tab.bac.root.14 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("U14R1.17", "U14R2.17","U14R3.17") & melted_tab.bac$Var1 %in% c("U0G1.17", "U0G2.17","U0G3.17"))

intact_melted_tab.bac.shoot.7 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("U7Sh1.17", "U7Sh2.17","U7Sh3.17") & melted_tab.bac$Var1 %in% c("U0G1.17", "U0G2.17","U0G3.17"))
intact_melted_tab.bac.root.7 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("U7R1.17", "U7R2.17","U7R3.17") & melted_tab.bac$Var1 %in% c("U0G1.17", "U0G2.17","U0G3.17"))

intact_melted_tab.bac.grain.4 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("U4G1.17", "U4G2.17","U4G3.17") & melted_tab.bac$Var1 %in% c("U0G1.17", "U0G2.17","U0G3.17"))
intact_melted_tab.bac.grain.1 <- subset(melted_tab.bac,melted_tab.bac$Var2 %in% c("U1G1.17", "U1G2.17","U1G3.17") & melted_tab.bac$Var1 %in% c("U0G1.17", "U0G2.17","U0G3.17"))


intact_melted_tab.bac.leaf.14$Compartment <- "Leaf_14_Intact"
intact_melted_tab.bac.stem.14$Compartment <- "Stem_14_Intact"
intact_melted_tab.bac.root.14$Compartment <- "Root_14_Intact"
intact_melted_tab.bac.root.7$Compartment <- "Root_7_Intact"
intact_melted_tab.bac.shoot.7$Compartment <- "Shoot_7_Intact"
intact_melted_tab.bac.grain.4$Compartment <- "Grain_4_Intact"
intact_melted_tab.bac.grain.1$Compartment <- "Grain_1_Intact"


melted_tab.merged_hulled<-rbind(melted_tab.bac.grain.1,melted_tab.bac.grain.4,melted_tab.bac.shoot.7,melted_tab.bac.root.7,
                         melted_tab.bac.leaf.14,melted_tab.bac.stem.14,melted_tab.bac.root.14,
                         melted_tab.fun.grain.1,melted_tab.fun.grain.4,melted_tab.fun.shoot.7,melted_tab.fun.root.7,
                         melted_tab.fun.leaf.14,melted_tab.fun.stem.14,melted_tab.fun.root.14)


intact_melted_tab.merged<-rbind(intact_melted_tab.bac.grain.1,intact_melted_tab.bac.grain.4,intact_melted_tab.bac.shoot.7,intact_melted_tab.bac.root.7,
                                intact_melted_tab.bac.leaf.14,intact_melted_tab.bac.stem.14,intact_melted_tab.bac.root.14,
                                intact_melted_tab.fun.grain.1,intact_melted_tab.fun.grain.4,intact_melted_tab.fun.shoot.7,intact_melted_tab.fun.root.7,
                                intact_melted_tab.fun.leaf.14,intact_melted_tab.fun.stem.14,intact_melted_tab.fun.root.14)


#### plotting
library(ggpubr)

melted_tab.merged_hulled.bac <- subset(melted_tab.merged_hulled, kingdom == "bacteria")
melted_tab.merged_hulled.fun <- subset(melted_tab.merged_hulled, kingdom == "fungi")
###Kruskal-Wallis test
max.value <- aggregate(melted_tab.merged_hulled.bac$value, by = list(melted_tab.merged_hulled.bac$Compartment), max)
colnames(max.value) <- c("Group", "MaxValue")

##Kruskal-Wallis test

melted_tab.merged_hulled.bac$Compartment <- as.factor(melted_tab.merged_hulled.bac$Compartment)

kw<-kruskal.test(value ~ Compartment, data = melted_tab.merged_hulled.bac)
kw$p.value

#library(FSA)
DT = dunnTest(value ~ Compartment, data = melted_tab.merged_hulled.bac,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.value, by = 'Group')


write.csv(PT, "Jaccard distance_statistical analyses_bacteria_hulled.csv")

p <- ggboxplot(data = subset(melted_tab.merged_hulled, kingdom == "bacteria"), x="Compartment", y="value", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("bacteria" = "#6699CC", "fungi" = "#CC9900"))+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+geom_point(size = 2.5, alpha = 0.5, position = "jitter", shape = 1)+
  geom_text(data=hsd1,aes(x=Group,y=MaxValue, label=Letter), vjust=-1)
p


max.value <- aggregate(melted_tab.merged_hulled.bac$value, by = list(melted_tab.merged_hulled.bac$Compartment), max)
colnames(max.value) <- c("Group", "MaxValue")



##Kruskal-Wallis test
max.value <- aggregate(melted_tab.merged_hulled.fun$value, by = list(melted_tab.merged_hulled.fun$Compartment), max)
colnames(max.value) <- c("Group", "MaxValue")

melted_tab.merged_hulled.fun$Compartment <- as.factor(melted_tab.merged_hulled.fun$Compartment)

kw<-kruskal.test(value ~ Compartment, data = melted_tab.merged_hulled.fun)
kw$p.value

#library(FSA)
DT = dunnTest(value ~ Compartment, data = melted_tab.merged_hulled.fun,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.value, by = 'Group')

write.csv(PT, "Jaccard distance_statistical analyses_fungi_hulled.csv")

p <- ggboxplot(data = subset(melted_tab.merged_hulled.fun, kingdom == "fungi"), x="Compartment", y="value", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("bacteria" = "#6699CC", "fungi" = "#CC9900"))+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+geom_point(size = 2.5, alpha = 0.5, position = "jitter", shape = 1)+
  geom_text(data=hsd1,aes(x=Group,y=MaxValue, label=Letter), vjust=-1)
p

##Kruskal-Wallis test
intact_melted_tab.merged.bac <- subset(intact_melted_tab.merged, kingdom == "bacteria")
intact_melted_tab.merged.fun <- subset(intact_melted_tab.merged, kingdom == "fungi")


max.value <- aggregate(intact_melted_tab.merged.bac$value, by = list(intact_melted_tab.merged.bac$Compartment), max)
colnames(max.value) <- c("Group", "MaxValue")

intact_melted_tab.merged.bac$Compartment <- as.factor(intact_melted_tab.merged.bac$Compartment)

kw<-kruskal.test(value ~ Compartment, data = intact_melted_tab.merged.bac)
kw$p.value

#library(FSA)
DT = dunnTest(value ~ Compartment, data = intact_melted_tab.merged.bac,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.value, by = 'Group')

write.csv(PT, "Jaccard distance_statistical analyses_bacteria_intact.csv")


p <- ggboxplot(data = subset(intact_melted_tab.merged.bac, kingdom == "bacteria"), x="Compartment", y="value", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("bacteria" = "#6699CC", "fungi" = "#CC9900"))+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+geom_point(size = 2.5, alpha = 0.5, position = "jitter", shape = 1)+
  geom_text(data=hsd1,aes(x=Group,y=MaxValue, label=Letter), vjust=-1)
p



max.value <- aggregate(intact_melted_tab.merged.fun$value, by = list(intact_melted_tab.merged.fun$Compartment), max)
colnames(max.value) <- c("Group", "MaxValue")

intact_melted_tab.merged.fun$Compartment <- as.factor(intact_melted_tab.merged.fun$Compartment)

kw<-kruskal.test(value ~ Compartment, data = intact_melted_tab.merged.fun)
kw$p.value

#library(FSA)
DT = dunnTest(value ~ Compartment, data = intact_melted_tab.merged.fun,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.value, by = 'Group')

write.csv(PT, "Jaccard distance_statistical analyses_fungi_intact.csv")


p <- ggboxplot(data = subset(intact_melted_tab.merged.fun, kingdom == "fungi"), x="Compartment", y="value", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("bacteria" = "#6699CC", "fungi" = "#CC9900"))+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+geom_point(size = 2.5, alpha = 0.5, position = "jitter", shape = 1)+
  geom_text(data=hsd1,aes(x=Group,y=MaxValue, label=Letter), vjust=-1)
p


dev.off()

write.csv(intact_melted_tab.merged, "Jaccard distance with seed 0_intact.csv")
write.csv(melted_tab.merged_hulled, "Jaccard distance with seed 0_hulled.csv")
