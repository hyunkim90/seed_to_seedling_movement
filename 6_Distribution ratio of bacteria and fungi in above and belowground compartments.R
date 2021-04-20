### Comparison of aboveground and belowground parts in terms of abundance and OTU richness
bac.clean.ss.U7R #19
bac.clean.ss.U7Sh #17

bac.clean.ss.U14R #22
bac.clean.ss.U14S #33
bac.clean.ss.U14L #28

bac.clean.ss.H7R #40
bac.clean.ss.H7Sh #20

bac.clean.ss.H14R #122
bac.clean.ss.H14S #17
bac.clean.ss.H14L #15


bac.clean.ss
fun.clean.ss

### abundance

bac.melt <- bac.clean.ss %>% psmelt() %>% group_by(Group)%>% summarise(SumAbund = sum(Abundance))
bac.melt.total <- bac.clean.ss %>% psmelt() %>% group_by(Hull, Age)%>% summarise(TotalAbund = sum(Abundance))
head(bac.melt)

bac.melt$TotalAbund <- 0
bac.melt$TotalAbund[bac.melt$Group == "H0"] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 0]
bac.melt$TotalAbund[bac.melt$Group == "H1"] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 1]
bac.melt$TotalAbund[bac.melt$Group == "H4"] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 4]
bac.melt$TotalAbund[bac.melt$Group %in% c("H7Sh", "H7R")] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 7]
bac.melt$TotalAbund[bac.melt$Group %in% c("H14L", "H14S","H14R")] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 14]

bac.melt$TotalAbund[bac.melt$Group == "U0"] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "intact" & bac.melt.total$Age == 0]
bac.melt$TotalAbund[bac.melt$Group == "U1"] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "intact" & bac.melt.total$Age == 1]
bac.melt$TotalAbund[bac.melt$Group == "U4"] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "intact" & bac.melt.total$Age == 4]
bac.melt$TotalAbund[bac.melt$Group %in% c("U7Sh", "U7R")] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "intact" & bac.melt.total$Age == 7]
bac.melt$TotalAbund[bac.melt$Group %in% c("U14L", "U14S","U14R")] <- bac.melt.total$TotalAbund[bac.melt.total$Hull == "intact" & bac.melt.total$Age == 14]

write.xlsx(bac.melt, "bacterial abundance in each sample.xlsx")
### Abundance ratio between aboveground and belowground parts
## hulled grain 7 
ratio.7.R <- bac.melt$SumAbund[which(bac.melt$Group == "H7R")]/bac.melt.total$TotalAbund[which(bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 7)]
ratio.7 <- bac.melt$SumAbund[which(bac.melt$Group == "H7Sh")]/bac.melt.total$TotalAbund[which(bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 7)]

ratio.14.R <- bac.melt$SumAbund[which(bac.melt$Group == "H14R")]/bac.melt.total$TotalAbund[which(bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 14)]
ratio.14 <- (bac.melt$SumAbund[which(bac.melt$Group == "H14S")]+bac.melt$SumAbund[which(bac.melt$Group == "H14L")])/bac.melt.total$TotalAbund[which(bac.melt.total$Hull == "hulled" & bac.melt.total$Age == 14)]


### unhulled grains
ratio.7.R.U <- bac.melt$SumAbund[which(bac.melt$Group == "U7R")]/bac.melt.total$TotalAbund[which(bac.melt.total$Hull == "intact" & bac.melt.total$Age == 7)]
ratio.7.U <- bac.melt$SumAbund[which(bac.melt$Group == "U7Sh")]/bac.melt.total$TotalAbund[which(bac.melt.total$Hull == "intact" & bac.melt.total$Age == 7)]

ratio.14.R.U <- bac.melt$SumAbund[which(bac.melt$Group == "U14R")]/bac.melt.total$TotalAbund[which(bac.melt.total$Hull == "intact" & bac.melt.total$Age == 14)]
ratio.14.U <- (bac.melt$SumAbund[which(bac.melt$Group == "U14S")]+bac.melt$SumAbund[which(bac.melt$Group == "U14L")])/bac.melt.total$TotalAbund[which(bac.melt.total$Hull == "intact" & bac.melt.total$Age == 14)]


### Diversity - Shannon and Chao1 or ACE
# bac.clean.ss.group <- merge_phyloseq(bac.clean.ss, "Group")
tab.bac.alpha<-microbiome::alpha(bac.clean.ss, index = "all")
b.meta <-b.map
b.meta$Richness <- tab.bac.alpha$chao1
b.meta$Shannon <- tab.bac.alpha$diversity_shannon


## Richness of below and aboveground compartments
## hulled
chao1.below.7 <- mean(b.meta$Richness[which(b.meta$Group == "H7R")])
chao1.above.7 <- mean(b.meta$Richness[which(b.meta$Group =="H7Sh")])

chao1.below.14 <- mean(b.meta$Richness[which(b.meta$Group == "H14R")])
chao1.above.14 <- mean(b.meta$Richness[which(b.meta$Group =="H14L" | b.meta$Group =="H14S")])

## unhulled
chao1.below.7 <- mean(b.meta$Richness[which(b.meta$Group == "U7R")])
chao1.above.7 <- mean(b.meta$Richness[which(b.meta$Group =="U7Sh")])

chao1.below.14 <- mean(b.meta$Richness[which(b.meta$Group == "U14R")])
chao1.above.14 <- mean(b.meta$Richness[which(b.meta$Group =="U14L" | b.meta$Group =="U14S")])


## Shannon of below and aboveground compartments
## hulled
Shannon.below.7 <- mean(b.meta$Shannon[which(b.meta$Group == "H7R")])
Shannon.above.7 <- mean(b.meta$Shannon[which(b.meta$Group =="H7Sh")])

Shannon.below.14 <- mean(b.meta$Shannon[which(b.meta$Group == "H14R")])
Shannon.above.14 <- mean(b.meta$Shannon[which(b.meta$Group =="H14L" | b.meta$Group =="H14S")])

## unhulled
Shannon.below.7 <- mean(b.meta$Shannon[which(b.meta$Group == "U7R")])
Shannon.above.7 <- mean(b.meta$Shannon[which(b.meta$Group =="U7Sh")])

Shannon.below.14 <- mean(b.meta$Shannon[which(b.meta$Group == "U14R")])
Shannon.above.14 <- mean(b.meta$Shannon[which(b.meta$Group =="U14L" | b.meta$Group =="U14S")])



### Rarefied data
bac.rarefy <- rarefy_even_depth(bac.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
fun.rarefy <- rarefy_even_depth(fun.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

### Diversity - Shannon and Chao1 or ACE
# bac.clean.ss.group <- merge_phyloseq(bac.clean.ss, "Group")
tab.bac.alpha.rarefy<-microbiome::alpha(bac.rarefy, index = "all")
b.meta$Richness2 <- tab.bac.alpha.rarefy$chao1
b.meta$Shannon2 <- tab.bac.alpha.rarefy$diversity_shannon


## Richness of below and aboveground compartments
## hulled
chao1.below.7 <- mean(b.meta$Richness2[which(b.meta$Group == "H7R")])
chao1.above.7 <- mean(b.meta$Richness2[which(b.meta$Group =="H7Sh")])

chao1.below.14 <- mean(b.meta$Richness2[which(b.meta$Group == "H14R")])
chao1.above.14 <- mean(b.meta$Richness2[which(b.meta$Group =="H14L" | b.meta$Group =="H14S")])

## unhulled
chao1.below.7 <- mean(b.meta$Richness2[which(b.meta$Group == "U7R")])
chao1.above.7 <- mean(b.meta$Richness2[which(b.meta$Group =="U7Sh")])

chao1.below.14 <- mean(b.meta$Richness2[which(b.meta$Group == "U14R")])
chao1.above.14 <- mean(b.meta$Richness2[which(b.meta$Group =="U14L" | b.meta$Group =="U14S")])


## Shannon of below and aboveground compartments
## hulled
Shannon.below.7 <- mean(b.meta$Shannon2[which(b.meta$Group == "H7R")])
Shannon.above.7 <- mean(b.meta$Shannon2[which(b.meta$Group =="H7Sh")])

Shannon.below.14 <- mean(b.meta$Shannon2[which(b.meta$Group == "H14R")])
Shannon.above.14 <- mean(b.meta$Shannon2[which(b.meta$Group =="H14L" | b.meta$Group =="H14S")])

## unhulled
Shannon.below.7 <- mean(b.meta$Shannon2[which(b.meta$Group == "U7R")])
Shannon.above.7 <- mean(b.meta$Shannon2[which(b.meta$Group =="U7Sh")])

Shannon.below.14 <- mean(b.meta$Shannon2[which(b.meta$Group == "U14R")])
Shannon.above.14 <- mean(b.meta$Shannon2[which(b.meta$Group =="U14L" | b.meta$Group =="U14S")])

# wilcoxon test
x <- b.meta$Shannon2[which(b.meta$Group == "U14R")]
y <- b.meta$Shannon2[which(b.meta$Group =="U14L" | b.meta$Group =="U14S")]
wilcox.test(x, y, conf.int = TRUE) #0.04762

x <- b.meta$Shannon2[which(b.meta$Group == "H14R")]
y <- b.meta$Shannon2[which(b.meta$Group =="H14L" | b.meta$Group =="H14S")]
wilcox.test(x, y, conf.int = TRUE) #0.02381

x <- b.meta$Shannon2[which(b.meta$Group == "U7R")]
y <- b.meta$Shannon2[which(b.meta$Group =="U7Sh")]
wilcox.test(x, y, conf.int = TRUE) #0.3333

x <- b.meta$Shannon2[which(b.meta$Group == "H7R")]
y <- b.meta$Shannon2[which(b.meta$Group =="H7Sh")]
wilcox.test(x, y, conf.int = TRUE) #0.1


### Fungi
sample_data(fun.clean.ss) <- sample_data(f.map)
fun.melt <- fun.clean.ss %>% psmelt() %>% group_by(Group)%>% summarise(SumAbund = sum(Abundance))
fun.melt.total <- fun.clean.ss %>% psmelt() %>% group_by(Hull, Age)%>% summarise(TotalAbund = sum(Abundance))
head(fun.melt)

fun.melt$TotalAbund <- 0
fun.melt$TotalAbund[fun.melt$Group == "H0"] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 0]
fun.melt$TotalAbund[fun.melt$Group == "H1"] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 1]
fun.melt$TotalAbund[fun.melt$Group == "H4"] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 4]
fun.melt$TotalAbund[fun.melt$Group %in% c("H7Sh", "H7R")] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 7]
fun.melt$TotalAbund[fun.melt$Group %in% c("H14L", "H14S","H14R")] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 14]

fun.melt$TotalAbund[fun.melt$Group == "U0"] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "intact" & fun.melt.total$Age == 0]
fun.melt$TotalAbund[fun.melt$Group == "U1"] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "intact" & fun.melt.total$Age == 1]
fun.melt$TotalAbund[fun.melt$Group == "U4"] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "intact" & fun.melt.total$Age == 4]
fun.melt$TotalAbund[fun.melt$Group %in% c("U7Sh", "U7R")] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "intact" & fun.melt.total$Age == 7]
fun.melt$TotalAbund[fun.melt$Group %in% c("U14L", "U14S","U14R")] <- fun.melt.total$TotalAbund[fun.melt.total$Hull == "intact" & fun.melt.total$Age == 14]

write.xlsx(fun.melt, "fungal abundance in each sample.xlsx")

### Abundance ratio between aboveground and belowground parts
## hulled grain 7 
ratio.7.R <- fun.melt$SumAbund[which(fun.melt$Group == "H7R")]/fun.melt.total$TotalAbund[which(fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 7)]
ratio.7 <- fun.melt$SumAbund[which(fun.melt$Group == "H7Sh")]/fun.melt.total$TotalAbund[which(fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 7)]

ratio.14.R <- fun.melt$SumAbund[which(fun.melt$Group == "H14R")]/fun.melt.total$TotalAbund[which(fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 14)]
ratio.14 <- (fun.melt$SumAbund[which(fun.melt$Group == "H14S")]+fun.melt$SumAbund[which(fun.melt$Group == "H14L")])/fun.melt.total$TotalAbund[which(fun.melt.total$Hull == "hulled" & fun.melt.total$Age == 14)]


### unhulled grains
ratio.7.R.U <- fun.melt$SumAbund[which(fun.melt$Group == "U7R")]/fun.melt.total$TotalAbund[which(fun.melt.total$Hull == "intact" & fun.melt.total$Age == 7)]
ratio.7.U <- fun.melt$SumAbund[which(fun.melt$Group == "U7Sh")]/fun.melt.total$TotalAbund[which(fun.melt.total$Hull == "intact" & fun.melt.total$Age == 7)]

ratio.14.R.U <- fun.melt$SumAbund[which(fun.melt$Group == "U14R")]/fun.melt.total$TotalAbund[which(fun.melt.total$Hull == "intact" & fun.melt.total$Age == 14)]
ratio.14.U <- (fun.melt$SumAbund[which(fun.melt$Group == "U14S")]+fun.melt$SumAbund[which(fun.melt$Group == "U14L")])/fun.melt.total$TotalAbund[which(fun.melt.total$Hull == "intact" & fun.melt.total$Age == 14)]



### Diversity - Shannon and Chao1 or ACE
f.meta <- f.map
tab.fun.alpha.rarefy<-microbiome::alpha(fun.rarefy, index = "chao1")
f.meta$Richness2 <- tab.fun.alpha.rarefy$chao1
tab.fun.alpha.rarefy<-microbiome::alpha(fun.rarefy, index = "shannon")

f.meta$Shannon2 <- tab.fun.alpha.rarefy$diversity_shannon


## Richness of below and aboveground compartments
## hulled
chao1.below.7 <- mean(f.meta$Richness2[which(f.meta$Group == "H7R")])
chao1.above.7 <- mean(f.meta$Richness2[which(f.meta$Group =="H7Sh")])

chao1.below.14 <- mean(f.meta$Richness2[which(f.meta$Group == "H14R")])
chao1.above.14 <- mean(f.meta$Richness2[which(f.meta$Group =="H14L" | f.meta$Group =="H14S")])

## unhulled
chao1.below.7 <- mean(f.meta$Richness2[which(f.meta$Group == "U7R")])
chao1.above.7 <- mean(f.meta$Richness2[which(f.meta$Group =="U7Sh")])

chao1.below.14 <- mean(f.meta$Richness2[which(f.meta$Group == "U14R")])
chao1.above.14 <- mean(f.meta$Richness2[which(f.meta$Group =="U14L" | f.meta$Group =="U14S")])


## Shannon of below and aboveground compartments
## hulled
Shannon.below.7 <- mean(f.meta$Shannon2[which(f.meta$Group == "H7R")])
Shannon.above.7 <- mean(f.meta$Shannon2[which(f.meta$Group =="H7Sh")])

Shannon.below.14 <- mean(f.meta$Shannon2[which(f.meta$Group == "H14R")])
Shannon.above.14 <- mean(f.meta$Shannon2[which(f.meta$Group =="H14L" | f.meta$Group =="H14S")])

## unhulled
Shannon.below.7 <- mean(f.meta$Shannon2[which(f.meta$Group == "U7R")])
Shannon.above.7 <- mean(f.meta$Shannon2[which(f.meta$Group =="U7Sh")])

Shannon.below.14 <- mean(f.meta$Shannon2[which(f.meta$Group == "U14R")])
Shannon.above.14 <- mean(f.meta$Shannon2[which(f.meta$Group =="U14L" | f.meta$Group =="U14S")])

# wilcoxon test
x <- f.meta$Shannon2[which(f.meta$Group == "U14R")]
y <- f.meta$Shannon2[which(f.meta$Group =="U14L" | f.meta$Group =="U14S")]
wilcox.test(x, y, conf.int = TRUE) #0.1667

x <- f.meta$Shannon2[which(f.meta$Group == "H14R")]
y <- f.meta$Shannon2[which(f.meta$Group =="H14L" | f.meta$Group =="H14S")]
wilcox.test(x, y, conf.int = TRUE) #0.5476

x <- f.meta$Shannon2[which(f.meta$Group == "U7R")]
y <- f.meta$Shannon2[which(f.meta$Group =="U7Sh")]
wilcox.test(x, y, conf.int = TRUE) #0.3333

x <- f.meta$Shannon2[which(f.meta$Group == "H7R")]
y <- f.meta$Shannon2[which(f.meta$Group =="H7Sh")]
wilcox.test(x, y, conf.int = TRUE) #1




#### plotting
bac.melt
bac.melt$Group2 <- c(rep("H0",1), rep("H1",1), rep("H14",3), rep("H4",1),rep("H7",2),
                       rep("U0",1), rep("U1",1), rep("U14",3), rep("U4",1),rep("U7",2))

bac.melt$Compartment <- c(rep("Grain",2), rep("Leaf",1), rep("Root",1), rep("Stem",1), rep("Grain",1),rep("Root",1),rep("Shoot",1),
                          rep("Grain",2), rep("Leaf",1), rep("Root",1), rep("Stem",1), rep("Grain",1),rep("Root",1),rep("Shoot",1))
bac.melt.rel <- bac.melt %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = SumAbund*100/TotalAbund)  # Transform to rel. abundance

bac.melt.rel.sub <- subset(bac.melt.rel, Group2 %in% c("H7", "H14","U7","U14"))


ord.b <-bac.melt.rel.sub %>% group_by(Compartment)%>% summarise(CompartAbund = sum(SumAbund)) %>%  arrange(CompartAbund)
vec.b <- ord.b$Compartment

bac.melt.rel.sub$Compartment <- factor(bac.melt.rel.sub$Compartment, levels = vec.b)
bac.melt.rel.sub$Group2 <- factor(bac.melt.rel.sub$Group2, levels=c("H7", "H14","U7","U14"))
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
bac.melt.rel.sub.p1 <- ggplot(bac.melt.rel.sub, aes(x=Group2, y = RelAbundance, fill = Compartment)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  # scale_fill_manual(values = c("Gammaproteobacteria" = "darkolivegreen3", "Alphaproteobacteria"= "darkolivegreen", 
  #                              "Actinobacteria"="indianred2","Deltaproteobacteria" ="darkolivegreen1", 
  #                              "Bacteroidetes"="steelblue1", "Firmicutes" ="tan1", "Acidobacteria"="lightsalmon4", 
  #                              "Chloroflexi"="gold1", "Verrucomicrobia"="orchid3", "Nitrospirae"="palevioletred2",
  #                              "Planctomycetes" = "seagreen3", "Cyanobacteria" = "chartreuse4", "Epsilonbacteraeota" = "darkslategray4", 
  #                              "Spirochaetes" = "bisque4", "Deinococcus-Thermus" = "lightpink3", "Patescibacteria" = "lightblue1", 
  #                              "Lentisphaerae" = "lightgoldenrod3", "Dependentiae" = "chocolate3",
  #                              "Fibrobacteres" = "sienna2", "Armatimonadetes" = "coral", "Nitrospinae" = "darkkhaki", "Chlamydiae" = "darkslateblue",
  #                              "Kiritimatiellaeota" = "hotpink3", "Rokubacteria" = "lightsteelblue3", "Elusimicrobia" = "turquoise3",
  #                              "Fusobacteria" = "mistyrose3", "BRC1" = "plum2", "Latescibacteria" = "darkorchid2", "FBP" = "burlywood2",
  #                              "WPS-2" = "cornflowerblue", "Tenericutes" = "paleturquoise3", "Zixibacteria" = "#00CC99",
  #                              "Gemmatimonadetes"= "peachpuff3", "Low abundance" = "light grey", "unidentified" = "black")) +
  # 
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.5)
bac.melt.rel.sub.p1


fun.melt
fun.melt$Group2 <- c(rep("H0",1), rep("H1",1), rep("H14",3), rep("H4",1),rep("H7",2),
                     rep("U0",1), rep("U1",1), rep("U14",3), rep("U4",1),rep("U7",2))

fun.melt$Compartment <- c(rep("Grain",2), rep("Leaf",1), rep("Root",1), rep("Stem",1), rep("Grain",1),rep("Root",1),rep("Shoot",1),
                          rep("Grain",2), rep("Leaf",1), rep("Root",1), rep("Stem",1), rep("Grain",1),rep("Root",1),rep("Shoot",1))
fun.melt.rel <- fun.melt %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = SumAbund*100/TotalAbund)  # Transform to rel. abundance

fun.melt.rel.sub <- subset(fun.melt.rel, Group2 %in% c("H7", "H14","U7","U14"))


ord.f <-fun.melt.rel.sub %>% group_by(Compartment)%>% summarise(CompartAbund = sum(SumAbund)) %>%  arrange(CompartAbund)
vec.f <- ord.f$Compartment

fun.melt.rel.sub$Compartment <- factor(fun.melt.rel.sub$Compartment, levels = vec.f)
fun.melt.rel.sub$Group2 <- factor(fun.melt.rel.sub$Group2, levels=c("H7", "H14","U7","U14"))
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
fun.melt.rel.sub.p1 <- ggplot(fun.melt.rel.sub, aes(x=Group2, y = RelAbundance, fill = Compartment)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  # scale_fill_manual(values = c("Gammaproteobacteria" = "darkolivegreen3", "Alphaproteobacteria"= "darkolivegreen", 
  #                              "Actinobacteria"="indianred2","Deltaproteobacteria" ="darkolivegreen1", 
  #                              "Bacteroidetes"="steelblue1", "Firmicutes" ="tan1", "Acidobacteria"="lightsalmon4", 
  #                              "Chloroflexi"="gold1", "Verrucomicrobia"="orchid3", "Nitrospirae"="palevioletred2",
  #                              "Planctomycetes" = "seagreen3", "Cyanobacteria" = "chartreuse4", "Epsilonbacteraeota" = "darkslategray4", 
  #                              "Spirochaetes" = "bisque4", "Deinococcus-Thermus" = "lightpink3", "Patescibacteria" = "lightblue1", 
  #                              "Lentisphaerae" = "lightgoldenrod3", "Dependentiae" = "chocolate3",
  #                              "Fibrobacteres" = "sienna2", "Armatimonadetes" = "coral", "Nitrospinae" = "darkkhaki", "Chlamydiae" = "darkslateblue",
  #                              "Kiritimatiellaeota" = "hotpink3", "Rokubacteria" = "lightsteelblue3", "Elusimicrobia" = "turquoise3",
  #                              "Fusobacteria" = "mistyrose3", "BRC1" = "plum2", "Latescibacteria" = "darkorchid2", "FBP" = "burlywood2",
#                              "WPS-2" = "cornflowerblue", "Tenericutes" = "paleturquoise3", "Zixibacteria" = "#00CC99",
#                              "Gemmatimonadetes"= "peachpuff3", "Low abundance" = "light grey", "unidentified" = "black")) +
# 
xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.5)
fun.melt.rel.sub.p1




sort(colSums(otu_table(bac.clean.ss)))
sort(colSums(otu_table(fun.clean.ss)))
