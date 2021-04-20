### hulled vs unhulled
###
## Defining daOTUs
#Phyloseq files
bac.clean.ss
fun.clean.ss

## Subgrouping
bac.clean.ss.grain <- subset_samples(bac.clean.ss, Compartment == "Grain")
bac.clean.ss.leaf <- subset_samples(bac.clean.ss, Compartment %in% c("Leaf","Stem","Shoot") )


zero.clean.filt <- phyloseq::filter_taxa(bac.clean.ss.leaf, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(zero.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(zero.clean.filt))

obj <- phyloseq_to_metagenomeSeq(zero.clean.filt)
obj 

## fitZig sample
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 30, verbose = TRUE)

Type  <-  pData(obj)$Hull

mod  <-  model.matrix(~Type)
colnames(mod)  <-  levels(Type)
colnames(mod)

res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)


zigFit = res$fit
finalMod= res$fit$design


contrast.matrix.1 =makeContrasts(hulled - intact, levels = finalMod)


fit2_hull = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 

fit3.1 = eBayes(fit2_hull)

topTable(fit3.1, coef="hulled - intact")

res.1 <- topTable(fit3.1,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =2)
head(res.1)
dim(res.1)

## Taxonomy
res.1$OTU <- rownames(res.1)
bac.list.otu_id


res.tax_bac <- merge(res.1,bac.list, by = "OTU")

write.xlsx(res.tax_bac, 'bac_daOTU_hull_leaf.xlsx')


## MA plot
log2AverageAbundance <- psmelt(bac.clean.ss) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
log2AverageAbundance
Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
Ta <- unique(Ta)
Ta
resSig = res[!is.na(res$adj.P.Val), ]
resSig = data.frame(resSig)
head(resSig)
resSig <- tibble::rownames_to_column(resSig, 'OTU')
resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))


#### fungal community ###
fun.clean.ss.grain <- subset_samples(fun.clean.ss, Organ == "Grain")
fun.clean.ss.leaf <- subset_samples(fun.clean.ss, Organ %in% c("Leaf","Shoot","Stem"))
zero.clean.filt <- phyloseq::filter_taxa(fun.clean.ss.leaf, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(zero.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(zero.clean.filt))

obj <- phyloseq_to_metagenomeSeq(zero.clean.filt)

## fitZig sample
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 30, verbose = TRUE)

Type  <-  pData(obj)$Hull

mod  <-  model.matrix(~Type)
colnames(mod)  <-  levels(Type)
colnames(mod)

res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)

class(res)

zigFit = res$fit
finalMod= res$fit$design


contrast.matrix.1 =makeContrasts(hulled - intact, levels = finalMod)


fit2_hull = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 

fit3.1 = eBayes(fit2_hull)

topTable(fit3.1, coef="hulled - intact")

res.1 <- topTable(fit3.1,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =2)
head(res.1)
dim(res.1)

## Taxonomy
res.1$OTU <- rownames(res.1)

res.tax_fun <- merge(res.1,fun.list, by = "OTU")

write.xlsx(res.tax_fun, 'fun_daOTU_hull_leaf.xlsx')


## MA plot
log2AverageAbundance <- psmelt(fun.clean.ss) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
log2AverageAbundance
Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
Ta <- unique(Ta)
Ta
resSig = res[!is.na(res$adj.P.Val), ]
resSig = data.frame(resSig)
head(resSig)
resSig <- tibble::rownames_to_column(resSig, 'OTU')
resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))





### common OTUs in hulled and unhulled seeds
taxa_names(bac.clean.ss.U0)
taxa_names(bac.clean.ss.H0)

intersect(taxa_names(bac.clean.ss.U0), taxa_names(bac.clean.ss.H0)) #16/52(0.3076923)   16/93(0.172043)
intersect(taxa_names(bac.clean.ss.U4), taxa_names(bac.clean.ss.H4))# 5/31(0.1612903)     5/28(0.1785714)
intersect(taxa_names(bac.clean.ss.U7R), taxa_names(bac.clean.ss.H7R))#5/19(0.2631579)    5/40(0.125)
intersect(taxa_names(bac.clean.ss.U7Sh), taxa_names(bac.clean.ss.H7Sh))# 7/17(0.4117647)  7/20(0.35)
intersect(taxa_names(bac.clean.ss.U14R), taxa_names(bac.clean.ss.H14R))#11/22(0.5)       11/122(0.09016393)
intersect(taxa_names(bac.clean.ss.U14S), taxa_names(bac.clean.ss.H14S))#12/33(0.3636364)  12/17(0.7058824)
intersect(taxa_names(bac.clean.ss.U14L), taxa_names(bac.clean.ss.H14L))#8/28(0.2857143)   8/15(0.5333333)


intersect(taxa_names(fun.clean.ss.U0), taxa_names(fun.clean.ss.H0))#11/26(0.4230769)   11/35(0.3142857)
intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.H4))#11/24(0.4583333)   11/20(0.55)
intersect(taxa_names(fun.clean.ss.U7R), taxa_names(fun.clean.ss.H7R))#8/20(0.4)       8/22(0.3636364)
intersect(taxa_names(fun.clean.ss.U7Sh), taxa_names(fun.clean.ss.H7Sh))#3/15(0.2)     3/11(0.2727273)
intersect(taxa_names(fun.clean.ss.U14R), taxa_names(fun.clean.ss.H14R))#8/26(0.3076923)   8/26(0.3076923)
intersect(taxa_names(fun.clean.ss.U14S), taxa_names(fun.clean.ss.H14S))#3/26(0.1153846)   3/9(0.3333333)
intersect(taxa_names(fun.clean.ss.U14L), taxa_names(fun.clean.ss.H14L))#4/16(0.25)   4/16(0.25)


### Relative abundance of common OTUs in each compartment
bac.melt <- bac.clean.ss %>% psmelt() %>% group_by(Group)
bac.melt.sub <- subset(bac.melt, Group %in% c("U0", "H0"))
totals<-bac.melt.sub %>% group_by(Group) %>% summarise(totalabund = sum(Abundance))
bac.melt.sub$TotalAbund <- 0
bac.melt.sub$TotalAbund <- ifelse(bac.melt.sub$Group == "U0", totals$totalabund[totals$Group == "U0"], totals$totalabund[totals$Group == "H0"])

bac.melt.sub.rel <- bac.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

bac.melt.sub.rel <- subset(bac.melt.sub.rel, OTU %in% intersect(taxa_names(bac.clean.ss.U0), taxa_names(bac.clean.ss.H0)))

bac.melt.sub.rel.p1 <- ggplot(bac.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
bac.melt.sub.rel.p1


bac.melt.sub <- subset(bac.melt, Group %in% c("U4", "H4"))
bac.melt.sub.rel <- bac.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

bac.melt.sub.rel <- subset(bac.melt.sub.rel, OTU %in% intersect(taxa_names(bac.clean.ss.U4), taxa_names(bac.clean.ss.H4)))

bac.melt.sub.rel.p1 <- ggplot(bac.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
bac.melt.sub.rel.p1



bac.melt.sub <- subset(bac.melt, Group %in% c("U7Sh", "H7Sh"))
bac.melt.sub.rel <- bac.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

bac.melt.sub.rel <- subset(bac.melt.sub.rel, OTU %in% intersect(taxa_names(bac.clean.ss.U7Sh), taxa_names(bac.clean.ss.H7Sh)))

bac.melt.sub.rel.p1 <- ggplot(bac.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
bac.melt.sub.rel.p1


bac.melt.sub <- subset(bac.melt, Group %in% c("U7R", "H7R"))
bac.melt.sub.rel <- bac.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

bac.melt.sub.rel <- subset(bac.melt.sub.rel, OTU %in% intersect(taxa_names(bac.clean.ss.U7R), taxa_names(bac.clean.ss.H7R)))

bac.melt.sub.rel.p1 <- ggplot(bac.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
bac.melt.sub.rel.p1


bac.melt.sub <- subset(bac.melt, Group %in% c("U14R", "H14R"))
bac.melt.sub.rel <- bac.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

bac.melt.sub.rel <- subset(bac.melt.sub.rel, OTU %in% intersect(taxa_names(bac.clean.ss.U14R), taxa_names(bac.clean.ss.H14R)))

bac.melt.sub.rel.p1 <- ggplot(bac.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
bac.melt.sub.rel.p1


bac.melt.sub <- subset(bac.melt, Group %in% c("U14S", "H14S"))
bac.melt.sub.rel <- bac.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

bac.melt.sub.rel <- subset(bac.melt.sub.rel, OTU %in% intersect(taxa_names(bac.clean.ss.U14S), taxa_names(bac.clean.ss.H14S)))

bac.melt.sub.rel.p1 <- ggplot(bac.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
bac.melt.sub.rel.p1


bac.melt.sub <- subset(bac.melt, Group %in% c("U14L", "H14L"))
bac.melt.sub.rel <- bac.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

bac.melt.sub.rel <- subset(bac.melt.sub.rel, OTU %in% intersect(taxa_names(bac.clean.ss.U14L), taxa_names(bac.clean.ss.H14L)))

bac.melt.sub.rel.p1 <- ggplot(bac.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
bac.melt.sub.rel.p1




### Fungi
fun.melt <- fun.clean.ss %>% psmelt() %>% group_by(Group)
fun.melt.sub <- subset(fun.melt, Group %in% c("U0", "H0"))
fun.melt.sub.rel <- fun.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

fun.melt.sub.rel <- subset(fun.melt.sub.rel, OTU %in% intersect(taxa_names(fun.clean.ss.U0), taxa_names(fun.clean.ss.H0)))

fun.melt.sub.rel.p1 <- ggplot(fun.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
fun.melt.sub.rel.p1


fun.melt.sub <- subset(fun.melt, Group %in% c("U4", "H4"))
fun.melt.sub.rel <- fun.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

fun.melt.sub.rel <- subset(fun.melt.sub.rel, OTU %in% intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.H4)))

fun.melt.sub.rel.p1 <- ggplot(fun.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
fun.melt.sub.rel.p1



fun.melt.sub <- subset(fun.melt, Group %in% c("U7Sh", "H7Sh"))
fun.melt.sub.rel <- fun.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

fun.melt.sub.rel <- subset(fun.melt.sub.rel, OTU %in% intersect(taxa_names(fun.clean.ss.U7Sh), taxa_names(fun.clean.ss.H7Sh)))

fun.melt.sub.rel.p1 <- ggplot(fun.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
fun.melt.sub.rel.p1


fun.melt.sub <- subset(fun.melt, Group %in% c("U7R", "H7R"))
fun.melt.sub.rel <- fun.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

fun.melt.sub.rel <- subset(fun.melt.sub.rel, OTU %in% intersect(taxa_names(fun.clean.ss.U7R), taxa_names(fun.clean.ss.H7R)))

fun.melt.sub.rel.p1 <- ggplot(fun.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
fun.melt.sub.rel.p1


fun.melt.sub <- subset(fun.melt, Group %in% c("U14R", "H14R"))
fun.melt.sub.rel <- fun.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

fun.melt.sub.rel <- subset(fun.melt.sub.rel, OTU %in% intersect(taxa_names(fun.clean.ss.U14R), taxa_names(fun.clean.ss.H14R)))

fun.melt.sub.rel.p1 <- ggplot(fun.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
fun.melt.sub.rel.p1


fun.melt.sub <- subset(fun.melt, Group %in% c("U14S", "H14S"))
fun.melt.sub.rel <- fun.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

fun.melt.sub.rel <- subset(fun.melt.sub.rel, OTU %in% intersect(taxa_names(fun.clean.ss.U14S), taxa_names(fun.clean.ss.H14S)))

fun.melt.sub.rel.p1 <- ggplot(fun.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
fun.melt.sub.rel.p1


fun.melt.sub <- subset(fun.melt, Group %in% c("U14L", "H14L"))
fun.melt.sub.rel <- fun.melt.sub %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))

fun.melt.sub.rel <- subset(fun.melt.sub.rel, OTU %in% intersect(taxa_names(fun.clean.ss.U14L), taxa_names(fun.clean.ss.H14L)))

fun.melt.sub.rel.p1 <- ggplot(fun.melt.sub.rel, aes(x=Group, y = RelAbundance, fill = Group)) + 
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
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+ylim(0,100)+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
fun.melt.sub.rel.p1
