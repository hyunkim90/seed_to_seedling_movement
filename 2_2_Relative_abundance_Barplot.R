
###### 2. make barplot with phyloseq

order.sample <- c("H0",'H1',"H4", "H7Sh","H7R", "H14L","H14S","H14R","U0",'U1',"U4", "U7Sh","U7R", "U14L","U14S","U14R")
## Phylum level
df.phylum <- bac.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.phylum$Phylum <- as.character(df.phylum$Phylum)
df.phylum$Phylum2 <- df.phylum$Phylum
df.phylum$Phylum2[which(df.phylum$Class=="Alphaproteobacteria")] <- "Alphaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Gammaproteobacteria")] <- "Gammaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Deltaproteobacteria")] <- "Deltaproteobacteria"

unique(df.phylum$Phylum2)
head(df.phylum)

library(forcats) 
df.phylum %<>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))

df.phylum.2 <- df.phylum %>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))

levels(df.phylum$Phylum2)
levels(df.phylum.2$Phylum2) = c(levels(df.phylum.2$Phylum2), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum.2 %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 0.5,]$Phylum2 <- 'Low abundance'
unique(df.phylum$Phylum2)

ord <- df.phylum.rel %>% group_by(Phylum2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Phylum2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac %in% c("Low abundance", "Alphaproteobacteria", "Gammaproteobacteria","Deltaproteobacteria","unidentified"))]
vec.proteo <- c("Deltaproteobacteria","Alphaproteobacteria", "Gammaproteobacteria")
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.order)
vec.reorder <- append(vec.reorder, vec.proteo)

df.phylum.rel$Phylum2 <- factor(df.phylum.rel$Phylum2, levels = vec.reorder) 
df.phylum.rel$Group <- factor(df.phylum.rel$Group, levels=order.sample)
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=Group, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Gammaproteobacteria" = "darkolivegreen3", "Alphaproteobacteria"= "darkolivegreen", 
                               "Actinobacteria"="indianred2","Deltaproteobacteria" ="darkolivegreen1", 
                               "Bacteroidetes"="steelblue1", "Firmicutes" ="tan1", "Acidobacteria"="lightsalmon4", 
                               "Chloroflexi"="gold1", "Verrucomicrobia"="orchid3", "Nitrospirae"="palevioletred2",
                               "Planctomycetes" = "seagreen3", "Cyanobacteria" = "chartreuse4", "Epsilonbacteraeota" = "darkslategray4", 
                               "Spirochaetes" = "bisque4", "Deinococcus-Thermus" = "lightpink3", "Patescibacteria" = "lightblue1", 
                               "Lentisphaerae" = "lightgoldenrod3", "Dependentiae" = "chocolate3",
                               "Fibrobacteres" = "sienna2", "Armatimonadetes" = "coral", "Nitrospinae" = "darkkhaki", "Chlamydiae" = "darkslateblue",
                               "Kiritimatiellaeota" = "hotpink3", "Rokubacteria" = "lightsteelblue3", "Elusimicrobia" = "turquoise3",
                               "Fusobacteria" = "mistyrose3", "BRC1" = "plum2", "Latescibacteria" = "darkorchid2", "FBP" = "burlywood2",
                               "WPS-2" = "cornflowerblue", "Tenericutes" = "paleturquoise3", "Zixibacteria" = "#00CC99",
                               "Gemmatimonadetes"= "peachpuff3", "Low abundance" = "light grey", "unidentified" = "black")) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,200))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.5)
df.phylum.rel.p1


##Class level
df.class <- bac.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.class %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))


levels(df.class$Class)
levels(df.class$Class) = c(levels(df.class$Class), 'Low abundance')

# we need to group by samples
df.class.rel <- df.class %>%  
  group_by(Diagnosis) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.class.rel[df.class.rel$RelAbundance < 0.05,]$Class <- 'Low abundance'
unique(df.class$Class)

ord <- df.class.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Class
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac =="Low abundance")]
vec.Low <- c("Low abundance")
vec.reorder <- append(vec.Low, vec.order)


df.class.rel$Class <- factor(df.class.rel$Class, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.class.rel.p1 <- ggplot(df.class.rel, aes(x=Diagnosis, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_Class) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+theme(aspect.ratio=2)
df.class.rel.p1


###Order
df.order <- bac.clean.ss %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.order %<>% mutate(Order = fct_explicit_na(Order, na_level = "unidentified"))


levels(df.order$Order)
levels(df.order$Order) = c(levels(df.order$Order), 'Low abundance')

# we need to group by samples
df.order.rel <- df.order %>%  
  group_by(Diagnosis) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.order.rel[df.order.rel$RelAbundance < 0.05,]$Order <- 'Low abundance'
unique(df.order$Order)

ord <- df.order.rel %>% group_by(Order) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Order
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.order)


df.order.rel$Order <- factor(df.order.rel$Order, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.order.rel.p1 <- ggplot(df.order.rel, aes(x=Diagnosis, y = RelAbundance, fill = Order)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_Order) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Order Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 2)

df.order.rel.p1


##Family
df.family <- phy.clean.ss %>%
  tax_glom(taxrank = "Family", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.family %<>% mutate(Family = fct_explicit_na(Family, na_level = "unidentified"))


levels(df.family$Family)
levels(df.family$Family) = c(levels(df.family$Family), 'Low abundance')

# we need to group by samples
df.family.rel <- df.family %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.family.rel[df.family.rel$RelAbundance < 1,]$Family <- 'Low abundance'
unique(df.family$Family)

ord <- df.family.rel %>% group_by(Family) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Family
vec.charac<-as.character(vec)
vec.family <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder.fam <- append(vec.Low, vec.family)


df.family.rel$Family <- factor(df.family.rel$Family, levels = vec.reorder.fam) 
df.family.rel$Group <- factor(df.family.rel$Group, levels = order.sample) 
df.family.rel <- subset(df.family.rel, Group %in% c("H0","H4","H14L","H14S","H14R","U0","U4","U14L","U14S","U14R"))

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.family.rel.p1 <- ggplot(df.family.rel, aes(x=Group, y = RelAbundance, fill = Family)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_Family) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Family Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 6,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 0.5)

df.family.rel.p1


### Genus
df.genus <- phy.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = F) %>% 
  psmelt()

df.genus %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))

# df.genus$Phylum<-as.character(df.genus$Phylum)
# df.genus$Class<-as.character(df.genus$Class)
# df.genus$Order<-as.character(df.genus$Order)
# df.genus$Family<-as.character(df.genus$Family)
# df.genus$Genus<-as.character(df.genus$Genus)
# 
# df.genus$Phylum[is.na(df.genus$Phylum)] <- "unidentified"
# df.genus$Class[is.na(df.genus$Class)] <- "unidentified"
# df.genus$Order[is.na(df.genus$Order)] <- "unidentified"
# df.genus$Family[is.na(df.genus$Family)] <- "unidentified"
# df.genus$Genus[is.na(df.genus$Genus)] <- "unidentified"
# 
# df.genus$Genus2 <- ifelse(df.genus$Genus == "unidentified",ifelse(df.genus$Family=="unidentified",paste0(df.genus$Order,'_',"Unidentified genus"),paste0(df.genus$Family,'_',"Unidentified genus")),paste0(df.genus$Genus))
# 

levels(df.genus$Genus) = c(levels(df.genus$Genus), 'Low abundance')

# we need to group by samples
df.genus.rel <- df.genus %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance < 1,]$Genus <- 'Low abundance'
unique(df.genus$Genus)

ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.genus <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.genus)
df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.reorder) 
df.genus.rel$Group <- factor(df.genus.rel$Group, levels = order.sample)
df.genus.rel <- subset(df.genus.rel, Group %in% c("H0","H4","H14L","H14S","H14R","U0","U4","U14L","U14S","U14R"))
df.genus.rel$Group_num <- ifelse(df.genus.rel$Group == "U0",1,ifelse(df.genus.rel$Group == "U4", 2, ifelse(df.genus.rel$Group == "U14L",3,ifelse(df.genus.rel$Group == "U14S",4,5))))

df.genus.rel.2<-df.genus.rel %>% group_by(Genus,Group_num) %>% summarise(CumRelAbundance = sum(RelAbundance))
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=Group, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_OTU) +
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 6,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  #scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.5)

df.genus.rel.p1


################ Fungal community #######################
## Phylum level
##Grouped by diagnosis
df.phylum.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.phylum.fun %<>% mutate(Phylum = fct_explicit_na(Phylum, na_level = "unidentified"))
levels(df.phylum.fun$Phylum) = c(levels(df.phylum.fun$Phylum), 'Low abundance')

# we need to group by samples
df.phylum.fun.rel <- df.phylum.fun %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.phylum.fun.rel[df.phylum.fun.rel$RelAbundance < 0.05,]$Phylum <- 'Low abundance'
unique(df.phylum.fun$Phylum)

ord.f <- df.phylum.fun.rel %>% group_by(Phylum) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Phylum
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.phylum.fun.rel$Phylum <- factor(df.phylum.fun.rel$Phylum, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.fun.rel.p1 <- ggplot(df.phylum.fun.rel, aes(x=Group, y = RelAbundance, fill = Phylum)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Ascomycota" = "#11335F", "Basidiomycota"= "#BE4146",
                               "Mortierellomycota"="#E4AF2C","Mucoromycota"="#4E734E","Olpidiomycota"="#4E734E","unidentified" ="#000000",
                               "Low abundance"="#BFBEBE")) +

  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
df.phylum.fun.rel.p1



###class
##Grouped by diagnosis
df.class.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.class.fun %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
levels(df.class.fun$Class) = c(levels(df.class.fun$Class), 'Low abundance')

# we need to group by samples
df.class.fun.rel <- df.class.fun %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.class.fun.rel[df.class.fun.rel$RelAbundance < 0.5,]$Class <- 'Low abundance'
unique(df.class.fun$Class)

ord.f <- df.class.fun.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Class
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.class.fun.rel$Class <- factor(df.class.fun.rel$Class, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
df.class.fun.rel$Group <- factor(df.class.fun.rel$Group, levels = order.sample)
## plot relative abundance
#bar plot
df.class.fun.rel.p1 <- ggplot(df.class.fun.rel, aes(x=Group, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
scale_fill_manual(values = c("Mortierellomycetes" = "#4E734E", "Eurotiomycetes"= "#6DA9DC",
                             "Mucoromycetes"="#E4AF2C","Tremellomycetes"="#BE4146", "Microbotryomycetes" ="#DC9A9E",
                             "Dothideomycetes" = "#5195D1", "Agaricomycetes" = "#CC6C71", "Blastocladiomycetes" = "#87AC88",
                             "Sordariomycetes"= "#1E63AF","Leotiomycetes"= "#11335F", "Cystobasidiomycetes" = "#A871AE", "Pezizomycetes" = "#C0DBF3",
                             "unidentified" ="#000000", "Endogonomycetes" ="#CC9900", "Pucciniomycetes" = "#0099FF", "Ustilaginomycetes"="#ffcc33","Arthoniomycetes" = "#cccc99",
                             "Saccharomycetes" = "#666699", "Low abundance" = "light grey", "Malasseziomycetes" = "#99cc99", "Exobasidiomycetes" = "#663366", 
                             "Rhizophydiomycetes"="#ffccff","Wallemiomycetes"="#333300","Agaricostilbomycetes"="#003333","Orbiliomycetes"="#cc9966","Chytridiomycetes" = "#006600",
                             "Spizellomycetes" = "#66cc99","Moniliellomycetes" = "#0099cc","Monoblepharidomycetes" = "#CC0000", "Dacrymycetes" = "#669999",
                             "Spiculogloeomycetes" = "#99cc00")) +

  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,200))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.5)
df.class.fun.rel.p1



##Order
##Grouped by diagnosis
df.order.fun <- fun.its1.clean.ss %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.order.fun %<>% mutate(Order = fct_explicit_na(Order, na_level = "unidentified"))
levels(df.order.fun$Order) = c(levels(df.order.fun$Order), 'Low abundance')

# we need to group by samples
df.order.fun.rel <- df.order.fun %>%  
  group_by(Diagnosis) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.order.fun.rel[df.order.fun.rel$RelAbundance < 0.05,]$Order <- 'Low abundance'
unique(df.order.fun$Order)

ord.f <- df.order.fun.rel %>% group_by(Order) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Order
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.order.fun.rel$Order <- factor(df.order.fun.rel$Order, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.order.fun.rel.p1 <- ggplot(df.order.fun.rel, aes(x=Diagnosis, y = RelAbundance, fill = Order)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_Order) +
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Order Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=2)
df.order.fun.rel.p1


##Family
df.family.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Family", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.family.fun %<>% mutate(Family = fct_explicit_na(Family, na_level = "unidentified"))
levels(df.family.fun$Family) = c(levels(df.family.fun$Family), 'Low abundance')

# we need to group by samples
df.family.fun.rel <- df.family.fun %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.family.fun.rel[df.family.fun.rel$RelAbundance < 1,]$Family <- 'Low abundance'
unique(df.family.fun$Family)

ord.f <- df.family.fun.rel %>% group_by(Family) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Family
vec.f.charac<-as.character(vec.f)
vec.family.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.family.f)


df.family.fun.rel$Family <- factor(df.family.fun.rel$Family, levels = vec.reorder.f) 
df.family.fun.rel$Group <- factor(df.family.fun.rel$Group, levels = order.sample)
df.family.fun.rel <- subset(df.family.fun.rel, Group %in% c("H0","H4","H14L","H14S","H14R","U0","U4","U14L","U14S","U14R"))

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.family.fun.rel.p1 <- ggplot(df.family.fun.rel, aes(x=Group, y = RelAbundance, fill = Family)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_Family) +
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Family Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 5,reverse = T))+
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
df.family.fun.rel.p1

write.csv(df.family.rel, "Source_bacterial_MS agar_abundance.csv")

write.csv(df.family.fun.rel, "Source_fungal_MS agar_abundance.csv")
##Genus
df.genus.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.genus.fun %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))
# df.genus.fun$Phylum<-as.character(df.genus.fun$Phylum)
# df.genus.fun$Class<-as.character(df.genus.fun$Class)
# df.genus.fun$Order<-as.character(df.genus.fun$Order)
# df.genus.fun$Family<-as.character(df.genus.fun$Family)
# df.genus.fun$Genus<-as.character(df.genus.fun$Genus)
# 
# df.genus.fun$Phylum[is.na(df.genus.fun$Phylum)] <- "unidentified"
# df.genus.fun$Class[is.na(df.genus.fun$Class)] <- "unidentified"
# df.genus.fun$Order[is.na(df.genus.fun$Order)] <- "unidentified"
# df.genus.fun$Family[is.na(df.genus.fun$Family)] <- "unidentified"
# df.genus.fun$Genus[is.na(df.genus.fun$Genus)] <- "unidentified"
# 
# df.genus.fun$Genus2 <- ifelse(df.genus.fun$Genus == "unidentified",ifelse(df.genus.fun$Family=="unidentified",paste0(df.genus.fun$Order,'_',"Unidentified genus"),paste0(df.genus.fun$Family,'_',"Unidentified genus")),paste0(df.genus.fun$Genus))
# 

levels(df.genus.fun$Genus) = c(levels(df.genus.fun$Genus), 'Low abundance')

# we need to group by samples
df.genus.fun.rel <- df.genus.fun %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.genus.fun.rel[df.genus.fun.rel$RelAbundance < 1,]$Genus <- 'Low abundance'
unique(df.genus.fun$Genus)

ord.f <- df.genus.fun.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Genus
vec.f.charac<-as.character(vec.f)
vec.genus.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.genus.f)


df.genus.fun.rel$Genus <- factor(df.genus.fun.rel$Genus, levels = vec.reorder.f) 
df.genus.fun.rel$Group<- factor(df.genus.fun.rel$Group, levels = c("H0","H1","H4","H7Sh","H7R","H14L","H14S","H14R",
                                                                    "U0","U1","U4","U7Sh","U7R","U14L","U14S","U14R")) 
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.fun.rel.p1 <- ggplot(df.genus.fun.rel, aes(x=Group, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_gen) +
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 10,reverse = T))+
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
df.genus.fun.rel.p1