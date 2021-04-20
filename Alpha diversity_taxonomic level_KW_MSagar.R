library(egg)
### Alpha diversity
###Rarefying
bac.rarefy
fun.rarefy

bac.rarefy.f <- prune_taxa(taxa_sums(bac.rarefy) > 0, bac.rarefy)
#tab_all <- microbiome::alpha(phy.rarefy.2017.t.oligo.t, index = "all")
#kable(head(tab_all))
#write.table(tab_all, "Alpha diversity_all.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")
ps1.meta <- data.frame(sample_data(bac.clean.ss))
ps1.meta

tab_shannon <- microbiome::alpha(bac.rarefy.f, index = "shannon")
tab_simpson <- microbiome::alpha(bac.rarefy.f, index = "evenness_simpson")
tab_observed <- microbiome::alpha(bac.rarefy.f, index = "observed")

ps1.meta$observed <- tab_observed$observed
ps1.meta$shannon <- tab_shannon$diversity_shannon
ps1.meta$simpson <- tab_simpson$evenness_simpson


ps1.meta$Group <- factor(ps1.meta$Group, levels = order.sample)

max.diversity <- aggregate(ps1.meta$observed, by = list(ps1.meta$Group), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(observed ~ Group, data = ps1.meta)
kw$p.value

#library(FSA)
DT = dunnTest(observed ~ Group,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=Group, y=observed)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 0.4)

p




write.csv(PT, "Statistical analysis_richness_Bacteria.csv")



##ANOVA
max.diversity <- aggregate(ps1.meta$shannon, by = list(ps1.meta$Group), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(shannon ~ Group, data = ps1.meta)
kw$p.value
#kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(shannon ~ Group,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=Group, y=shannon)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 0.4)

p

write.csv(PT, "Statistical analysis_Shannon_Bacteria.csv")

##ANOVA
max.diversity <- aggregate(ps1.meta$simpson, by = list(ps1.meta$Group), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simpson ~ Group, data = ps1.meta)
kw$p.value
#kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simpson ~ Group,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=Group, y=simpson)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("simpson\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")

p



###Fungi
fun.rarefy.f <- prune_taxa(taxa_sums(fun.rarefy) > 0, fun.rarefy)
#tab_all <- microbiome::alpha(phy.rarefy.2017.t.oligo.t, index = "all")
#kable(head(tab_all))
#write.table(tab_all, "Alpha diversity_all.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")
ps2.meta <- data.frame(sample_data(fun.clean.ss))
ps2.meta

tab_shannon <- microbiome::alpha(fun.rarefy.f, index = "shannon")
tab_simpson <- microbiome::alpha(fun.rarefy.f, index = "evenness_simpson")
tab_observed <- microbiome::alpha(fun.rarefy.f, index = "observed")

ps2.meta$observed <- tab_observed$observed
ps2.meta$shannon <- tab_shannon$diversity_shannon
ps2.meta$simpson <- tab_simpson$evenness_simpson

ps2.meta$Group <-  factor(ps2.meta$Group, levels = order.sample)

max.diversity <- aggregate(ps2.meta$observed, by = list(ps2.meta$Group), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(observed ~ Group, data = ps2.meta)
kw$p.value

#library(FSA)
DT = dunnTest(observed ~ Group,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=Group, y=observed)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Observed OTU\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 0.4)

p

write.csv(PT, "Statistical analysis_richness_Fungi.csv")



##ANOVA
max.diversity <- aggregate(ps2.meta$shannon, by = list(ps2.meta$Group), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(shannon ~ Group, data = ps2.meta)
kw$p.value
#kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(shannon ~ Group,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=Group, y=shannon)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Shannon\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 0.4)

p



write.csv(PT, "Statistical analysis_Shannon_Fungi.csv")

##ANOVA
max.diversity <- aggregate(ps2.meta$simpson, by = list(ps2.meta$Group), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simpson ~ Group, data = ps2.meta)
kw$p.value
#kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(simpson ~ Group,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=Group, y=simpson)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("simpson\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")

p


shapiro.test(ps1.meta$observed)
hist(ps1.meta$observed)

shapiro.test(ps1.meta$shannon)
hist(ps1.meta$shannon)

shapiro.test(ps2.meta$observed)
hist(ps2.meta$observed)

shapiro.test(ps2.meta$shannon)
hist(ps2.meta$shannon)
