### Read counts of each sample
order.sample <- c("H0",'H1',"H4", "H7Sh","H7R", "H14L","H14S","H14R","U0",'U1',"U4", "U7Sh","U7R", "U14L","U14S","U14R")
## Bacteria
df.sample <- bac.clean.ss %>%
  psmelt()
head(df.sample)
# we need to group by samples
df.sample.total <- df.sample %>%  
  group_by(Sample, Group) %>%                         # Filter out at absolute read of 20       
  summarise(Total = sum(Abundance))  # Transform to rel. abundance
df.sample.total$log_total <- log(df.sample.total$Total)

df.sample.total$Group <- factor(df.sample.total$Group, levels = order.sample)

### Plotting
library(ggpubr)
ggbarplot(df.sample.total, x = "Group", y = "Total",
          add = c("mean_se", "jitter"),
          position = position_dodge(0.8))+theme(aspect.ratio = 0.5)
write.csv(df.sample.total, "df.sample.total.csv")

library(data.table) 

read.tab<-dcast(setDT(df.sample.total), Phylum2 ~ Group, value.var = "RelAbundance")


##Fungi
df.sample <- fun.clean.ss %>%
  psmelt()
head(df.sample)
# we need to group by samples
df.sample.total <- df.sample %>%  
  group_by(Sample, Group) %>%                         # Filter out at absolute read of 20       
  summarise(Total = sum(Abundance))  # Transform to rel. abundance
df.sample.total$log_total <- log(df.sample.total$Total)

df.sample.total$Group <- factor(df.sample.total$Group, levels = order.sample)

### Plotting
library(ggpubr)
ggbarplot(df.sample.total, x = "Group", y = "Total",
          add = c("mean_se", "jitter"),
          position = position_dodge(0.8))+theme(aspect.ratio = 0.5)
write.csv(df.sample.total, "df.sample.total_fungi.csv")




### Relative abundance table
##Bacteria
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

# we need to group by samples
df.phylum.rel <- df.phylum.2 %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

ord <- df.phylum.rel %>% group_by(Phylum2, Group) %>% summarise(RelAbundance = sum(RelAbundance))
df.phylum.rel$RelAbundance

library(data.table) 
ord.tab<-dcast(setDT(ord), Phylum2 ~ Group, value.var = "RelAbundance")

write.csv(ord.tab, "Relative abundance of bacterial phyla_Group.csv")
##Fungi
df.class.fun<- fun.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.class.fun %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
levels(df.class.fun$Class) = c(levels(df.class.fun$Class), 'Low abundance')

# we need to group by samples
df.class.fun.rel <- df.class.fun %>%  
  group_by(Group) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

ord.f <- df.class.fun.rel %>% group_by(Class, Group) %>% summarise(RelAbundance = sum(RelAbundance))

ord.f.tab<-dcast(setDT(ord.f), Class ~ Group, value.var = "RelAbundance")

write.csv(ord.f.tab, "Relative abundance of fungal class_Group.csv")







