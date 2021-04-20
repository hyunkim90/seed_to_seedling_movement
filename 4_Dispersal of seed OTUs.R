### Dispersal of seed OTUs to compartments of seedlings
## OTUs in intact seeds

fun.clean.ss.U0 <- subset_samples(fun.clean.ss, Group == "U0")
fun.clean.ss.U0<-phyloseq::filter_taxa(fun.clean.ss.U0, function(x) sum(x) != 0, TRUE)

fun.clean.ss.U4 <- subset_samples(fun.clean.ss, Group == "U4")
fun.clean.ss.U4<-phyloseq::filter_taxa(fun.clean.ss.U4, function(x) sum(x) != 0, TRUE)

fun.clean.ss.U7Sh <- subset_samples(fun.clean.ss, Group == "U7Sh")
fun.clean.ss.U7Sh<-phyloseq::filter_taxa(fun.clean.ss.U7Sh, function(x) sum(x) != 0, TRUE)
fun.clean.ss.U7R <- subset_samples(fun.clean.ss, Group == "U7R")
fun.clean.ss.U7R<-phyloseq::filter_taxa(fun.clean.ss.U7R, function(x) sum(x) != 0, TRUE)

fun.clean.ss.U14L <- subset_samples(fun.clean.ss, Group == "U14L")
fun.clean.ss.U14L<-phyloseq::filter_taxa(fun.clean.ss.U14L, function(x) sum(x) != 0, TRUE)
fun.clean.ss.U14S <- subset_samples(fun.clean.ss, Group == "U14S")
fun.clean.ss.U14S<-phyloseq::filter_taxa(fun.clean.ss.U14S, function(x) sum(x) != 0, TRUE)
fun.clean.ss.U14R <- subset_samples(fun.clean.ss, Group == "U14R")
fun.clean.ss.U14R<-phyloseq::filter_taxa(fun.clean.ss.U14R, function(x) sum(x) != 0, TRUE)

Sh7<-intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.U7Sh)) #8
R7<-intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.U7R)) #8

fun.list.Sh7 <- subset(fun.list, OTU %in% Sh7)
fun.list.R7 <- subset(fun.list, OTU %in% R7)

Sh7<-intersect(taxa_names(fun.clean.ss.U0), taxa_names(fun.clean.ss.U7Sh)) #9
R7<-intersect(taxa_names(fun.clean.ss.U0), taxa_names(fun.clean.ss.U7R)) #3

Sh7.4<-intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.U7Sh)) #6
R7.4<-intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.U7R)) #5


zero_to_7.Sh<-intersect(Sh7, Sh7.4) #5
zero_to_7.R<-intersect(R7, R7.4) #3

fun.list.Sh7 <- subset(fun.list, OTU %in% zero_to_7.Sh)
fun.list.R7 <- subset(fun.list, OTU %in% zero_to_7.R)

### Seedlings 14-day-old
L14<-intersect(taxa_names(fun.clean.ss.U0), taxa_names(fun.clean.ss.U14L)) #9
S14<-intersect(taxa_names(fun.clean.ss.U0), taxa_names(fun.clean.ss.U14S)) #13
R14<-intersect(taxa_names(fun.clean.ss.U0), taxa_names(fun.clean.ss.U14R)) #3

L14.4<-intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.U14L)) #6
S14.4<-intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.U14S)) #8
R14.4<-intersect(taxa_names(fun.clean.ss.U4), taxa_names(fun.clean.ss.U14R)) #5


fun.list.L14 <- subset(fun.list, OTU %in% L14)
fun.list.S14 <- subset(fun.list, OTU %in% S14)
fun.list.R14 <- subset(fun.list, OTU %in% R14)

sort(colSums(otu_table(fun.clean.ss)))

zero_to_14.L<-intersect(L14, L14.4) #5
zero_to_14.S <-intersect(S14, S14.4) #7
zero_to_14.R<-intersect(R14, R14.4) #3

fun.list.L14 <- subset(fun.list, OTU %in% zero_to_14.L)
fun.list.S14 <- subset(fun.list, OTU %in% zero_to_14.S)
fun.list.R14 <- subset(fun.list, OTU %in% zero_to_14.R)


dispersal.L<-intersect(zero_to_7.Sh, zero_to_14.L)
dispersal.S<-intersect(zero_to_7.Sh, zero_to_14.S)
dispersal.R<-intersect(zero_to_7.R, zero_to_14.R)

compartment = c("U0", "U4", "U7Sh", "U14L")
dispersal = dispersal.L
df.ridge_2.fun.dispersal.leaf<- get_df_ridge(fun.clean.ss, dispersal = dispersal.L)

compartment = c("U0", "U4", "U7Sh", "U14S")
df.ridge_2.fun.dispersal.stem<- get_df_ridge(fun.clean.ss, dispersal = dispersal.S)

compartment = c("U0", "U4", "U7R", "U14R")
df.ridge_2.fun.dispersal.root<- get_df_ridge(fun.clean.ss, dispersal = dispersal.R)



df.ridge_2.fun.dispersal.leaf$Sample <- factor(df.ridge_2.fun.dispersal.leaf$Sample, levels = c("U0", "U4", "U7Sh", "U14L"))
df.ridge_2.fun.dispersal.stem$Sample <- factor(df.ridge_2.fun.dispersal.stem$Sample, levels = c("U0", "U4", "U7Sh", "U14S"))
df.ridge_2.fun.dispersal.root$Sample <- factor(df.ridge_2.fun.dispersal.root$Sample, levels = c("U0", "U4", "U7R", "U14R"))


fun.list.otu_id <- fun.list[,c("OTU","OTU_id")]
df.ridge_2.fun.dispersal.stem <- merge(df.ridge_2.fun.dispersal.stem, fun.list.otu_id, by = "OTU")
df.ridge_2.fun.dispersal.root <- merge(df.ridge_2.fun.dispersal.root, fun.list.otu_id, by = "OTU")
df.ridge_2.fun.dispersal.leaf <- merge(df.ridge_2.fun.dispersal.leaf, fun.list.otu_id, by = "OTU")

df.ridge_2.fun.dispersal.stem$OTU_id <- factor(df.ridge_2.fun.dispersal.stem$OTU_id, levels = c("F8_Cladosporium", "F7_o_NA", "F2_f_Didymellaceae"))

df.ridge_2.fun.dispersal.root$OTU_id <- factor(df.ridge_2.fun.dispersal.root$OTU_id, levels = c("F3_Cladosporium", "F4_Alternaria", "F2_f_Didymellaceae"))

ggplot(df.ridge_2.fun.dispersal.leaf, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=OTU_id)) + 
  # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
  geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.8)+
  ylab("Relative abundance (%) \n") +
  theme(legend.text=element_text(size=12)) + 
  # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
  theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
  theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
  #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 0.5)
geom_vline(xintercept=18.5, color="slategray3", linetype='dashed',size=1)



### hulled seed
fun.clean.ss.H0 <- subset_samples(fun.clean.ss, Group == "H0")
fun.clean.ss.H0<-phyloseq::filter_taxa(fun.clean.ss.H0, function(x) sum(x) != 0, TRUE)

fun.clean.ss.H4 <- subset_samples(fun.clean.ss, Group == "H4")
fun.clean.ss.H4<-phyloseq::filter_taxa(fun.clean.ss.H4, function(x) sum(x) != 0, TRUE)

fun.clean.ss.H7Sh <- subset_samples(fun.clean.ss, Group == "H7Sh")
fun.clean.ss.H7Sh<-phyloseq::filter_taxa(fun.clean.ss.H7Sh, function(x) sum(x) != 0, TRUE)
fun.clean.ss.H7R <- subset_samples(fun.clean.ss, Group == "H7R")
fun.clean.ss.H7R<-phyloseq::filter_taxa(fun.clean.ss.H7R, function(x) sum(x) != 0, TRUE)

fun.clean.ss.H14L <- subset_samples(fun.clean.ss, Group == "H14L")
fun.clean.ss.H14L<-phyloseq::filter_taxa(fun.clean.ss.H14L, function(x) sum(x) != 0, TRUE)
fun.clean.ss.H14S <- subset_samples(fun.clean.ss, Group == "H14S")
fun.clean.ss.H14S<-phyloseq::filter_taxa(fun.clean.ss.H14S, function(x) sum(x) != 0, TRUE)
fun.clean.ss.H14R <- subset_samples(fun.clean.ss, Group == "H14R")
fun.clean.ss.H14R<-phyloseq::filter_taxa(fun.clean.ss.H14R, function(x) sum(x) != 0, TRUE)

Sh7.H<-intersect(taxa_names(fun.clean.ss.H0), taxa_names(fun.clean.ss.H7Sh)) #2
R7.H<-intersect(taxa_names(fun.clean.ss.H0), taxa_names(fun.clean.ss.H7R)) #8

fun.list.Sh7 <- subset(fun.list, OTU %in% Sh7.H)
fun.list.R7 <- subset(fun.list, OTU %in% R7.H)


Sh7.4.H<-intersect(taxa_names(fun.clean.ss.H4), taxa_names(fun.clean.ss.H7Sh)) #2
R7.4.H<-intersect(taxa_names(fun.clean.ss.H4), taxa_names(fun.clean.ss.H7R)) #8


zero_to_7.Sh.H<-intersect(Sh7.H, Sh7.4.H) #1
zero_to_7.R.H<-intersect(R7.H, R7.4.H) #7

fun.list.Sh7.H <- subset(fun.list, OTU %in% zero_to_7.Sh.H)
fun.list.R7.H <- subset(fun.list, OTU %in% zero_to_7.R.H)

### Seedlings 14-day-old
L14.H<-intersect(taxa_names(fun.clean.ss.H0), taxa_names(fun.clean.ss.H14L)) #5
S14.H<-intersect(taxa_names(fun.clean.ss.H0), taxa_names(fun.clean.ss.H14S)) #3
R14.H<-intersect(taxa_names(fun.clean.ss.H0), taxa_names(fun.clean.ss.H14R)) #7

L14.4.H<-intersect(taxa_names(fun.clean.ss.H4), taxa_names(fun.clean.ss.H14L)) #6
S14.4.H<-intersect(taxa_names(fun.clean.ss.H4), taxa_names(fun.clean.ss.H14S)) #9
R14.4.H<-intersect(taxa_names(fun.clean.ss.H4), taxa_names(fun.clean.ss.H14R)) #5


fun.list.L14.H <- subset(fun.list, OTU %in% L14.H)
fun.list.S14.H <- subset(fun.list, OTU %in% S14.H)
fun.list.R14.H <- subset(fun.list, OTU %in% R14.H)

sort(colSums(otu_table(fun.clean.ss)))

zero_to_14.L.H<-intersect(L14.H, L14.4.H) #4
zero_to_14.S.H <-intersect(S14.H, S14.4.H) #2
zero_to_14.R.H<-intersect(R14.H, R14.4.H) #3

fun.list.L14.H <- subset(fun.list, OTU %in% zero_to_14.L.H)
fun.list.S14.H <- subset(fun.list, OTU %in% zero_to_14.S.H)
fun.list.R14.H <- subset(fun.list, OTU %in% zero_to_14.R.H)

zero_to_14.R.fun.H <-zero_to_14.R.H
zero_to_14.S.fun.H <-zero_to_14.S.H
zero_to_14.L.fun.H <-zero_to_14.L.H


dispersal.L<-intersect(zero_to_7.Sh.H, zero_to_14.L.H)
dispersal.S<-intersect(zero_to_7.Sh.H, zero_to_14.S.H)
dispersal.R<-intersect(zero_to_7.R.H, zero_to_14.R.H)




### Distribution of abundance
library(ggridges)


get_df_ridge <- function(phy.seqs, compartment, dispersal){
  phy.seqs.tissue <- subset_samples(phy.seqs, Group %in% compartment)
  phy.seqs.m <-merge_samples(phy.seqs.tissue, "Group")
  df.otu <- phy.seqs.m %>% tax_glom(Genus)%>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance
  
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% dispersal)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  str(df.ridge)
  
  return(df.ridge)
}



compartment = c("H0", "H4", "H7Sh", "H14L")
df.ridge_2.fun.dispersal.leaf<- get_df_ridge(fun.clean.ss, dispersal = dispersal.L)

compartment = c("H0", "H4", "H7Sh", "H14S")
df.ridge_2.fun.dispersal.stem<- get_df_ridge(fun.clean.ss, dispersal = dispersal.S)

compartment = c("H0", "H4", "H7R", "H14R")
df.ridge_2.fun.dispersal.root<- get_df_ridge(fun.clean.ss, dispersal = dispersal.R)



df.ridge_2.fun.dispersal.leaf$Sample <- factor(df.ridge_2.fun.dispersal.leaf$Sample, levels = c("H0", "H4", "H7Sh", "H14L"))
df.ridge_2.fun.dispersal.stem$Sample <- factor(df.ridge_2.fun.dispersal.stem$Sample, levels = c("H0", "H4", "H7Sh", "H14S"))
df.ridge_2.fun.dispersal.root$Sample <- factor(df.ridge_2.fun.dispersal.root$Sample, levels = c("H0", "H4", "H7R", "H14R"))


fun.list.otu_id <- fun.list[,c("OTU","OTU_id")]
df.ridge_2.fun.dispersal.stem <- merge(df.ridge_2.fun.dispersal.stem, fun.list.otu_id, by = "OTU")
df.ridge_2.fun.dispersal.root <- merge(df.ridge_2.fun.dispersal.root, fun.list.otu_id, by = "OTU")
df.ridge_2.fun.dispersal.leaf <- merge(df.ridge_2.fun.dispersal.leaf, fun.list.otu_id, by = "OTU")

df.ridge_2.fun.dispersal.root$OTU_id <- factor(df.ridge_2.fun.dispersal.root$OTU_id, levels = c("F4_Alternaria", "F3_Cladosporium", "F1_Pyricularia"))

ggplot(df.ridge_2.fun.dispersal.root, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=OTU_id)) + 
  # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
  geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.8)+
  ylab("Relative abundance (%) \n") +
  theme(legend.text=element_text(size=12)) + 
  # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
  theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
  theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
  #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 0.5)
  geom_vline(xintercept=18.5, color="slategray3", linetype='dashed',size=1)

  
  
  
### Bacteria
  ### Dispersal of seed OTUs to compartments of seedlings
  ## OTUs in intact seeds
  
  bac.clean.ss.U0 <- subset_samples(bac.clean.ss, Group == "U0")
  bac.clean.ss.U0<-phyloseq::filter_taxa(bac.clean.ss.U0, function(x) sum(x) != 0, TRUE)
  
  bac.clean.ss.U4 <- subset_samples(bac.clean.ss, Group == "U4")
  bac.clean.ss.U4<-phyloseq::filter_taxa(bac.clean.ss.U4, function(x) sum(x) != 0, TRUE)
  
  
  intersect(taxa_names(bac.clean.ss.U0), taxa_names(bac.clean.ss.U4))
  
  bac.clean.ss.U7Sh <- subset_samples(bac.clean.ss, Group == "U7Sh")
  bac.clean.ss.U7Sh<-phyloseq::filter_taxa(bac.clean.ss.U7Sh, function(x) sum(x) != 0, TRUE)
  bac.clean.ss.U7R <- subset_samples(bac.clean.ss, Group == "U7R")
  bac.clean.ss.U7R<-phyloseq::filter_taxa(bac.clean.ss.U7R, function(x) sum(x) != 0, TRUE)
  
  bac.clean.ss.U14L <- subset_samples(bac.clean.ss, Group == "U14L")
  bac.clean.ss.U14L<-phyloseq::filter_taxa(bac.clean.ss.U14L, function(x) sum(x) != 0, TRUE)
  bac.clean.ss.U14S <- subset_samples(bac.clean.ss, Group == "U14S")
  bac.clean.ss.U14S<-phyloseq::filter_taxa(bac.clean.ss.U14S, function(x) sum(x) != 0, TRUE)
  bac.clean.ss.U14R <- subset_samples(bac.clean.ss, Group == "U14R")
  bac.clean.ss.U14R<-phyloseq::filter_taxa(bac.clean.ss.U14R, function(x) sum(x) != 0, TRUE)
  
  Sh7<-intersect(taxa_names(bac.clean.ss.U0), taxa_names(bac.clean.ss.U7Sh)) #9
  R7<-intersect(taxa_names(bac.clean.ss.U0), taxa_names(bac.clean.ss.U7R)) #3
  
  Sh7.4<-intersect(taxa_names(bac.clean.ss.U4), taxa_names(bac.clean.ss.U7Sh)) #6
  R7.4<-intersect(taxa_names(bac.clean.ss.U4), taxa_names(bac.clean.ss.U7R)) #5
  
  
  zero_to_7.Sh<-intersect(Sh7, Sh7.4) #5
  zero_to_7.R<-intersect(R7, R7.4) #3
  
  bac.list.Sh7 <- subset(bac.list, OTU %in% zero_to_7.Sh)
  bac.list.R7 <- subset(bac.list, OTU %in% zero_to_7.R)
  
  ### Seedlings 14-day-old
  L14<-intersect(taxa_names(bac.clean.ss.U0), taxa_names(bac.clean.ss.U14L)) #9
  S14<-intersect(taxa_names(bac.clean.ss.U0), taxa_names(bac.clean.ss.U14S)) #13
  R14<-intersect(taxa_names(bac.clean.ss.U0), taxa_names(bac.clean.ss.U14R)) #3
  
  L14.4<-intersect(taxa_names(bac.clean.ss.U4), taxa_names(bac.clean.ss.U14L)) #6
  S14.4<-intersect(taxa_names(bac.clean.ss.U4), taxa_names(bac.clean.ss.U14S)) #8
  R14.4<-intersect(taxa_names(bac.clean.ss.U4), taxa_names(bac.clean.ss.U14R)) #5
  
  
  bac.list.L14 <- subset(bac.list, OTU %in% L14)
  bac.list.S14 <- subset(bac.list, OTU %in% S14)
  bac.list.R14 <- subset(bac.list, OTU %in% R14)
  
  sort(colSums(otu_table(bac.clean.ss)))
  
  zero_to_14.L<-intersect(L14, L14.4) #5
  zero_to_14.S <-intersect(S14, S14.4) #7
  zero_to_14.R<-intersect(R14, R14.4) #3
  
  bac.list.L14 <- subset(bac.list, OTU %in% zero_to_14.L)
  bac.list.S14 <- subset(bac.list, OTU %in% zero_to_14.S)
  bac.list.R14 <- subset(bac.list, OTU %in% zero_to_14.R)
  
  
  L7.14<-intersect(taxa_names(bac.clean.ss.U14L), taxa_names(bac.clean.ss.U7Sh)) #6
  S7.14<-intersect(taxa_names(bac.clean.ss.U14S), taxa_names(bac.clean.ss.U7Sh)) #5
  R7.14<-intersect(taxa_names(bac.clean.ss.U14R), taxa_names(bac.clean.ss.U7R)) #5
  
  bac.list.L7.14 <- subset(bac.list, OTU %in% L7.14)
  bac.list.S7.14 <- subset(bac.list, OTU %in%  S7.14)
  bac.list.R7.14 <- subset(bac.list, OTU %in% R7.14)
  
  dispersal.L<-intersect(zero_to_7.Sh, zero_to_14.L)
  dispersal.S<-intersect(zero_to_7.Sh, zero_to_14.S)
  dispersal.R<-intersect(zero_to_7.R, zero_to_14.R)
  
  zero_to_14.R.bac <-zero_to_14.R
  zero_to_14.S.bac <-zero_to_14.S
  zero_to_14.L.bac <-zero_to_14.L
  
  compartment = c("U0", "U4", "U7Sh", "U14L")
  dispersal = L14.4
  df.ridge_2.bac.dispersal.leaf<- get_df_ridge(bac.clean.ss, dispersal = L14.4)
  
  compartment = c("U0", "U4", "U7Sh", "U14S")
  df.ridge_2.bac.dispersal.stem<- get_df_ridge(bac.clean.ss, dispersal = S14.4)
  
  compartment = c("U0", "U4", "U7R", "U14R")
  df.ridge_2.bac.dispersal.root<- get_df_ridge(bac.clean.ss, dispersal = R14.4)
  
  
  
  df.ridge_2.bac.dispersal.leaf$Sample <- factor(df.ridge_2.bac.dispersal.leaf$Sample, levels = c("U0", "U4", "U7Sh", "U14L"))
  df.ridge_2.bac.dispersal.stem$Sample <- factor(df.ridge_2.bac.dispersal.stem$Sample, levels = c("U0", "U4", "U7Sh", "U14S"))
  df.ridge_2.bac.dispersal.root$Sample <- factor(df.ridge_2.bac.dispersal.root$Sample, levels = c("U0", "U4", "U7R", "U14R"))
  
  
  bac.list.otu_id <- bac.list[,c("OTU","OTU_id")]
  df.ridge_2.bac.dispersal.stem <- merge(df.ridge_2.bac.dispersal.stem, bac.list.otu_id, by = "OTU")
  df.ridge_2.bac.dispersal.root <- merge(df.ridge_2.bac.dispersal.root, bac.list.otu_id, by = "OTU")
  df.ridge_2.bac.dispersal.leaf <- merge(df.ridge_2.bac.dispersal.leaf, bac.list.otu_id, by = "OTU")
  
  df.ridge_2.bac.dispersal.stem$OTU_id <- factor(df.ridge_2.bac.dispersal.stem$OTU_id, levels = c("F8_Cladosporium", "F7_o_NA", "F2_f_Didymellaceae"))
  
  df.ridge_2.bac.dispersal.root$OTU_id <- factor(df.ridge_2.bac.dispersal.root$OTU_id, levels = c("F3_Cladosporium", "F4_Alternaria", "F2_f_Didymellaceae"))
  
  
  library(ggridges)
  ggplot(df.ridge_2.bac.dispersal.leaf, aes(x=Sample, y=OTU, height=RelAbundance, group = OTU_id, fill=OTU_id)) + 
    # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
    geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.8)+
    ylab("Relative abundance (%) \n") +
    theme(legend.text=element_text(size=12)) + 
    # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    theme(aspect.ratio = 0.5)+
  geom_vline(xintercept=18.5, color="slategray3", linetype='dashed',size=1)
  
  
  
  ### hulled seed
  bac.clean.ss.H0 <- subset_samples(bac.clean.ss, Group == "H0")
  bac.clean.ss.H0<-phyloseq::filter_taxa(bac.clean.ss.H0, function(x) sum(x) != 0, TRUE)
  
  bac.clean.ss.H4 <- subset_samples(bac.clean.ss, Group == "H4")
  bac.clean.ss.H4<-phyloseq::filter_taxa(bac.clean.ss.H4, function(x) sum(x) != 0, TRUE)
  
  bac.clean.ss.H7Sh <- subset_samples(bac.clean.ss, Group == "H7Sh")
  bac.clean.ss.H7Sh<-phyloseq::filter_taxa(bac.clean.ss.H7Sh, function(x) sum(x) != 0, TRUE)
  bac.clean.ss.H7R <- subset_samples(bac.clean.ss, Group == "H7R")
  bac.clean.ss.H7R<-phyloseq::filter_taxa(bac.clean.ss.H7R, function(x) sum(x) != 0, TRUE)
  
  bac.clean.ss.H14L <- subset_samples(bac.clean.ss, Group == "H14L")
  bac.clean.ss.H14L<-phyloseq::filter_taxa(bac.clean.ss.H14L, function(x) sum(x) != 0, TRUE)
  bac.clean.ss.H14S <- subset_samples(bac.clean.ss, Group == "H14S")
  bac.clean.ss.H14S<-phyloseq::filter_taxa(bac.clean.ss.H14S, function(x) sum(x) != 0, TRUE)
  bac.clean.ss.H14R <- subset_samples(bac.clean.ss, Group == "H14R")
  bac.clean.ss.H14R<-phyloseq::filter_taxa(bac.clean.ss.H14R, function(x) sum(x) != 0, TRUE)
  
  Sh7.H<-intersect(taxa_names(bac.clean.ss.H0), taxa_names(bac.clean.ss.H7Sh)) #2
  R7.H<-intersect(taxa_names(bac.clean.ss.H0), taxa_names(bac.clean.ss.H7R)) #8
  
  bac.list.Sh7 <- subset(bac.list, OTU %in% Sh7.H)
  bac.list.R7 <- subset(bac.list, OTU %in% R7.H)
  
  
  Sh7.4.H<-intersect(taxa_names(bac.clean.ss.H4), taxa_names(bac.clean.ss.H7Sh)) #2
  R7.4.H<-intersect(taxa_names(bac.clean.ss.H4), taxa_names(bac.clean.ss.H7R)) #8
  
  L7.14.H<-intersect(taxa_names(bac.clean.ss.H14L), taxa_names(bac.clean.ss.H7Sh)) #2
  S7.14.H<-intersect(taxa_names(bac.clean.ss.H14S), taxa_names(bac.clean.ss.H7Sh)) #8
  
  L7.14.H<-intersect(taxa_names(bac.clean.ss.H14L), taxa_names(bac.clean.ss.H7Sh)) #2
  S7.14.H<-intersect(taxa_names(bac.clean.ss.H14S), taxa_names(bac.clean.ss.H7Sh)) #8
  R7.14.H<-intersect(taxa_names(bac.clean.ss.H14R), taxa_names(bac.clean.ss.H7R)) #8
  
  bac.list.L7.14H <- subset(bac.list, OTU %in%  L7.14.H)
  bac.list.S7.14H <- subset(bac.list, OTU %in% S7.14.H)
  bac.list.R7.14H <- subset(bac.list, OTU %in% R7.14.H)
  
  
  
  zero_to_7.Sh.H<-intersect(Sh7.H, Sh7.4.H) #1
  zero_to_7.R.H<-intersect(R7.H, R7.4.H) #7
  
  bac.list.Sh7.H <- subset(bac.list, OTU %in% zero_to_7.Sh.H)
  bac.list.R7.H <- subset(bac.list, OTU %in% zero_to_7.R.H)
  
  ### Seedlings 14-day-old
  L14.H<-intersect(taxa_names(bac.clean.ss.H0), taxa_names(bac.clean.ss.H14L)) #5
  S14.H<-intersect(taxa_names(bac.clean.ss.H0), taxa_names(bac.clean.ss.H14S)) #3
  R14.H<-intersect(taxa_names(bac.clean.ss.H0), taxa_names(bac.clean.ss.H14R)) #7
  
  L14.4.H<-intersect(taxa_names(bac.clean.ss.H4), taxa_names(bac.clean.ss.H14L)) #6
  S14.4.H<-intersect(taxa_names(bac.clean.ss.H4), taxa_names(bac.clean.ss.H14S)) #9
  R14.4.H<-intersect(taxa_names(bac.clean.ss.H4), taxa_names(bac.clean.ss.H14R)) #5
  
  
  bac.list.L14.H <- subset(bac.list, OTU %in% L14.H)
  bac.list.S14.H <- subset(bac.list, OTU %in% S14.H)
  bac.list.R14.H <- subset(bac.list, OTU %in% R14.H)
  
  zero_to_14.L.H<-intersect(L14.H, L14.4.H) #4
  zero_to_14.S.H <-intersect(S14.H, S14.4.H) #2
  zero_to_14.R.H<-intersect(R14.H, R14.4.H) #3
  
  bac.list.L14.H <- subset(bac.list, OTU %in% zero_to_14.L.H)
  bac.list.S14.H <- subset(bac.list, OTU %in% zero_to_14.S.H)
  bac.list.R14.H <- subset(bac.list, OTU %in% zero_to_14.R.H)
  

  
  zero_to_14.L.H.bac <- zero_to_14.L.H
  zero_to_14.S.H.bac <- zero_to_14.S.H
  zero_to_14.R.H.bac <- zero_to_14.R.H
  
  dispersal.L.H<-intersect(zero_to_7.Sh, zero_to_14.L.H)
  dispersal.S.H<-intersect(zero_to_7.Sh, zero_to_14.S.H)
  dispersal.R.H<-intersect(zero_to_7.R, zero_to_14.R.H)
  
  
  ### Distribution of abundance
  library(ggridges)
  
  
  get_df_ridge <- function(bac.seqs, compartment, dispersal){
    bac.seqs.tissue <- subset_samples(bac.seqs, Group %in% compartment)
    bac.seqs.m <-merge_samples(bac.seqs.tissue, "Group")
    df.otu <- bac.seqs.m %>% psmelt()
    head(df.otu)
    # we need to group by samples
    df.otu.rel <- df.otu %>%  
      group_by(Sample) %>%                         # Filter out at absolute read of 20       
      mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance
    
    
    df.selected.rel <- df.otu.rel %>% filter(OTU %in% dispersal)
    df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
    str(df.ridge)
    
    return(df.ridge)
  }
  
  
  
  compartment = c("H0", "H4", "H7Sh", "H14L")
  df.ridge_2.bac.dispersal.leaf<- get_df_ridge(bac.clean.ss, dispersal = dispersal.L)
  
  compartment = c("H0", "H4", "H7Sh", "H14S")
  df.ridge_2.bac.dispersal.stem<- get_df_ridge(bac.clean.ss, dispersal = dispersal.S)
  
  compartment = c("H0", "H4", "H7R", "H14R")
  df.ridge_2.bac.dispersal.root<- get_df_ridge(bac.clean.ss, dispersal = dispersal.R)
  
  
  
  df.ridge_2.bac.dispersal.leaf$Sample <- factor(df.ridge_2.bac.dispersal.leaf$Sample, levels = c("H0", "H4", "H7Sh", "H14L"))
  df.ridge_2.bac.dispersal.stem$Sample <- factor(df.ridge_2.bac.dispersal.stem$Sample, levels = c("H0", "H4", "H7Sh", "H14S"))
  df.ridge_2.bac.dispersal.root$Sample <- factor(df.ridge_2.bac.dispersal.root$Sample, levels = c("H0", "H4", "H7R", "H14R"))
  
  
  bac.list.otu_id <- bac.list[,c("OTU","OTU_id")]
  df.ridge_2.bac.dispersal.stem <- merge(df.ridge_2.bac.dispersal.stem, bac.list.otu_id, by = "OTU")
  df.ridge_2.bac.dispersal.root <- merge(df.ridge_2.bac.dispersal.root, bac.list.otu_id, by = "OTU")
  df.ridge_2.bac.dispersal.leaf <- merge(df.ridge_2.bac.dispersal.leaf, bac.list.otu_id, by = "OTU")
  
  df.ridge_2.bac.dispersal.root$OTU_id <- factor(df.ridge_2.bac.dispersal.root$OTU_id, levels = c("F4_Alternaria", "F3_Cladosporium", "F1_Pyricularia"))
  
  ggplot(df.ridge_2.bac.dispersal.stem, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=OTU_id)) + 
    # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
    geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.8)+
    ylab("Relative abundance (%) \n") +
    theme(legend.text=element_text(size=12)) + 
    # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    theme(aspect.ratio = 0.5)
  
  ##### Source
  bac.list.L14.H$Compartment <- "Hulled L14"
  bac.list.S14.H$Compartment <- "Hulled S14"
  bac.list.R14.H$Compartment <- "Hulled R14"
  bac.list.seed.origin.H<- rbind(bac.list.L14.H,bac.list.S14.H, bac.list.R14.H)
  
  fun.list.L14.H$Compartment <- "Hulled L14"
  fun.list.S14.H$Compartment <- "Hulled S14"
  fun.list.R14.H$Compartment <- "Hulled R14"
  fun.list.seed.origin.H<- rbind(fun.list.L14.H,fun.list.S14.H, fun.list.R14.H)
  
  
  
  bac.list.L14$Compartment <- "Intact L14"
  bac.list.S14$Compartment <- "Intact S14"
  bac.list.R14$Compartment <- "Intact R14"
  bac.list.seed.origin<- rbind(bac.list.L14,bac.list.S14, bac.list.R14)
  
  fun.list.L14$Compartment <- "Intact L14"
  fun.list.S14$Compartment <- "Intact S14"
  fun.list.R14$Compartment <- "Intact R14"
  fun.list.seed.origin<- rbind(fun.list.L14,fun.list.S14, fun.list.R14)
  
  write.csv(bac.list.seed.origin.H, "seed origin-hulled-bacteria.csv")
  write.csv(fun.list.seed.origin.H, "seed origin-hulled-fungi.csv")
  
  write.csv(bac.list.seed.origin, "seed origin-intact-bacteria.csv")
  write.csv(fun.list.seed.origin, "seed origin-intact-fungi.csv")