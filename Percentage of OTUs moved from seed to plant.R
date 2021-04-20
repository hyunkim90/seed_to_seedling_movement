### Distribution of seed-borne OTUs in above and belowground tissues
### Using samples collected at the tillering (48), heading (90), and harvest (141)

bac.clean.ss.18.sub<-subset_samples(bac.clean.ss.18, Days %in% c("0",'48','90','141'))
fun.clean.ss.18.sub<-subset_samples(fun.clean.ss.18, Days %in% c("0",'48','90','141'))


b.meta.sub<-sample_data(bac.clean.ss.18.sub)
f.meta.sub<-sample_data(fun.clean.ss.18.sub)

df.b.meta.sub <- data.frame(b.meta.sub)
df.f.meta.sub <- data.frame(f.meta.sub)

df.b.meta.sub$Compart_age <- paste(df.b.meta.sub$Compartment, df.b.meta.sub$Days, sep = "_")
head(df.b.meta.sub)


df.f.meta.sub$Compart_age <- paste(df.f.meta.sub$Compartment, df.f.meta.sub$Days, sep = "_")
head(df.f.meta.sub)


sample_data(bac.clean.ss.18.sub) <- sample_data(df.b.meta.sub)
sample_data(fun.clean.ss.18.sub) <- sample_data(df.f.meta.sub)


bac.clean.ss.18.sub.merged<-merge_samples(bac.clean.ss.18.sub, "Compart_age")
sample_data(bac.clean.ss.18.sub.merged)
t(otu_table(bac.clean.ss.18.sub.merged))

fun.clean.ss.18.sub.merged<-merge_samples(fun.clean.ss.18.sub, "Compart_age")
sample_data(fun.clean.ss.18.sub.merged)
t(otu_table(fun.clean.ss.18.sub.merged))

###remove OTUs with 0 read
bac.clean.ss.18.sub.merged <- phyloseq::filter_taxa(bac.clean.ss.18.sub.merged, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.sub.merged <- phyloseq::filter_taxa(fun.clean.ss.18.sub.merged, function(x) sum(x) != 0, TRUE)


### Percentage of OTUs moved from seed to plant
## Bacteria
G0.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Seed_0')]
head(G0.bac)
G0.bac<-data.frame(G0.bac) 
G0.bac.over0 <-subset(G0.bac, Seed_0 > 0)

G0.bac.list<-rownames(G0.bac.over0)

### G141
G141.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Seed_141')]
head(G141.bac)
G141.bac<-data.frame(G141.bac) 
G141.bac.over141 <-subset(G141.bac, Seed_141 > 0)

G141.bac.list<-rownames(G141.bac.over141)


G0_G141_common.bac<-intersect(G0.bac.list, G141.bac.list) #36


### L48
L48.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Leaf_48')]
head(L48.bac)
L48.bac<-data.frame(L48.bac) 
L48.bac.over0 <-subset(L48.bac, Leaf_48 > 0)

L48.bac.list<-rownames(L48.bac.over0)


G0_L48_common.bac<-intersect(G0.bac.list, L48.bac.list) #11


### L90
L90.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Leaf_90')]
head(L90.bac)
L90.bac<-data.frame(L90.bac) 
L90.bac.over0 <-subset(L90.bac, Leaf_90 > 0)

L90.bac.list<-rownames(L90.bac.over0)


G0_L90_common.bac<-intersect(G0.bac.list, L90.bac.list) #22

### L141
L141.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Leaf_141')]
head(L141.bac)
L141.bac<-data.frame(L141.bac) 
L141.bac.over0 <-subset(L141.bac, Leaf_141 > 0)

L141.bac.list<-rownames(L141.bac.over0)


G0_L141_common.bac<-intersect(G0.bac.list, L141.bac.list) #49


### S48
S48.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Stem_48')]
head(S48.bac)
S48.bac<-data.frame(S48.bac) 
S48.bac.over0 <-subset(S48.bac, Stem_48 > 0)

S48.bac.list<-rownames(S48.bac.over0)


G0_S48_common.bac<-intersect(G0.bac.list, S48.bac.list) #25


### S90
S90.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Stem_90')]
head(S90.bac)
S90.bac<-data.frame(S90.bac) 
S90.bac.over0 <-subset(S90.bac, Stem_90 > 0)

S90.bac.list<-rownames(S90.bac.over0)


G0_S90_common.bac<-intersect(G0.bac.list, S90.bac.list) #51

### S141
S141.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Stem_141')]
head(S141.bac)
S141.bac<-data.frame(S141.bac) 
S141.bac.over0 <-subset(S141.bac, Stem_141 > 0)

S141.bac.list<-rownames(S141.bac.over0)


G0_S141_common.bac<-intersect(G0.bac.list, S141.bac.list) #65



### R48
R48.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Root_48')]
head(R48.bac)
R48.bac<-data.frame(R48.bac) 
R48.bac.over0 <-subset(R48.bac, Root_48 > 0)

R48.bac.list<-rownames(R48.bac.over0)


G0_R48_common.bac<-intersect(G0.bac.list, R48.bac.list) #7


### R90
R90.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Root_90')]
head(R90.bac)
R90.bac<-data.frame(R90.bac) 
R90.bac.over0 <-subset(R90.bac, Root_90 > 0)

R90.bac.list<-rownames(R90.bac.over0)


G0_R90_common.bac<-intersect(G0.bac.list, R90.bac.list) #6

### R141
R141.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Root_141')]
head(R141.bac)
R141.bac<-data.frame(R141.bac) 
R141.bac.over0 <-subset(R141.bac, Root_141 > 0)

R141.bac.list<-rownames(R141.bac.over0)


G0_R141_common.bac<-intersect(G0.bac.list, R141.bac.list) #27


### RS48
RS48.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Rhizosphere_48')]
head(RS48.bac)
RS48.bac<-data.frame(RS48.bac) 
RS48.bac.over0 <-subset(RS48.bac, Rhizosphere_48 > 0)

RS48.bac.list<-rownames(RS48.bac.over0)


G0_RS48_common.bac<-intersect(G0.bac.list, RS48.bac.list) #1


### RS90
RS90.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Rhizosphere_90')]
head(RS90.bac)
RS90.bac<-data.frame(RS90.bac) 
RS90.bac.over0 <-subset(RS90.bac, Rhizosphere_90 > 0)

RS90.bac.list<-rownames(RS90.bac.over0)


G0_RS90_common.bac<-intersect(G0.bac.list, RS90.bac.list) #1

### RS141
RS141.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Rhizosphere_141')]
head(RS141.bac)
RS141.bac<-data.frame(RS141.bac) 
RS141.bac.over0 <-subset(RS141.bac, Rhizosphere_141 > 0)

RS141.bac.list<-rownames(RS141.bac.over0)


G0_RS141_common.bac<-intersect(G0.bac.list, RS141.bac.list) #0


### BS48
BS48.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Bulk_soil_48')]
head(BS48.bac)
BS48.bac<-data.frame(BS48.bac) 
BS48.bac.over0 <-subset(BS48.bac, Bulk_soil_48 > 0)

BS48.bac.list<-rownames(BS48.bac.over0)


G0_BS48_common.bac<-intersect(G0.bac.list, BS48.bac.list) #0


### BS90
BS90.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Bulk_soil_90')]
head(BS90.bac)
BS90.bac<-data.frame(BS90.bac) 
BS90.bac.over0 <-subset(BS90.bac, Bulk_soil_90 > 0)

BS90.bac.list<-rownames(BS90.bac.over0)


G0_BS90_common.bac<-intersect(G0.bac.list, BS90.bac.list) #0

### BS141
BS141.bac<-t(otu_table(bac.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(bac.clean.ss.18.sub.merged))) == 'Bulk_soil_141')]
head(BS141.bac)
BS141.bac<-data.frame(BS141.bac) 
BS141.bac.over0 <-subset(BS141.bac, Bulk_soil_141 > 0)

BS141.bac.list<-rownames(BS141.bac.over0)


G0_BS141_common.bac<-intersect(G0.bac.list, BS141.bac.list) #1

G0_G141_common.bac.list<-subset(bac.list, OTU %in% G0_G141_common.bac)
G0_G141_common.bac.list$Group <- "G0_G141"

G0_L48_common.bac.list<-subset(bac.list, OTU %in% G0_L48_common.bac)
G0_L48_common.bac.list$Group <- "G0_L48"

G0_L90_common.bac.list<-subset(bac.list, OTU %in% G0_L90_common.bac)
G0_L90_common.bac.list$Group <- "G0_L90"

G0_L141_common.bac.list<-subset(bac.list, OTU %in% G0_L141_common.bac)
G0_L141_common.bac.list$Group <- "G0_L141"

G0_S48_common.bac.list<-subset(bac.list, OTU %in% G0_S48_common.bac)
G0_S48_common.bac.list$Group <- "G0_S48"

G0_S90_common.bac.list<-subset(bac.list, OTU %in% G0_S90_common.bac)
G0_S90_common.bac.list$Group <- "G0_S90"

G0_S141_common.bac.list<-subset(bac.list, OTU %in% G0_S141_common.bac)
G0_S141_common.bac.list$Group <- "G0_S141"

G0_R48_common.bac.list<-subset(bac.list, OTU %in% G0_R48_common.bac)
G0_R48_common.bac.list$Group <- "G0_R48"

G0_R90_common.bac.list<-subset(bac.list, OTU %in% G0_R90_common.bac)
G0_R90_common.bac.list$Group <- "G0_R90"

G0_R141_common.bac.list<-subset(bac.list, OTU %in% G0_R141_common.bac)
G0_R141_common.bac.list$Group <- "G0_R141"

G0_RS48_common.bac.list<-subset(bac.list, OTU %in% G0_RS48_common.bac)
G0_RS48_common.bac.list$Group <- "G0_RS48"

G0_RS90_common.bac.list<-subset(bac.list, OTU %in% G0_RS90_common.bac)
G0_RS90_common.bac.list$Group <- "G0_RS90"

G0_RS141_common.bac.list<-subset(bac.list, OTU %in% G0_RS141_common.bac)
G0_RS141_common.bac.list$Group <- "G0_RS141"

G0_BS48_common.bac.list<-subset(bac.list, OTU %in% G0_BS48_common.bac)
G0_BS48_common.bac.list$Group <- "G0_BS48"

G0_BS90_common.bac.list<-subset(bac.list, OTU %in% G0_BS90_common.bac)
G0_BS90_common.bac.list$Group <- "G0_BS90"

G0_BS141_common.bac.list<-subset(bac.list, OTU %in% G0_BS141_common.bac)
G0_BS141_common.bac.list$Group <- "G0_BS141"


moveOTU.list.bac <- rbind(G0_G141_common.bac.list, G0_L48_common.bac.list, G0_L90_common.bac.list, G0_L141_common.bac.list,
                          G0_S48_common.bac.list, G0_S90_common.bac.list, G0_S141_common.bac.list, 
                          G0_R48_common.bac.list, G0_R90_common.bac.list, G0_R141_common.bac.list,
                          G0_RS48_common.bac.list, G0_RS90_common.bac.list, G0_RS141_common.bac.list,
                          G0_BS48_common.bac.list, G0_BS90_common.bac.list, G0_BS141_common.bac.list)

write.csv(moveOTU.list.bac, "Seed to plant OTUs_bacteria.csv")



##Fungi
G0.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Grain_0')]
head(G0.fun)
G0.fun<-data.frame(G0.fun) 
G0.fun.over0 <-subset(G0.fun, Grain_0 > 0)

G0.fun.list<-rownames(G0.fun.over0)

### G141
G141.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Grain_141')]
head(G141.fun)
G141.fun<-data.frame(G141.fun) 
G141.fun.over141 <-subset(G141.fun, Grain_141 > 0)

G141.fun.list<-rownames(G141.fun.over141)


G0_G141_common.fun<-intersect(G0.fun.list, G141.fun.list) #63


### L48
L48.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Leaf_48')]
head(L48.fun)
L48.fun<-data.frame(L48.fun) 
L48.fun.over0 <-subset(L48.fun, Leaf_48 > 0)

L48.fun.list<-rownames(L48.fun.over0)


G0_L48_common.fun<-intersect(G0.fun.list, L48.fun.list) #11


### L90
L90.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Leaf_90')]
head(L90.fun)
L90.fun<-data.frame(L90.fun) 
L90.fun.over0 <-subset(L90.fun, Leaf_90 > 0)

L90.fun.list<-rownames(L90.fun.over0)


G0_L90_common.fun<-intersect(G0.fun.list, L90.fun.list) #22

### L141
L141.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Leaf_141')]
head(L141.fun)
L141.fun<-data.frame(L141.fun) 
L141.fun.over0 <-subset(L141.fun, Leaf_141 > 0)

L141.fun.list<-rownames(L141.fun.over0)


G0_L141_common.fun<-intersect(G0.fun.list, L141.fun.list) #49


### S48
S48.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Stem_48')]
head(S48.fun)
S48.fun<-data.frame(S48.fun) 
S48.fun.over0 <-subset(S48.fun, Stem_48 > 0)

S48.fun.list<-rownames(S48.fun.over0)


G0_S48_common.fun<-intersect(G0.fun.list, S48.fun.list) #25


### S90
S90.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Stem_90')]
head(S90.fun)
S90.fun<-data.frame(S90.fun) 
S90.fun.over0 <-subset(S90.fun, Stem_90 > 0)

S90.fun.list<-rownames(S90.fun.over0)


G0_S90_common.fun<-intersect(G0.fun.list, S90.fun.list) #51

### S141
S141.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Stem_141')]
head(S141.fun)
S141.fun<-data.frame(S141.fun) 
S141.fun.over0 <-subset(S141.fun, Stem_141 > 0)

S141.fun.list<-rownames(S141.fun.over0)


G0_S141_common.fun<-intersect(G0.fun.list, S141.fun.list) #65



### R48
R48.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Root_48')]
head(R48.fun)
R48.fun<-data.frame(R48.fun) 
R48.fun.over0 <-subset(R48.fun, Root_48 > 0)

R48.fun.list<-rownames(R48.fun.over0)


G0_R48_common.fun<-intersect(G0.fun.list, R48.fun.list) #7


### R90
R90.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Root_90')]
head(R90.fun)
R90.fun<-data.frame(R90.fun) 
R90.fun.over0 <-subset(R90.fun, Root_90 > 0)

R90.fun.list<-rownames(R90.fun.over0)


G0_R90_common.fun<-intersect(G0.fun.list, R90.fun.list) #6

### R141
R141.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Root_141')]
head(R141.fun)
R141.fun<-data.frame(R141.fun) 
R141.fun.over0 <-subset(R141.fun, Root_141 > 0)

R141.fun.list<-rownames(R141.fun.over0)


G0_R141_common.fun<-intersect(G0.fun.list, R141.fun.list) #27


### RS48
RS48.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Rhizosphere_48')]
head(RS48.fun)
RS48.fun<-data.frame(RS48.fun) 
RS48.fun.over0 <-subset(RS48.fun, Rhizosphere_48 > 0)

RS48.fun.list<-rownames(RS48.fun.over0)


G0_RS48_common.fun<-intersect(G0.fun.list, RS48.fun.list) #1


### RS90
RS90.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Rhizosphere_90')]
head(RS90.fun)
RS90.fun<-data.frame(RS90.fun) 
RS90.fun.over0 <-subset(RS90.fun, Rhizosphere_90 > 0)

RS90.fun.list<-rownames(RS90.fun.over0)


G0_RS90_common.fun<-intersect(G0.fun.list, RS90.fun.list) #1

### RS141
RS141.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Rhizosphere_141')]
head(RS141.fun)
RS141.fun<-data.frame(RS141.fun) 
RS141.fun.over0 <-subset(RS141.fun, Rhizosphere_141 > 0)

RS141.fun.list<-rownames(RS141.fun.over0)


G0_RS141_common.fun<-intersect(G0.fun.list, RS141.fun.list) #0


### BS48
BS48.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Bulk_soil_48')]
head(BS48.fun)
BS48.fun<-data.frame(BS48.fun) 
BS48.fun.over0 <-subset(BS48.fun, Bulk_soil_48 > 0)

BS48.fun.list<-rownames(BS48.fun.over0)


G0_BS48_common.fun<-intersect(G0.fun.list, BS48.fun.list) #0


### BS90
BS90.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Bulk_soil_90')]
head(BS90.fun)
BS90.fun<-data.frame(BS90.fun) 
BS90.fun.over0 <-subset(BS90.fun, Bulk_soil_90 > 0)

BS90.fun.list<-rownames(BS90.fun.over0)


G0_BS90_common.fun<-intersect(G0.fun.list, BS90.fun.list) #0

### BS141
BS141.fun<-t(otu_table(fun.clean.ss.18.sub.merged))[,which(colnames(t(otu_table(fun.clean.ss.18.sub.merged))) == 'Bulk_soil_141')]
head(BS141.fun)
BS141.fun<-data.frame(BS141.fun) 
BS141.fun.over0 <-subset(BS141.fun, Bulk_soil_141 > 0)

BS141.fun.list<-rownames(BS141.fun.over0)


G0_BS141_common.fun<-intersect(G0.fun.list, BS141.fun.list) #1

G0_G141_common.fun.list<-subset(fun.list, OTU %in% G0_G141_common.fun)
G0_G141_common.fun.list$Group <- "G0_G141"

G0_L48_common.fun.list<-subset(fun.list, OTU %in% G0_L48_common.fun)
G0_L48_common.fun.list$Group <- "G0_L48"

G0_L90_common.fun.list<-subset(fun.list, OTU %in% G0_L90_common.fun)
G0_L90_common.fun.list$Group <- "G0_L90"

G0_L141_common.fun.list<-subset(fun.list, OTU %in% G0_L141_common.fun)
G0_L141_common.fun.list$Group <- "G0_L141"

G0_S48_common.fun.list<-subset(fun.list, OTU %in% G0_S48_common.fun)
G0_S48_common.fun.list$Group <- "G0_S48"

G0_S90_common.fun.list<-subset(fun.list, OTU %in% G0_S90_common.fun)
G0_S90_common.fun.list$Group <- "G0_S90"

G0_S141_common.fun.list<-subset(fun.list, OTU %in% G0_S141_common.fun)
G0_S141_common.fun.list$Group <- "G0_S141"

G0_R48_common.fun.list<-subset(fun.list, OTU %in% G0_R48_common.fun)
G0_R48_common.fun.list$Group <- "G0_R48"

G0_R90_common.fun.list<-subset(fun.list, OTU %in% G0_R90_common.fun)
G0_R90_common.fun.list$Group <- "G0_R90"

G0_R141_common.fun.list<-subset(fun.list, OTU %in% G0_R141_common.fun)
G0_R141_common.fun.list$Group <- "G0_R141"

G0_RS48_common.fun.list<-subset(fun.list, OTU %in% G0_RS48_common.fun)
G0_RS48_common.fun.list$Group <- "G0_RS48"

G0_RS90_common.fun.list<-subset(fun.list, OTU %in% G0_RS90_common.fun)
G0_RS90_common.fun.list$Group <- "G0_RS90"

G0_RS141_common.fun.list<-subset(fun.list, OTU %in% G0_RS141_common.fun)
G0_RS141_common.fun.list$Group <- "G0_RS141"

G0_BS48_common.fun.list<-subset(fun.list, OTU %in% G0_BS48_common.fun)
G0_BS48_common.fun.list$Group <- "G0_BS48"

G0_BS90_common.fun.list<-subset(fun.list, OTU %in% G0_BS90_common.fun)
G0_BS90_common.fun.list$Group <- "G0_BS90"

G0_BS141_common.fun.list<-subset(fun.list, OTU %in% G0_BS141_common.fun)
G0_BS141_common.fun.list$Group <- "G0_BS141"


moveOTU.list.fun <- rbind(G0_G141_common.fun.list, G0_L48_common.fun.list, G0_L90_common.fun.list, G0_L141_common.fun.list,
                          G0_S48_common.fun.list, G0_S90_common.fun.list, G0_S141_common.fun.list, 
                          G0_R48_common.fun.list, G0_R90_common.fun.list, G0_R141_common.fun.list,
                          G0_RS48_common.fun.list, G0_RS90_common.fun.list, G0_RS141_common.fun.list,
                          G0_BS48_common.fun.list, G0_BS90_common.fun.list, G0_BS141_common.fun.list)

write.csv(moveOTU.list.fun, "Grain to plant OTUs_fungi.csv")


### Calculation of the percentage of moved OTUs
length(G0.bac.list)
length(G0_G141_common.bac)
length(G141.bac.list)


###Bacteria
percent.G141.bac <- length(G0_G141_common.bac)/length(G141.bac.list)
percent.L48.bac <- length(G0_L48_common.bac)/length(L48.bac.list)
percent.L90.bac <- length(G0_L90_common.bac)/length(L90.bac.list)
percent.L141.bac <- length(G0_L141_common.bac)/length(L141.bac.list)
percent.S48.bac <- length(G0_S48_common.bac)/length(S48.bac.list)
percent.S90.bac <- length(G0_S90_common.bac)/length(S90.bac.list)
percent.S141.bac <- length(G0_S141_common.bac)/length(S141.bac.list)
percent.R48.bac <- length(G0_R48_common.bac)/length(R48.bac.list)
percent.R90.bac <- length(G0_R90_common.bac)/length(R90.bac.list)
percent.R141.bac <- length(G0_R141_common.bac)/length(R141.bac.list)
percent.RS48.bac <- length(G0_RS48_common.bac)/length(RS48.bac.list)
percent.RS90.bac <- length(G0_RS90_common.bac)/length(RS90.bac.list)
percent.RS141.bac <- length(G0_RS141_common.bac)/length(RS141.bac.list)
percent.BS48.bac <- length(G0_BS48_common.bac)/length(BS48.bac.list)
percent.BS90.bac <- length(G0_BS90_common.bac)/length(BS90.bac.list)
percent.BS141.bac <- length(G0_BS141_common.bac)/length(BS141.bac.list)


percent.bac <-c(percent.G141.bac, percent.L48.bac, percent.L90.bac, percent.L141.bac,
                percent.S48.bac, percent.S90.bac, percent.S141.bac, percent.R48.bac, percent.R90.bac,
                percent.R141.bac, percent.RS48.bac, percent.RS90.bac, percent.RS141.bac,
                percent.BS48.bac, percent.BS90.bac, percent.BS141.bac) 

###fungi
percent.G141.fun <- length(G0_G141_common.fun)/length(G141.fun.list)
percent.L48.fun <- length(G0_L48_common.fun)/length(L48.fun.list)
percent.L90.fun <- length(G0_L90_common.fun)/length(L90.fun.list)
percent.L141.fun <- length(G0_L141_common.fun)/length(L141.fun.list)
percent.S48.fun <- length(G0_S48_common.fun)/length(S48.fun.list)
percent.S90.fun <- length(G0_S90_common.fun)/length(S90.fun.list)
percent.S141.fun <- length(G0_S141_common.fun)/length(S141.fun.list)
percent.R48.fun <- length(G0_R48_common.fun)/length(R48.fun.list)
percent.R90.fun <- length(G0_R90_common.fun)/length(R90.fun.list)
percent.R141.fun <- length(G0_R141_common.fun)/length(R141.fun.list)
percent.RS48.fun <- length(G0_RS48_common.fun)/length(RS48.fun.list)
percent.RS90.fun <- length(G0_RS90_common.fun)/length(RS90.fun.list)
percent.RS141.fun <- length(G0_RS141_common.fun)/length(RS141.fun.list)
percent.BS48.fun <- length(G0_BS48_common.fun)/length(BS48.fun.list)
percent.BS90.fun <- length(G0_BS90_common.fun)/length(BS90.fun.list)
percent.BS141.fun <- length(G0_BS141_common.fun)/length(BS141.fun.list)


percent.fun <-c(percent.G141.fun, percent.L48.fun, percent.L90.fun, percent.L141.fun,
                percent.S48.fun, percent.S90.fun, percent.S141.fun, percent.R48.fun, percent.R90.fun,
                percent.R141.fun, percent.RS48.fun, percent.RS90.fun, percent.RS141.fun,
                percent.BS48.fun, percent.BS90.fun, percent.BS141.fun) 



bac.sample.list <- c('G141', 'L48', 'L90','L141','S48','S90','S141','R48','R90','R141',
                     'RS48', 'RS90','RS141','BS48', 'BS90','BS141')

df.move.percent<-data.frame(bac.sample.list, percent.bac, percent.fun)

write.xlsx(df.move.percent, "Percentage of moved OTUs in each tissue over time.xlsx")

####Bacteria
df.move.percent.bac<-data.frame(bac.sample.list, percent.bac)
names(df.move.percent.bac)[2] <- "moved"
df.move.percent.bac$non_moved <- 1-df.move.percent.bac$moved

df.move.percent<-read.xlsx("Percentage of moved OTUs in each tissue over time.xlsx",1)
df.move.percent$Stage <- ifelse(df.move.percent$Day == 48, "Tillering", ifelse(df.move.percent$Day == 90, "Heading", "Harvest"))

df.move.percent.bac$Stage <- df.move.percent$Stage
df.move.percent.bac$Compartment <- df.move.percent$Compartment

melt.df.bac<- melt(df.move.percent.bac)
names(melt.df.bac)[1] <- "Sample"

##Fungi
df.move.percent.fun<-data.frame(bac.sample.list, percent.fun)
names(df.move.percent.fun)[2] <- "moved"
df.move.percent.fun$non_moved <- 1-df.move.percent.fun$moved

df.move.percent.fun$Stage <- df.move.percent$Stage
df.move.percent.fun$Compartment <- df.move.percent$Compartment

melt.df.fun<- melt(df.move.percent.fun)
names(melt.df.fun)[1] <- "Sample"

melt.df.bac$Stage <- factor(melt.df.bac$Stage,levels = c('Tillering', 'Heading', 'Harvest'))
melt.df.fun$Stage <- factor(melt.df.fun$Stage,levels = c('Tillering', 'Heading', 'Harvest'))
### bar plot
df.phylum.rel.p1 <- ggplot(melt.df.fun, aes(x=Stage, y = value, fill = variable)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  xlab('')+ theme(aspect.ratio = 1)+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  facet_wrap(~Compartment, scales = "free")+
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1



## Taxonomic composition of seed-borne OTUs in each plant tissue
bac.clean.ss.18
fun.clean.ss.18

G0_G141_common.bac



bac.clean.ss.18.sub.merged.rel <-transform(bac.clean.ss.18.sub.merged, 'compositional')
bac.rel.table<-t(otu_table(bac.clean.ss.18.sub.merged.rel))
df.bac.rel.table <- data.frame(bac.rel.table)

fun.clean.ss.18.sub.merged.rel <-transform(fun.clean.ss.18.sub.merged, 'compositional')
fun.rel.table<-t(otu_table(fun.clean.ss.18.sub.merged.rel))
df.fun.rel.table <- data.frame(fun.rel.table)

### Add taxonomic information
df.bac.rel.table$OTU <- rownames(df.bac.rel.table)
df.bac.rel.tax.table<-merge(df.bac.rel.table, bac.list, by = "OTU")
head(df.bac.rel.tax.table)

df.bac.rel.tax.table$Genus <- as.character(df.bac.rel.tax.table$Genus)

df.bac.rel.tax.table$Genus[is.na(df.bac.rel.tax.table$Genus)] <- "Unassigned"
rownames(df.bac.rel.tax.table) <- df.bac.rel.tax.table$OTU


df.fun.rel.table$OTU <- rownames(df.fun.rel.table)
df.fun.rel.tax.table<-merge(df.fun.rel.table, fun.list, by = "OTU")
head(df.fun.rel.tax.table)

df.fun.rel.tax.table$Genus <- as.character(df.fun.rel.tax.table$Genus)

df.fun.rel.tax.table$Genus[is.na(df.fun.rel.tax.table$Genus)] <- "Unassigned"
rownames(df.fun.rel.tax.table) <- df.fun.rel.tax.table$OTU
###G141
df.bac.rel.G141.table<-subset(df.bac.rel.tax.table[,c("Seed_141","Genus")], rownames(df.bac.rel.tax.table) %in% G0_G141_common.bac)
df.bac.rel.G141.table.gen <- df.bac.rel.G141.table %>% group_by(Genus) %>% summarise(Abund.G.141 = sum(Seed_141))
df.bac.rel.G141.table.gen$Compartment <- "Seed"
df.bac.rel.G141.table.gen$Stage <- "Harvest"

##Leaf
#L48
df.bac.rel.L48.table<-subset(df.bac.rel.tax.table[,c("Leaf_48","Genus")], rownames(df.bac.rel.tax.table) %in% G0_L48_common.bac)
df.bac.rel.L48.table.gen <- df.bac.rel.L48.table %>% group_by(Genus) %>% summarise(Abund.L.48 = sum(Leaf_48))
df.bac.rel.L48.table.gen$Compartment <- "Leaf"
df.bac.rel.L48.table.gen$Stage <- "Tillering"
#L90
df.bac.rel.L90.table<-subset(df.bac.rel.tax.table[,c("Leaf_90","Genus")], rownames(df.bac.rel.tax.table) %in% G0_L90_common.bac)
df.bac.rel.L90.table.gen <- df.bac.rel.L90.table %>% group_by(Genus) %>% summarise(Abund.L.90 = sum(Leaf_90))
df.bac.rel.L90.table.gen$Compartment <- "Leaf"
df.bac.rel.L90.table.gen$Stage <- "Heading"
#L141
df.bac.rel.L141.table<-subset(df.bac.rel.tax.table[,c("Leaf_141","Genus")], rownames(df.bac.rel.tax.table) %in% G0_L141_common.bac)
df.bac.rel.L141.table.gen <- df.bac.rel.L141.table %>% group_by(Genus) %>% summarise(Abund.L.141 = sum(Leaf_141))
df.bac.rel.L141.table.gen$Compartment <- "Leaf"
df.bac.rel.L141.table.gen$Stage <- "Harvest"


##Stem
#S48
df.bac.rel.S48.table<-subset(df.bac.rel.tax.table[,c("Stem_48","Genus")], rownames(df.bac.rel.tax.table) %in% G0_S48_common.bac)
df.bac.rel.S48.table.gen <- df.bac.rel.S48.table %>% group_by(Genus) %>% summarise(Abund.S.48 = sum(Stem_48))
df.bac.rel.S48.table.gen$Compartment <- "Stem"
df.bac.rel.S48.table.gen$Stage <- "Tillering"
#S90
df.bac.rel.S90.table<-subset(df.bac.rel.tax.table[,c("Stem_90","Genus")], rownames(df.bac.rel.tax.table) %in% G0_S90_common.bac)
df.bac.rel.S90.table.gen <- df.bac.rel.S90.table %>% group_by(Genus) %>% summarise(Abund.S.90 = sum(Stem_90))
df.bac.rel.S90.table.gen$Compartment <- "Stem"
df.bac.rel.S90.table.gen$Stage <- "Heading"
#S141
df.bac.rel.S141.table<-subset(df.bac.rel.tax.table[,c("Stem_141","Genus")], rownames(df.bac.rel.tax.table) %in% G0_S141_common.bac)
df.bac.rel.S141.table.gen <- df.bac.rel.S141.table %>% group_by(Genus) %>% summarise(Abund.S.141 = sum(Stem_141))
df.bac.rel.S141.table.gen$Compartment <- "Stem"
df.bac.rel.S141.table.gen$Stage <- "Harvest"


##Root
#R48
df.bac.rel.R48.table<-subset(df.bac.rel.tax.table[,c("Root_48","Genus")], rownames(df.bac.rel.tax.table) %in% G0_R48_common.bac)
df.bac.rel.R48.table.gen <- df.bac.rel.R48.table %>% group_by(Genus) %>% summarise(Abund.R.48 = sum(Root_48))
df.bac.rel.R48.table.gen$Compartment <- "Root"
df.bac.rel.R48.table.gen$Stage <- "Tillering"
#R90
df.bac.rel.R90.table<-subset(df.bac.rel.tax.table[,c("Root_90","Genus")], rownames(df.bac.rel.tax.table) %in% G0_R90_common.bac)
df.bac.rel.R90.table.gen <- df.bac.rel.R90.table %>% group_by(Genus) %>% summarise(Abund.R.90 = sum(Root_90))
df.bac.rel.R90.table.gen$Compartment <- "Root"
df.bac.rel.R90.table.gen$Stage <- "Heading"
#R141
df.bac.rel.R141.table<-subset(df.bac.rel.tax.table[,c("Root_141","Genus")], rownames(df.bac.rel.tax.table) %in% G0_R141_common.bac)
df.bac.rel.R141.table.gen <- df.bac.rel.R141.table %>% group_by(Genus) %>% summarise(Abund.R.141 = sum(Root_141))
df.bac.rel.R141.table.gen$Compartment <- "Root"
df.bac.rel.R141.table.gen$Stage <- "Harvest"


##Rhizosphere
#RS48
df.bac.rel.RS48.table<-subset(df.bac.rel.tax.table[,c("Rhizosphere_48","Genus")], rownames(df.bac.rel.tax.table) %in% G0_RS48_common.bac)
df.bac.rel.RS48.table.gen <- df.bac.rel.RS48.table %>% group_by(Genus) %>% summarise(Abund.RS.48 = sum(Rhizosphere_48))
df.bac.rel.RS48.table.gen$Compartment <- "Rhizosphere"
df.bac.rel.RS48.table.gen$Stage <- "Tillering"
#RS90
df.bac.rel.RS90.table<-subset(df.bac.rel.tax.table[,c("Rhizosphere_90","Genus")], rownames(df.bac.rel.tax.table) %in% G0_RS90_common.bac)
df.bac.rel.RS90.table.gen <- df.bac.rel.RS90.table %>% group_by(Genus) %>% summarise(Abund.RS.90 = sum(Rhizosphere_90))
df.bac.rel.RS90.table.gen$Compartment <- "Rhizosphere"
df.bac.rel.RS90.table.gen$Stage <- "Heading"
#RS141
df.bac.rel.RS141.table<-subset(df.bac.rel.tax.table[,c("Rhizosphere_141","Genus")], rownames(df.bac.rel.tax.table) %in% G0_RS141_common.bac)
df.bac.rel.RS141.table.gen <- df.bac.rel.RS141.table %>% group_by(Genus) %>% summarise(Abund.RS.141 = sum(Rhizosphere_141))
df.bac.rel.RS141.table.gen$Compartment <- "Rhizosphere"
df.bac.rel.RS141.table.gen$Stage <- "Harvest"


##Bulk_soil
#BS48
df.bac.rel.BS48.table<-subset(df.bac.rel.tax.table[,c("Bulk_soil_48","Genus")], rownames(df.bac.rel.tax.table) %in% G0_BS48_common.bac)
df.bac.rel.BS48.table.gen <- df.bac.rel.BS48.table %>% group_by(Genus) %>% summarise(Abund.BS.48 = sum(Bulk_soil_48))
df.bac.rel.BS48.table.gen$Compartment <- "Bulk_soil"
df.bac.rel.BS48.table.gen$Stage <- "Tillering"
#BS90
df.bac.rel.BS90.table<-subset(df.bac.rel.tax.table[,c("Bulk_soil_90","Genus")], rownames(df.bac.rel.tax.table) %in% G0_BS90_common.bac)
df.bac.rel.BS90.table.gen <- df.bac.rel.BS90.table %>% group_by(Genus) %>% summarise(Abund.BS.90 = sum(Bulk_soil_90))
df.bac.rel.BS90.table.gen$Compartment <- "Bulk_soil"
df.bac.rel.BS90.table.gen$Stage <- "Heading"
#BS141
df.bac.rel.BS141.table<-subset(df.bac.rel.tax.table[,c("Bulk_soil_141","Genus")], rownames(df.bac.rel.tax.table) %in% G0_BS141_common.bac)
df.bac.rel.BS141.table.gen <- df.bac.rel.BS141.table %>% group_by(Genus) %>% summarise(Abund.BS.141 = sum(Bulk_soil_141))
df.bac.rel.BS141.table.gen$Compartment <- "Bulk_soil"
df.bac.rel.BS141.table.gen$Stage <- "Harvest"

melt.bac.rel.G141.table.gen <- melt(df.bac.rel.G141.table.gen)
melt.bac.rel.L48.table.gen <- melt(df.bac.rel.L48.table.gen)
melt.bac.rel.L90.table.gen <- melt(df.bac.rel.L90.table.gen)
melt.bac.rel.L141.table.gen <- melt(df.bac.rel.L141.table.gen)
melt.bac.rel.S48.table.gen <- melt(df.bac.rel.S48.table.gen)
melt.bac.rel.S90.table.gen <- melt(df.bac.rel.S90.table.gen)
melt.bac.rel.S141.table.gen <- melt(df.bac.rel.S141.table.gen)
melt.bac.rel.R48.table.gen <- melt(df.bac.rel.R48.table.gen)
melt.bac.rel.R90.table.gen <- melt(df.bac.rel.R90.table.gen)
melt.bac.rel.R141.table.gen <- melt(df.bac.rel.R141.table.gen)
melt.bac.rel.RS48.table.gen <- melt(df.bac.rel.RS48.table.gen)
melt.bac.rel.RS90.table.gen <- melt(df.bac.rel.RS90.table.gen)
melt.bac.rel.RS141.table.gen <- melt(df.bac.rel.RS141.table.gen)
melt.bac.rel.BS48.table.gen <- melt(df.bac.rel.BS48.table.gen)
melt.bac.rel.BS90.table.gen <- melt(df.bac.rel.BS90.table.gen)
melt.bac.rel.BS141.table.gen <- melt(df.bac.rel.BS141.table.gen)

melt.bac.rel.table.gen<-rbind(melt.bac.rel.G141.table.gen,melt.bac.rel.L48.table.gen,
      melt.bac.rel.L90.table.gen,melt.bac.rel.L141.table.gen,
      melt.bac.rel.S48.table.gen,melt.bac.rel.S90.table.gen,
      melt.bac.rel.S141.table.gen,melt.bac.rel.R48.table.gen,
      melt.bac.rel.R90.table.gen, melt.bac.rel.R141.table.gen,
      melt.bac.rel.RS48.table.gen,melt.bac.rel.RS90.table.gen,
      melt.bac.rel.RS141.table.gen, melt.bac.rel.BS48.table.gen,
      melt.bac.rel.BS90.table.gen,melt.bac.rel.BS141.table.gen)

melt.bac.rel.table.gen.ord <- melt.bac.rel.table.gen %>% group_by(Genus) %>% summarise(total = sum(value)) %>% arrange(total)
melt.bac.rel.table.gen.ord$Genus[which(melt.bac.rel.table.gen.ord$total < 0.03)] <- "Others"
vec.gen<-melt.bac.rel.table.gen.ord$Genus
vec.order <- vec.gen[-which(vec.gen%in%c("Unassigned", "Others"))]
vec.uniden.Low <- c("Unassigned", "Others")
vec.reorder <- append(vec.uniden.Low,vec.order)

melt.bac.rel.table.gen$Genus[which(!(melt.bac.rel.table.gen$Genus %in% vec.reorder))] <- "Others"


melt.bac.rel.table.gen$Genus <- factor(melt.bac.rel.table.gen$Genus, levels = vec.reorder) 
melt.bac.rel.table.gen$Stage <- factor(melt.bac.rel.table.gen$Stage, levels = c("Tillering","Heading", "Harvest")) 

### bar plot
df.phylum.rel.p1 <- ggplot(melt.bac.rel.table.gen, aes(x=Stage, y = value, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  xlab('')+ theme(aspect.ratio = 1)+
  ylab("Relative abundance\n") +
  facet_wrap(~Compartment, scales = "free")+
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = F))+
  theme(legend.position="bottom") +
  coord_cartesian(ylim = c(0, 1))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1


### Cumulative abundance table
melt.bac.rel.table.gen.tab <-melt.bac.rel.table.gen%>% group_by(Compartment, Stage) %>% summarise(sum_abund = sum(value))


###Fungi
###G141
df.fun.rel.G141.table<-subset(df.fun.rel.tax.table[,c("Grain_141","Genus")], rownames(df.fun.rel.tax.table) %in% G0_G141_common.fun)
df.fun.rel.G141.table.gen <- df.fun.rel.G141.table %>% group_by(Genus) %>% summarise(Abund.G.141 = sum(Grain_141))
df.fun.rel.G141.table.gen$Compartment <- "Seed"
df.fun.rel.G141.table.gen$Stage <- "Harvest"

##Leaf
#L48
df.fun.rel.L48.table<-subset(df.fun.rel.tax.table[,c("Leaf_48","Genus")], rownames(df.fun.rel.tax.table) %in% G0_L48_common.fun)
df.fun.rel.L48.table.gen <- df.fun.rel.L48.table %>% group_by(Genus) %>% summarise(Abund.L.48 = sum(Leaf_48))
df.fun.rel.L48.table.gen$Compartment <- "Leaf"
df.fun.rel.L48.table.gen$Stage <- "Tillering"
#L90
df.fun.rel.L90.table<-subset(df.fun.rel.tax.table[,c("Leaf_90","Genus")], rownames(df.fun.rel.tax.table) %in% G0_L90_common.fun)
df.fun.rel.L90.table.gen <- df.fun.rel.L90.table %>% group_by(Genus) %>% summarise(Abund.L.90 = sum(Leaf_90))
df.fun.rel.L90.table.gen$Compartment <- "Leaf"
df.fun.rel.L90.table.gen$Stage <- "Heading"
#L141
df.fun.rel.L141.table<-subset(df.fun.rel.tax.table[,c("Leaf_141","Genus")], rownames(df.fun.rel.tax.table) %in% G0_L141_common.fun)
df.fun.rel.L141.table.gen <- df.fun.rel.L141.table %>% group_by(Genus) %>% summarise(Abund.L.141 = sum(Leaf_141))
df.fun.rel.L141.table.gen$Compartment <- "Leaf"
df.fun.rel.L141.table.gen$Stage <- "Harvest"


##Stem
#S48
df.fun.rel.S48.table<-subset(df.fun.rel.tax.table[,c("Stem_48","Genus")], rownames(df.fun.rel.tax.table) %in% G0_S48_common.fun)
df.fun.rel.S48.table.gen <- df.fun.rel.S48.table %>% group_by(Genus) %>% summarise(Abund.S.48 = sum(Stem_48))
df.fun.rel.S48.table.gen$Compartment <- "Stem"
df.fun.rel.S48.table.gen$Stage <- "Tillering"
#S90
df.fun.rel.S90.table<-subset(df.fun.rel.tax.table[,c("Stem_90","Genus")], rownames(df.fun.rel.tax.table) %in% G0_S90_common.fun)
df.fun.rel.S90.table.gen <- df.fun.rel.S90.table %>% group_by(Genus) %>% summarise(Abund.S.90 = sum(Stem_90))
df.fun.rel.S90.table.gen$Compartment <- "Stem"
df.fun.rel.S90.table.gen$Stage <- "Heading"
#S141
df.fun.rel.S141.table<-subset(df.fun.rel.tax.table[,c("Stem_141","Genus")], rownames(df.fun.rel.tax.table) %in% G0_S141_common.fun)
df.fun.rel.S141.table.gen <- df.fun.rel.S141.table %>% group_by(Genus) %>% summarise(Abund.S.141 = sum(Stem_141))
df.fun.rel.S141.table.gen$Compartment <- "Stem"
df.fun.rel.S141.table.gen$Stage <- "Harvest"


##Root
#R48
df.fun.rel.R48.table<-subset(df.fun.rel.tax.table[,c("Root_48","Genus")], rownames(df.fun.rel.tax.table) %in% G0_R48_common.fun)
df.fun.rel.R48.table.gen <- df.fun.rel.R48.table %>% group_by(Genus) %>% summarise(Abund.R.48 = sum(Root_48))
df.fun.rel.R48.table.gen$Compartment <- "Root"
df.fun.rel.R48.table.gen$Stage <- "Tillering"
#R90
df.fun.rel.R90.table<-subset(df.fun.rel.tax.table[,c("Root_90","Genus")], rownames(df.fun.rel.tax.table) %in% G0_R90_common.fun)
df.fun.rel.R90.table.gen <- df.fun.rel.R90.table %>% group_by(Genus) %>% summarise(Abund.R.90 = sum(Root_90))
df.fun.rel.R90.table.gen$Compartment <- "Root"
df.fun.rel.R90.table.gen$Stage <- "Heading"
#R141
df.fun.rel.R141.table<-subset(df.fun.rel.tax.table[,c("Root_141","Genus")], rownames(df.fun.rel.tax.table) %in% G0_R141_common.fun)
df.fun.rel.R141.table.gen <- df.fun.rel.R141.table %>% group_by(Genus) %>% summarise(Abund.R.141 = sum(Root_141))
df.fun.rel.R141.table.gen$Compartment <- "Root"
df.fun.rel.R141.table.gen$Stage <- "Harvest"


##Rhizosphere
#RS48
df.fun.rel.RS48.table<-subset(df.fun.rel.tax.table[,c("Rhizosphere_48","Genus")], rownames(df.fun.rel.tax.table) %in% G0_RS48_common.fun)
df.fun.rel.RS48.table.gen <- df.fun.rel.RS48.table %>% group_by(Genus) %>% summarise(Abund.RS.48 = sum(Rhizosphere_48))
df.fun.rel.RS48.table.gen$Compartment <- "Rhizosphere"
df.fun.rel.RS48.table.gen$Stage <- "Tillering"
#RS90
df.fun.rel.RS90.table<-subset(df.fun.rel.tax.table[,c("Rhizosphere_90","Genus")], rownames(df.fun.rel.tax.table) %in% G0_RS90_common.fun)
df.fun.rel.RS90.table.gen <- df.fun.rel.RS90.table %>% group_by(Genus) %>% summarise(Abund.RS.90 = sum(Rhizosphere_90))
df.fun.rel.RS90.table.gen$Compartment <- "Rhizosphere"
df.fun.rel.RS90.table.gen$Stage <- "Heading"
#RS141
df.fun.rel.RS141.table<-subset(df.fun.rel.tax.table[,c("Rhizosphere_141","Genus")], rownames(df.fun.rel.tax.table) %in% G0_RS141_common.fun)
df.fun.rel.RS141.table.gen <- df.fun.rel.RS141.table %>% group_by(Genus) %>% summarise(Abund.RS.141 = sum(Rhizosphere_141))
df.fun.rel.RS141.table.gen$Compartment <- "Rhizosphere"
df.fun.rel.RS141.table.gen$Stage <- "Harvest"


##Bulk_soil
#BS48
df.fun.rel.BS48.table<-subset(df.fun.rel.tax.table[,c("Bulk_soil_48","Genus")], rownames(df.fun.rel.tax.table) %in% G0_BS48_common.fun)
df.fun.rel.BS48.table.gen <- df.fun.rel.BS48.table %>% group_by(Genus) %>% summarise(Abund.BS.48 = sum(Bulk_soil_48))
df.fun.rel.BS48.table.gen$Compartment <- "Bulk_soil"
df.fun.rel.BS48.table.gen$Stage <- "Tillering"
#BS90
df.fun.rel.BS90.table<-subset(df.fun.rel.tax.table[,c("Bulk_soil_90","Genus")], rownames(df.fun.rel.tax.table) %in% G0_BS90_common.fun)
df.fun.rel.BS90.table.gen <- df.fun.rel.BS90.table %>% group_by(Genus) %>% summarise(Abund.BS.90 = sum(Bulk_soil_90))
df.fun.rel.BS90.table.gen$Compartment <- "Bulk_soil"
df.fun.rel.BS90.table.gen$Stage <- "Heading"
#BS141
df.fun.rel.BS141.table<-subset(df.fun.rel.tax.table[,c("Bulk_soil_141","Genus")], rownames(df.fun.rel.tax.table) %in% G0_BS141_common.fun)
df.fun.rel.BS141.table.gen <- df.fun.rel.BS141.table %>% group_by(Genus) %>% summarise(Abund.BS.141 = sum(Bulk_soil_141))
df.fun.rel.BS141.table.gen$Compartment <- "Bulk_soil"
df.fun.rel.BS141.table.gen$Stage <- "Harvest"

melt.fun.rel.G141.table.gen <- melt(df.fun.rel.G141.table.gen)
melt.fun.rel.L48.table.gen <- melt(df.fun.rel.L48.table.gen)
melt.fun.rel.L90.table.gen <- melt(df.fun.rel.L90.table.gen)
melt.fun.rel.L141.table.gen <- melt(df.fun.rel.L141.table.gen)
melt.fun.rel.S48.table.gen <- melt(df.fun.rel.S48.table.gen)
melt.fun.rel.S90.table.gen <- melt(df.fun.rel.S90.table.gen)
melt.fun.rel.S141.table.gen <- melt(df.fun.rel.S141.table.gen)
melt.fun.rel.R48.table.gen <- melt(df.fun.rel.R48.table.gen)
melt.fun.rel.R90.table.gen <- melt(df.fun.rel.R90.table.gen)
melt.fun.rel.R141.table.gen <- melt(df.fun.rel.R141.table.gen)
melt.fun.rel.RS48.table.gen <- melt(df.fun.rel.RS48.table.gen)
melt.fun.rel.RS90.table.gen <- melt(df.fun.rel.RS90.table.gen)
melt.fun.rel.RS141.table.gen <- melt(df.fun.rel.RS141.table.gen)
melt.fun.rel.BS48.table.gen <- melt(df.fun.rel.BS48.table.gen)
melt.fun.rel.BS90.table.gen <- melt(df.fun.rel.BS90.table.gen)
melt.fun.rel.BS141.table.gen <- melt(df.fun.rel.BS141.table.gen)

melt.fun.rel.table.gen<-rbind(melt.fun.rel.G141.table.gen,melt.fun.rel.L48.table.gen,
                              melt.fun.rel.L90.table.gen,melt.fun.rel.L141.table.gen,
                              melt.fun.rel.S48.table.gen,melt.fun.rel.S90.table.gen,
                              melt.fun.rel.S141.table.gen,melt.fun.rel.R48.table.gen,
                              melt.fun.rel.R90.table.gen, melt.fun.rel.R141.table.gen,
                              melt.fun.rel.RS48.table.gen,melt.fun.rel.RS90.table.gen,
                              melt.fun.rel.RS141.table.gen, melt.fun.rel.BS48.table.gen,
                              melt.fun.rel.BS90.table.gen,melt.fun.rel.BS141.table.gen)

melt.fun.rel.table.gen.ord <- melt.fun.rel.table.gen %>% group_by(Genus) %>% summarise(total = sum(value)) %>% arrange(total)
melt.fun.rel.table.gen.ord$Genus[which(melt.fun.rel.table.gen.ord$total < 0.03)] <- "Others"
melt.fun.rel.table.gen.ord$Genus[which(melt.fun.rel.table.gen.ord$Genus == "unidentified")] <- "Unassigned"
vec.gen<-melt.fun.rel.table.gen.ord$Genus
vec.order <- vec.gen[-which(vec.gen%in%c("Unassigned", "Others"))]
vec.uniden.Low <- c("Unassigned", "Others")
vec.reorder <- append(vec.uniden.Low,vec.order)

melt.fun.rel.table.gen$Genus[which(melt.fun.rel.table.gen$Genus == "unidentified")] <- "Unassigned"
melt.fun.rel.table.gen$Genus[which(!(melt.fun.rel.table.gen$Genus %in% vec.reorder))] <- "Others"


melt.fun.rel.table.gen$Genus <- factor(melt.fun.rel.table.gen$Genus, levels = vec.reorder) 
melt.fun.rel.table.gen$Stage <- factor(melt.fun.rel.table.gen$Stage, levels = c("Tillering","Heading", "Harvest")) 

### bar plot
df.phylum.rel.p1 <- ggplot(melt.fun.rel.table.gen, aes(x=Stage, y = value, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  xlab('')+ theme(aspect.ratio = 1)+
  ylab("Relative abundance\n") +
  facet_wrap(~Compartment, scales = "free")+
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = F))+
  theme(legend.position="bottom") +
  coord_cartesian(ylim = c(0, 1))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1


### Cumulative abundance table
melt.fun.rel.table.gen.tab <-melt.fun.rel.table.gen%>% group_by(Compartment, Stage) %>% summarise(sum_abund = sum(value))


write.xlsx(melt.fun.rel.table.gen, 'raw data table for fungal relative abundance.xlsx')
write.xlsx(melt.bac.rel.table.gen, 'raw data table for bacterial relative abundance.xlsx')
