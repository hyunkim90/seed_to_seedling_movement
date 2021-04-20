### 

# Set colors for plotting

### Bacteria

#### get unassigned vectors
#### get CP and MT phyloseq obj and vectors
# (1) CP
phy.cp <- subset_taxa(phy, Class == "D_2__Chloroplast") ## just confirming
phy.cp <- subset_taxa(phy, Order == "D_3__Chloroplast")
vec.cp <- rownames(otu_table(phy.cp))
length(rownames(otu_table(phy.cp))) ## 8 otus of CP
vec.cp

# (2) MT
phy.mt <- subset_taxa(phy, Family == "D_4__Mitochondria")
phy.mt <- subset_taxa(phy, Order == "D_3__Rickettsiales") #337
vec.mt <- rownames(otu_table(phy.mt))
tax_table(phy.mt)
length(rownames(otu_table(phy.mt))) ## 8 otus of CP

# (3) Unassigned
unique(tax_table(phy)[,'Kingdom']) ## only bacteria, then no need to exclude
phy.un <- subset_taxa(phy, Kingdom == "Unassigned")
vec.un <- rownames(otu_table(phy.un))
tax_table(phy.un)
length(rownames(otu_table(phy.un))) ## 1 otus

### exclude those vectors
phy #1346 taxa and 129 samples

sample_names(phy)
sample_variables(phy)

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### let's do it!!!

phy.clean <- pop_taxa(phy, c(vec.mt,vec.cp))

##before clean-up
phy 
sum(otu_table(phy)) ##926870

##after clean-up
phy.clean
sum(otu_table(phy.clean)) #21148

# checking procedure of whether the MT and CP otus are cleaned
# taxa_names(subset_taxa(phy, Order == "D_3__Rickettsiales"))  # before 337
# taxa_names(subset_taxa(phy, Order=="D_3__Chloroplast")) # before 130
# 
# taxa_names(subset_taxa(phy.clean , Order == "D_3__Rickettsiales")) # after 0
# taxa_names(subset_taxa(phy.clean , Family=="D_4__Mitochondria")) # after 0
# taxa_names(subset_taxa(phy.clean , Order == "D_3__Chloroplast")) # after 0

tax_table(phy.clean)

#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(phy.two)[,colnames(tax_table(phy.two))] <- gsub(tax_table(phy.two)[,colnames(tax_table(phy.two))],pattern="[A-Z]_[0-9]__",replacement="")
# phy.test4 <- phy.two %>% psmelt() 
# phy.test4

tax_table(phy.clean)[,colnames(tax_table(phy.clean))] <- gsub(tax_table(phy.clean)[,colnames(tax_table(phy.clean))],pattern="[A-Z]_[0-9]__",replacement="")

#' sample_data(phy.clean)$SampleID <- factor(sample_data(phy.clean)$SampleID, levels =target_PAB)


tax_table(phy.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

phy.clean    ## 28330
str(phy.clean)
otu_table(phy.clean)
summarize_phyloseq(phy.clean)
phy.clean  ## 28330


# ## fix it in phy.clean object!!! pop_taxa does the work
# phy.clean <- pop_taxa(phy.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(phy.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement 
phy.clean.otu <- otu_table(phy.clean)
head(phy.clean.otu)
df.clean.otu <- data.frame(phy.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


sample_names(phy.clean)


## Exclude false-positive OTUs
negative <- df.clean.otu[,c('OTU','Negative')]
negative_0 <- subset(negative, Negative > 0)
neg.otu <- negative_0$OTU  ### 69 otu to be eliminated
neg.otu
length(neg.otu)
tax_table(phy.clean)[neg.otu]
library(xlsx)
write.xlsx(tax_table(phy.clean)[neg.otu], 'neg.otu_for_fungi_taxonomy.xlsx')

phy.clean.ss <- pop_taxa(phy.clean, c(neg.otu))
summarize_phyloseq(phy.clean.ss)
## Remove reads with over 464 bp
bac.seq <- read.fasta(file = "dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")

colSums(otu_table(phy.clean.ss))
#length
getLength(bac.seq)

## max length and min length
# getLength(bac.seq)
# min(getLength(bac.seq)) #39
# max(getLength(bac.seq)) #447
# 
# ##
# bac.seq[which(getLength(bac.seq)<400)]
# otu_less_400bp <- attr(bac.seq[which(getLength(bac.seq)<400)], "names")
# 
# phy.clean
# phy.clean.ss <- pop_taxa(phy.clean,otu_less_400bp)
# phy.clean.ss  ## 1343
# sum(otu_table(phy.clean)) #1835114
# sum(otu_table(phy.clean.ss)) #1832965


## Remove Negative (number of read is 0) from OTU table
sample_names(phy.clean.ss)
samples.wo.neg<- sample_names(phy.clean.ss)[-which(sample_names(phy.clean.ss)=="Negative")]

otu_table(phy.clean.ss) <- otu_table(phy.clean.ss)[,c(samples.wo.neg)]


sample_names(phy.clean)
samples.wo.neg<- sample_names(phy.clean)[-which(sample_names(phy.clean)=="Negative")]

otu_table(phy.clean) <- otu_table(phy.clean)[,c(samples.wo.neg)]

sum(otu_table(phy.clean.ss)) #1831386
sum(otu_table(phy.clean))


bac.clean.ss<-phy.clean

summarize_phyloseq(bac.clean.ss)
colSums(otu_table(phy.clean.ss))
colSums(otu_table(phy.clean))

###Fungal community
#### get unassigned vectors

# (3) Unassigned
unique(tax_table(fun)[,'Kingdom']) ## "k__Chromista", "k__Plantae", "Unassigned"
tax_table(fun)[,'Kingdom']
fun.un <- subset_taxa(fun, Kingdom %in% c("Unassigned","k__Chromista","k__Plantae"))
vec.un <- rownames(otu_table(fun.un))
tax_table(fun.un)
length(rownames(otu_table(fun.un))) ##  43

### exclude those vectors
fun  # 135 samples, 1629 taxa/ 1641 taxa

sample_names(fun)
sample_variables(fun)

### pop taxa application
## get rid of CP and MT otus
## pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### Clean it up!!!
fun.clean <- pop_taxa(fun, c(vec.un))

fun.clean #135 samples 1598 OTUs

#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(fun.two)[,colnames(tax_table(fun.two))] <- gsub(tax_table(fun.two)[,colnames(tax_table(fun.two))],pattern="[A-Z]_[0-9]__",replacement="")
# fun.test4 <- fun.two %>% psmelt() 
# fun.test4


tax_table(fun.clean)[,colnames(tax_table(fun.clean))] <- gsub(tax_table(fun.clean)[,colnames(tax_table(fun.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.clean)$SampleID <- factor(sample_data(fun.clean)$SampleID, levels =target_PAB)

tax_table(fun.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.clean 
str(fun.clean)
otu_table(fun.clean)

# ## fix it in fun.clean object!!! pop_taxa does the work
# fun.clean <- pop_taxa(fun.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement 
fun.clean.otu <- otu_table(fun.clean)
head(fun.clean.otu)
df.clean.otu <- data.frame(fun.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


# # how about we delete major taxa in the negative sequence?
# negative <- df.clean.otu[,c("OTU",'F204','F205','F206')]
# negative
# negative$total <- apply(negative[,-1], 1, sum)
# negative_0 <- subset(negative, negative$total > 0) 
# neg.otu <- negative_0$OTU  ### 58 otu to be eliminated
# neg.otu
# length(neg.otu)
# tax_table(fun.clean)[neg.otu]
# library(xlsx)
# write.xlsx(tax_table(fun.clean)[neg.otu], 'fungal neg.otu_taxonomy_old DB.xlsx')
# 

sum(otu_table(fun)) ##3154736
sum(otu_table(fun.clean)) ## 2983621 # calibrated 3154052 (2020 DB)/3051164(2017 DB)
sum(otu_table(fun.clean.ss)) #365409


##### get rid of otu of less than 100 reads

fun.seq <- read.fasta(file = "dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")

otu_less_than_100bp <- attr(fun.seq[which(getLength(fun.seq)<100)], "names")
fun.clean.ss
fun.clean.ss <- pop_taxa(fun.clean,otu_less_than_100bp)
fun.clean.ss ## 9437
sum(otu_table(fun.clean.ss)) ## 12191641

## Fungal OTUs which are 'unidentifed' at the phylum level
tax.fun<-tax_table(fun.clean.ss)
tax.fun <-data.frame(tax.fun)
unidentified.phylum<-rownames(tax.fun)[is.na(tax.fun$Phylum)]

unidentified.phylum.seq<-fun.seq[which(attr(fun.seq, "names")%in% unidentified.phylum)]
write.fasta(unidentified.phylum.seq, names(unidentified.phylum.seq), file.out = "unidentified.phylum.seq.fasta")

## Remove unidentified OTUs based on BLAST result
remove.fungal.otu <-c("20822ac1d45e6ed3dbb18a574b874b3e","30714dcecb4d0ded7d98793b44b7517d","357677e51a57c2a835a9440aa438ceb5",
                      "357677e51a57c2a835a9440aa438ceb5","5cb6a221d2477ea04dd9a9c8a8c9f848","695dfd8e4019d0187fda0d6d99060f13",
                      "76f3cf5a45b1ce2d24d68e89e0a79395","951962f1064faf3a655c974f57c85692","951962f1064faf3a655c974f57c85692",
                      "9f20d7f91556a373e8c15e6e3396be04","a69bd69c71d1ffee1304b62e6859a74a","b0c1eb3e17c35d916401e15149b1b40a",
                      "c616ecad6cd99147100d9b44346553c4","d41bb218d5b6766fc0903db1f1ab7cdb") 
fun.clean.ss <- pop_taxa(fun.clean.ss, remove.fungal.otu)
 ## Designating OTU id
bac.list <- bac.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
fun.list <- fun.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

bac.list$number <- paste0('B',1:dim(bac.list)[1])
bac.list

bac.list$OTU_id <- ifelse(is.na(bac.list$Genus),ifelse(is.na(bac.list$Family),paste0(bac.list$number,'_o_',bac.list$Order),paste0(bac.list$number,'_f_',bac.list$Family)),paste0(bac.list$number,'_',bac.list$Genus))
bac.list$OTU_id

fun.list$number <- paste0('F',1:dim(fun.list)[1])
fun.list

fun.list$OTU_id <- ifelse(is.na(fun.list$Genus),ifelse(is.na(fun.list$Family),paste0(fun.list$number,'_o_',fun.list$Order),paste0(fun.list$number,'_f_',fun.list$Family)),paste0(fun.list$number,'_',fun.list$Genus))
fun.list$OTU_id


bac.list

fun.list

otu.list <- rbind(bac.list, arch.list, fun.list)
dim(otu.list)
# write.xlsx(otu.list,'otu_id_book.xlsx')
fun.list.its1.OTU_id <-fun.list.its1[c('OTU','OTU_id')]


OTU_id.list <- rbind(bac.list[c('OTU','OTU_id')],arch.list[c('OTU','OTU_id')],fun.list[c('OTU','OTU_id')])
OTU_id.list$OTU_id
