###PERMANOVA
## ALL
bac.clean.log
b.otu <- otu_table(bac.clean.log)
b.meta <- sample_data(bac.clean.log)
b.meta <- data.frame(b.meta)

b.meta$Age <- as.factor(b.meta$Age)
b.permanova <- adonis(formula = t(b.otu) ~ (Hull+Compartment+Age), data = b.meta, permutations=9999, method = "bray")
b.permanova


fun.clean.log
f.otu <- otu_table(fun.clean.log)
f.meta <- sample_data(fun.clean.log)
f.meta <- data.frame(f.meta)
f.meta$Organ[f.meta$Organ == "Seed"] <- "Grain"
names(f.meta)[3] <- "Compartment"

f.meta$Age <- as.factor(f.meta$Age)
f.permanova <- adonis(formula = t(f.otu) ~ (Hull+Compartment+Age), data = f.meta, permutations=9999, method = "bray")
f.permanova
