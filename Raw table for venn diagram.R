### Raw table for venn diagram
## Add taxonomic information
## Raw table
bac.tab.unhulled<-read.csv("OTU_list_bac_0.7.14.csv")
fun.tab.unhulled<-read.csv("OTU_list_fun_0.7.14.csv")

bac.tab.hulled<-read.csv("OTH_list_bac_0.7.14_hulled.csv")
fun.tab.hulled<-read.csv("OTH_list_fun_0.7.14_hulled.csv")


names(bac.tab.hulled)[2] <- "OTU"
names(bac.tab.hulled)[3] <- "Compartment"

names(fun.tab.hulled)[2] <- "OTU"
names(fun.tab.hulled)[3] <- "Compartment"

names(bac.tab.unhulled)[2] <- "OTU"
names(bac.tab.unhulled)[3] <- "Compartment"

names(fun.tab.unhulled)[2] <- "OTU"
names(fun.tab.unhulled)[3] <- "Compartment"


bac.list
fun.list

bac.tab.hulled.merged <- merge(bac.tab.hulled, bac.list, by = "OTU")
bac.tab.unhulled.merged <- merge(bac.tab.unhulled, bac.list, by = "OTU")


fun.tab.hulled.merged <- merge(fun.tab.hulled, fun.list, by = "OTU")
fun.tab.unhulled.merged <- merge(fun.tab.unhulled, fun.list, by = "OTU")


write.csv(bac.tab.hulled.merged, "bac.tab.hulled.merged.csv")
write.csv(bac.tab.unhulled.merged, "bac.tab.unhulled.merged.csv")

write.csv(fun.tab.hulled.merged, "fun.tab.hulled.merged.csv")
write.csv(fun.tab.unhulled.merged, "fun.tab.unhulled.merged.csv")




### Trimmed table

library(xlsx)
OTU.list <- rbind(bac.list, fun.list)

trimmed.table<-read.xlsx("Venn diagram taxa table.xlsx",1)

trimmed.table.taxa<-merge(trimmed.table, OTU.list, by = "OTU")
write.xlsx(trimmed.table.taxa,"trimmed.table.taxa.xlsx")
