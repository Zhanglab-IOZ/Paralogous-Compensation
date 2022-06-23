
###Travor Hart's latest GenomeBiology data. CORUM same protein complex
library(reshape2);library(dplyr)
library(gridExtra);library(ggbeeswarm);library(ggsignif)

corum <- read.table("../datafiles/CORUM_allComplexes.txt",sep ="\t")
corum <- corum %>% inner_join(corum, by="V1")
names(corum) <- c("complex", "genename1", "genename2")
corum <- corum %>% filter(genename1 != genename2) %>% mutate(pairs = paste(genename1, genename2, sep="_"))

Hart2 <- read.table("../datafiles/Hart_Supplementary_Table2.txt", header=T)
Hart2$gene1 <- sub("_.*", "", row.names(Hart2))
Hart2$gene2 <- sub(".*_", "", row.names(Hart2))
Hart2 <- Hart2 %>% rename( A549.gi = A549, HT29.gi = HT29, OVCAR8.gi = OVCAR8)
Hart2 <- Hart2 %>% mutate(gi = ifelse(A549.gi < -3 | HT29.gi < -3 | OVCAR8.gi < -3, "Negative", "others")) %>% mutate( sum = A549.gi + HT29.gi + OVCAR8.gi ) %>% arrange(gi, sum )
Hart2$corum <- ifelse(paste(Hart2$gene1, Hart2$gene2, sep="_") %in% corum$pairs, "Same","Not.same")


###NBT's latest data. Ensembl paralogous family size
compare.gi    <- read.csv("../datafiles/Compare.gi.20200912novel.csv", header=T, row.names=1)
compare.gi    <- compare.gi %>% filter(V1 == "choose") %>% rename(gene1 = genename.x, gene2 = genename.y)
compare.gi    <- compare.gi %>% mutate(gi = ifelse(gi == "negative GI", "Negative", "others"))

Parrish.raw    <- read.csv("../datafiles/Phoebe.2021BioRxiv.csv", header=T, row.names=1)
Parrish.raw    <- Parrish.raw %>% select(target1, target2, PC9_GI_flag, HeLa_GI_flag, GI_flag) %>% rename(gene1 = target1, gene2 = target2) %>% mutate(significant.x = ifelse(PC9_GI_flag == "SL_in_PC9", "Y", "N"), significant.y = ifelse(HeLa_GI_flag == "SL_in_HeLa", "Y", "N")) %>% mutate(gi = ifelse(GI_flag == "synthetic_lethal", "Negative", "others"))
Parrish        <- Parrish.raw %>% select(gene1, gene2, gi)

nbt           <- read.csv("../datafiles/2020-NatureBiotech.raw.data.csv", header = T, row.names = 1)
nbt           <- nbt %>% mutate(gi = ifelse(hap1.negative.gi == "Negative" | rpe1.negative.gi == "Negative", "Negative", "others"))

for(i in 1:dim(compare.gi)[2]){if( class(compare.gi[,i]) == "factor" ){ compare.gi[,i] <- as.character(compare.gi[,i]) }}
for(i in 1:dim(Hart2)[2]){if( class(Hart2[,i]) == "factor" ){ Hart2[,i] <- as.character(Hart2[,i]) }}
for(i in 1:dim(Parrish)[2]){if( class(Parrish[,i]) == "factor" ){ Parrish[,i] <- as.character(Parrish[,i]) }}
#for(i in 1:dim(nbt)[2]){if( class(nbt[,i]) == "factor" ){ nbt[,i] <- as.character(nbt[,i]) }}

Thompson <- read.csv("../datafiles/2021.NC.SData2.csv", stringsAsFactors = F)
Thompson <- Thompson[Thompson[,3] == "Paralogous_gene_pair", ]
Thompson <- data.frame("pair" = Thompson[,2], "gene1" = sub("_.*", "", Thompson[,2]), "gene2" = sub(".*_", "", Thompson[,2]), stringsAsFactors = F)
Thompson.gi <- read.csv("../datafiles/2021.NC.SData5.csv", stringsAsFactors = F)
Thompson <- Thompson %>% mutate(A375.gi = ifelse(pair %in% Thompson.gi[,1], "Y", "N") ) %>% mutate(Mewo.gi = ifelse(pair %in% Thompson.gi[,2], "Y", "N") ) %>% mutate(RPE1.gi = ifelse(pair %in% Thompson.gi[,3], "Y", "N") )
Thompson <- Thompson %>% mutate("gi" = ifelse( A375.gi == "Y" | Mewo.gi == "Y" | RPE1.gi == "Y", "Negative", "others" ))


###Public data. Two in silico predicted pairs.
genename <- read.table("../datafiles/Homo_sapiens.GRCh37.73.genename", stringsAsFactors = F)

plosG <- read.table("../datafiles/journal.pgen.1008466.s005.csv", header=T, sep=",") 
for(i in 1:dim(plosG)[2]){if( class(plosG[,i]) == "factor" ){ plosG[,i] <- as.character(plosG[,i]) }}
plosG <- plosG %>% rename(gene1 = A1, gene2 = A2)
plosG <- plosG %>% mutate( pairname = ifelse(gene1 > gene2, paste(gene1, gene2, sep="."), paste(gene2, gene1, sep=".") ) )
plosG <- plosG[!(duplicated( plosG$pairname )), ] %>% mutate(pairname = NULL)

Hart  <- read.table("../datafiles/Hart_Supplementary_Table1.txt", header=T)
for(i in 1:dim(Hart)[2]){if( class(Hart[,i]) == "factor" ){ Hart[,i] <- as.character(Hart[,i]) }}
Hart  <- Hart %>% rename(gene1 = Gene_A, gene2 = Gene_B)
Hart  <- Hart %>% mutate( pairname = ifelse(gene1 > gene2, paste(gene1, gene2, sep="."), paste(gene2, gene1, sep=".") ) )
Hart  <- Hart[!(duplicated( Hart$pairname )), ] %>% mutate(pairname = NULL)

union.pair <- rbind.data.frame(plosG %>% select(gene1, gene2), Hart %>% select(gene1, gene2)) %>% unique()
for(i in 1:dim(union.pair)[2]){if( class(union.pair[,i]) == "factor" ){ union.pair[,i] <- as.character(union.pair[,i]) }}
union.pair <- union.pair %>% mutate( pairname = ifelse(gene1 > gene2, paste(gene1, gene2, sep="."), paste(gene2, gene1, sep=".") ) )
union.pair <- union.pair[!(duplicated( union.pair$pairname )), ] %>% mutate(pairname = NULL)



###GI and Genes
compare.gi.gene <- compare.gi  %>% select(gene1, gene2, gi) %>% melt(id.vars="gi") %>% select(gi, value) %>% unique()
Hart2.gene      <- Hart2       %>% select(gene1, gene2, gi) %>% melt(id.vars="gi") %>% select(gi, value) %>% unique()
Parrish.gene     <- Parrish      %>% select(gene1, gene2, gi) %>% melt(id.vars="gi") %>% select(gi, value) %>% unique()
nbt.gene        <- nbt         %>% select(gene1, gene2, gi) %>% melt(id.vars="gi") %>% select(gi, value) %>% unique()
Thompson.gene      <- Thompson     %>% select(gene1, gene2, gi) %>% melt(id.vars="gi") %>% select(gi, value) %>% unique()



###combing four dataset
Ours.pairs      <- compare.gi %>% mutate(pair1 = paste(gene1, gene2, sep="."), pair2 = paste(gene2, gene1, sep=".")) %>% select(pair1, pair2, gi, significant.x, significant.y) %>% melt(id.var=c("gi", "significant.x", "significant.y")) %>% mutate(variable = NULL) %>% rename(pair = value) %>% rename(HCT116.gi = significant.x, GM12878.gi = significant.y)
Hart2.pairs     <- Hart2      %>% mutate(pair1 = paste(gene1, gene2, sep="."), pair2 = paste(gene2, gene1, sep=".")) %>% select(pair1, pair2, gi, A549.gi, HT29.gi, OVCAR8.gi)  %>% melt(id.var=c("gi", "A549.gi", "HT29.gi", "OVCAR8.gi")) %>% mutate(variable = NULL) %>% rename(pair = value) %>% mutate(A549.gi = ifelse(A549.gi < -3, "Y", "N"), HT29.gi = ifelse(HT29.gi < -3, "Y", "N"), OVCAR8.gi = ifelse(OVCAR8.gi < -3, "Y", "N"))
Parrish.pairs    <- Parrish.raw %>% mutate(pair1 = paste(gene1, gene2, sep="."), pair2 = paste(gene2, gene1, sep=".")) %>% select(pair1, pair2, gi, significant.x, significant.y) %>% melt(id.var=c("gi", "significant.x", "significant.y")) %>% mutate(variable = NULL) %>% rename(pair = value) %>% rename(PC9.gi = significant.x, HeLa.gi = significant.y)
union.pair.list <- c(paste(union.pair$gene1, union.pair$gene2, sep='.'), paste(union.pair$gene2, union.pair$gene1, sep='.') )
Thompson.pairs    <- Thompson %>% mutate(pair1 = paste(gene1, gene2, sep="."), pair2 = paste(gene2, gene1, sep=".")) %>% select(pair1, pair2, gi, "A375.gi", "Mewo.gi", "RPE1.gi") %>% melt(id.var=c("gi", "A375.gi", "Mewo.gi", "RPE1.gi")) %>% mutate(variable = NULL) %>% rename(pair = value)


all.pair <- rbind.data.frame( compare.gi %>% select(gene1,gene2), compare.gi %>% select(gene2,gene1) %>% rename(x = gene2, gene2 = gene1) %>% rename(gene1 = x),  Hart2 %>% select(gene1,gene2), Hart2 %>% select(gene2,gene1) %>% rename(x = gene2, gene2 = gene1) %>% rename(gene1 = x), Parrish %>% select(gene1,gene2), Parrish %>% select(gene2,gene1) %>% rename(x = gene2, gene2 = gene1) %>% rename(gene1 = x), Thompson %>% select(gene1,gene2), Thompson %>% select(gene2,gene1) %>% rename(x = gene2, gene2 = gene1) %>% rename(gene1 = x) )  %>% unique()
all.pair <- all.pair %>% filter(gene1 > gene2)

all.pair.gi <- all.pair %>% mutate(pair = paste(gene1, gene2, sep='.')) %>% left_join(Ours.pairs, by="pair") %>% rename(Ours = gi) %>% left_join(Hart2.pairs, by="pair") %>% rename(Dede = gi) %>% left_join(Parrish.pairs, by="pair") %>% rename(Parrish = gi) %>% left_join(Thompson.pairs, by="pair") %>% rename(Thompson = gi)
all.pair.gi <- all.pair.gi %>% mutate(report = ifelse(pair %in% union.pair.list, "Reported", "novel") ) 



###Plots of negative GI pairs. Pairs are arranged according to the number of cell lines supporting negative GI.
all.pair.gi.count <- all.pair.gi %>% select(pair, HCT116.gi, GM12878.gi, A549.gi, HT29.gi, OVCAR8.gi, PC9.gi, HeLa.gi, A375.gi, Mewo.gi, RPE1.gi) %>% melt(id.var = "pair") %>% rename(gi = value) %>% select(pair, gi) %>% table() %>% as.data.frame() %>% dcast(pair~gi,value.var="Freq") %>% arrange(desc(Y),desc(N))
all.pair.gi.order <- rbind.data.frame(all.pair.gi.count %>% filter(Y > 0) %>% inner_join( all.pair.gi, by="pair") %>% arrange(desc(Y)), all.pair.gi.count %>% filter(Y == 0) %>% inner_join( all.pair.gi, by="pair") %>% arrange(desc(Ours), desc(Dede), desc(Parrish), desc(Thompson) ) ) %>% select(pair) %>% unlist() %>% as.character()

all.pair.gi.tile      <- all.pair.gi %>% select(pair, HCT116.gi, GM12878.gi, A549.gi, HT29.gi, OVCAR8.gi, PC9.gi, HeLa.gi, A375.gi, Mewo.gi, RPE1.gi) %>% melt(id.vars="pair") %>% na.omit() %>% mutate(gi=ifelse(value == "Y", "Negative", "others"))
all.pair.gi.tile$pair <-  factor(all.pair.gi.tile$pair, levels = all.pair.gi.order)

all.pair.gi.tile <- all.pair.gi.tile %>% inner_join(data.frame("pair" = all.pair.gi.order, "order" = 1:length(all.pair.gi.order) ), by="pair")
p.all.tile <- all.pair.gi.tile %>% ggplot() + geom_tile(aes(x=variable, y=order, fill=gi)) + theme_bw() + xlab("Cell lines") + ylab("Duplicate pairs") + scale_x_discrete(labels=c("HCT116 (Ours)", "GM12878 (Ours)", "A549 (Dede et al.)", "HT29 (Dede et al.)", "OVCAR8 (Dede et al.)", "PC9 (Parrish et al.)", "HeLa (Parrish et al.)", "A375 (Thompson et al.)", "Mewo (Thompson et al.)", "RPE1 (Thompson et al.)" )) + theme(axis.title=element_text(size=20, face="bold"), axis.text = element_text(size=17), legend.title=element_text(size=20, face="bold"), legend.text=element_text(size=18),  plot.tag=element_text(size=28, face="bold"), panel.background = element_rect(fill = NA), legend.position=c(0.8, 0.2) ) + scale_fill_manual(values=c("#F2CA66", "#3BA6D9"), labels=c("negative GI", "others"), guide=guide_legend( title = "GI class") ) + coord_flip() + scale_y_continuous(breaks=c(0:11)*200, labels=c(0:11)*200)


all.pair.gi.rate  <- all.pair.gi.count %>% mutate(cells = ifelse(N+Y <=3, "2-3 cell lines", "4-6 cell lines")) %>% mutate(cells = ifelse(N+Y >= 7, ">=7 cell lines", cells)) %>% mutate(gi = ifelse(Y > 0, "Negative", "others") ) %>% select(cells, gi) %>% table() %>% as.data.frame() %>% dcast(cells ~ gi, value.var="Freq") %>% mutate(sum = Negative + others) %>% mutate(ratio = Negative/sum) %>% mutate(sd = sqrt(ratio*(1-ratio)/sum)) %>% mutate(proportion = sum/nrow(all.pair.gi)) %>% mutate("Pair" = "Pair")
all.pair.gi.rate$cells <- factor(all.pair.gi.rate$cells, levels = c("2-3 cell lines" , "4-6 cell lines", ">=7 cell lines") )

p1 <- all.pair.gi.rate %>% ggplot() + geom_bar(aes(x=Pair, y=proportion, fill=cells), stat="identity") + geom_text(aes(x=Pair, y=y, label=sum), data=all.pair.gi.rate %>% mutate(y = ifelse(proportion < 0.5, proportion/2, 1-proportion/2 + 0.05)) %>% mutate(y = ifelse(y < 0.05, y, y+0.05) ), size=5.8) + scale_y_continuous(breaks=c(0:4 * 0.2), labels=c("0%", "20%", "40%", "60%", "80%")) + theme(axis.title.y=element_blank(), axis.title.x=element_text(size=20, face="bold"), axis.text = element_text(size=17), legend.title=element_text(size=20, face="bold"), legend.text=element_text(size=17), plot.tag=element_text(size=28, face="bold"), legend.position="top") + scale_fill_manual(values=c("#466E50", "#5D936A", "#AEC9B4"), labels=c("2-3 cell lines", "4-6 cell lines", ">=7 cell lines"), guide=guide_legend( title = "Assayed in", nrow=3) ) + ylab("Duplicate pairs")  + xlab("") + coord_polar("y")


k1  <-  fisher.test(all.pair.gi.rate %>% filter(cells != ">=7 cell lines") %>% select(Negative, others),alternative='less')$p.value %>% format(digits=2)
k2  <-  fisher.test(all.pair.gi.rate %>% filter(cells != "2-3 cell lines") %>% select(Negative, others),alternative='greater')$p.value %>% format(digits=2)
p2  <-  all.pair.gi.rate %>% ggplot(aes(x=cells, y=ratio)) + geom_bar(aes(fill=cells), stat="identity", alpha=0.9) + geom_errorbar(aes(x=cells, ymax=ratio + sd,ymin=ratio - sd), width=0.2) + theme_bw() + theme(axis.title=element_text(size=20, face="bold"), axis.text = element_text(size=17), legend.title=element_text(size=20, face="bold"), legend.text=element_text(size=18), legend.position="none", plot.tag=element_text(size=28, face="bold")) + geom_text(aes(x=cells, y=y, label=label), data=all.pair.gi.rate %>% mutate(y = 0.50, label=paste(Negative, sum, sep="/")), size=5.8) + scale_fill_manual(values=c("#466E50", "#5D936A", "#AEC9B4"), labels=c("2-3 cell lines", "4-6 cell lines", ">=7 cell lines", guide=guide_legend( title = "Assayed in") )) + xlab("Number of assayed cell lines") + ylab("The proportion of duplicate pairs\n with negative GI")
p2  <- p2 + geom_path(aes(x=x, y=y), data= data.frame("x"=c(1, 1, 2, 2), "y"=c(0.25, 0.26, 0.26, 0.25)), lwd=0.7) + geom_text(aes(x=x, y=y, label=label), data=data.frame("x"=1.5, "y"=0.28, label = k1 ), size=5.8) + geom_path(aes(x=x, y=y), data= data.frame("x"=c(2, 2, 3, 3), "y"=c(0.42, 0.43, 0.43, 0.42)), lwd=0.7) + geom_text(aes(x=x, y=y, label=label), data=data.frame("x"=2.5, "y"=0.45, label = k2 ), size=5.8)


all.pair.gi.bar       <- all.pair.gi.count %>% filter(Y > 0) %>% melt(id.var="pair")
all.pair.gi.bar.order <- all.pair.gi.count %>% filter(Y > 0) %>% select(pair) %>% unlist() %>% as.character()
all.pair.gi.bar$pair  <- factor(all.pair.gi.bar$pair, levels = all.pair.gi.bar.order )

all.pair.gi.bar       <- all.pair.gi.bar %>% inner_join(data.frame("pair" = all.pair.gi.bar.order, "order" = 1:length(all.pair.gi.bar.order) ), by="pair")
p.all.gi  <- ggplot(all.pair.gi.bar) + geom_bar(aes(x=order, y=value, fill=variable), stat="identity") + theme_bw() + xlab("Duplicate pairs") + ylab("Number of cell lines") + theme(axis.title=element_text(size=20, face="bold"), axis.text = element_text(size=17), legend.title=element_text(size=20, face="bold"), legend.text=element_text(size=18),  plot.tag=element_text(size=28, face="bold"), legend.position=c(0.8, 0.8) ) + scale_fill_manual(values=c("#3BA6D9", "#F2CA66"), labels=c("others", "negative GI"), guide=guide_legend( title = "GI class") ) + scale_y_continuous(breaks=0:10, labels=0:10) + scale_x_continuous(breaks=c(0:10)*20, labels=c(0:10)*20)

library(ggrepel)
p.all.gi <- p.all.gi + geom_text_repel(aes(x=order, y=value,label=pair), data= (all.pair.gi.bar %>% filter(order <= 10) %>% filter(variable == "Y")) )


all.pair.gi.ratio <- all.pair.gi.count %>% filter(Y > 0) %>% mutate(sum=N+Y) %>% mutate(ratio = Y/sum) %>% mutate(cells = ifelse(N+Y <=3, "2-3 cell lines", "4-6 cell lines")) %>% mutate(cells = ifelse(N+Y >= 7, ">=7 cell lines", cells)) %>% mutate(ratio = floor(ratio * 10 - 0.001) ) %>%  select(cells, ratio) %>% table() %>% as.data.frame()
all.pair.gi.ratio$cells <- factor(all.pair.gi.ratio$cells, levels = c("2-3 cell lines", "4-6 cell lines", ">=7 cell lines") )
all.pair.gi.half  <- all.pair.gi.count %>% filter(Y > 0) %>%  mutate(sum=N+Y) %>% mutate(cells = ifelse(N+Y <=3, "2-3 cell lines", "4-6 cell lines")) %>% mutate(cells = ifelse(N+Y >= 7, ">=7 cell lines", cells)) %>% mutate(half = ifelse(Y/sum > 0.5, "Major", "Minor")) %>% select(cells, half) %>% table() %>% as.data.frame() %>% dcast(cells~half, value.var="Freq") %>% mutate(sum = Minor + Major) %>% mutate(ratio = Major/sum)

p.all.ratio <- all.pair.gi.ratio %>% ggplot() + geom_bar(aes(x=ratio, y=Freq), stat="identity", width=0.7, alpha=0.8) + xlab( "Sharing frequency across cell lines" ) + ylab( "Number of negative GI pairs" ) +  facet_wrap(~cells, scale="free_y", nrow=3) + theme_bw() + theme(axis.title=element_text(size=20, face="bold"), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), legend.title=element_text(size=20, face="bold"), legend.text=element_text(size=17), strip.background = element_rect(fill = NA), strip.text=element_text(size=18), plot.tag=element_text(size=28, face="bold") ) + geom_vline(xintercept = 4.5, lty="dashed", col="blue") + scale_x_discrete(labels = c("(0,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.5]", "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]", "(0.9,1.0]") )
p.all.ratio <- p.all.ratio + geom_text(aes(x=x, y=y, label=label), data = all.pair.gi.half %>% mutate(x = 5.6, y = ifelse(cells == ">=7 cell lines", 7, 38), label = paste("(", format(ratio * 100, digits=3), "%", " > 0.5)", sep='') ), size=5.8)


all.pair.gi.report <- all.pair.gi %>% inner_join(all.pair.gi.count, by="pair") %>% filter(Y > 0) %>% mutate(sum=N+Y) %>% mutate(half = ifelse(Y/sum > 0.5, "Major", "Minor")) %>% select(half, report) %>% table() %>% as.data.frame() %>% dcast(half~report, value.var="Freq") %>% mutate(sum=novel+Reported) %>% mutate(ratio = Reported/sum) %>% mutate(sd = sqrt(ratio*(1-ratio)/sum))
all.pair.gi.report$half <- factor(all.pair.gi.report$half, levels=c("Minor", "Major") )
k <- fisher.test(all.pair.gi.report %>% select(Reported, novel),alternative = 'greater')$p.value %>% format(digits=2)

p.all.report  <- all.pair.gi.report %>% ggplot() + geom_bar(aes(x=half, y=ratio), stat="identity", alpha=0.9) + theme_bw() + xlab("Duplicate pairs") + ylab("Number of cell lines") + theme(axis.title=element_text(size=20, face="bold"), axis.text = element_text(size=17), legend.title=element_text(size=20, face="bold"), legend.text=element_text(size=18),  plot.tag=element_text(size=28, face="bold"), legend.position=c(0.8, 0.8) ) + geom_errorbar(aes(x=half, ymin=ratio-sd, ymax=ratio+sd), stat="identity", width=0.2) + geom_text(aes(x=half, y=y, label=label), data=all.pair.gi.report %>% mutate(y = 0.83, label=paste(Reported, sum, sep="/")), size=5.8) + xlab("Sharing frequency \nacross cell lines") + ylab("Negative GI pairs identified \nby the association dataset") + scale_x_discrete(labels = c("<=0.5", ">0.5") )
p.all.report  <- p.all.report + geom_path(aes(x=x, y=y), data= data.frame("x"=c(1, 1, 2, 2), "y"=c(0.68, 0.70, 0.70, 0.68)), lwd=0.7) + geom_text(aes(x=x, y=y, label=label), data=data.frame("x"=1.5, "y"=0.73, label = k ), size=5.8)



###Robust pairs
robust.pair.gi <- all.pair.gi %>% inner_join(all.pair.gi.count, by="pair") %>% filter(Y > 0) %>%  mutate(sum=N+Y) %>% mutate(cells = ifelse(N+Y <=3, "2-3 cell lines", "4-6 cell lines")) %>% mutate(cells = ifelse(N+Y >= 7, ">=7 cell lines", cells)) %>% mutate(half = ifelse(Y/sum > 0.5, "Major", "Minor")) %>% filter(half == "Major" | report == "Reported")

write.csv(all.pair.gi.count %>% inner_join(all.pair.gi, by=c("pair")) %>% mutate(sum=N+Y) %>% mutate(cells = ifelse(N+Y <=3, "2-3 cell lines", "4-6 cell lines")) %>% mutate(cells = ifelse(N+Y >= 7, ">=7 cell lines", cells)) %>% mutate(half = ifelse(Y/sum > 0.5, "Major", "Minor")) %>% mutate(half = ifelse(Y>0, half, "") ) %>% mutate(robust = ifelse(pair %in% robust.pair.gi$pair, "Robust", "")) %>% mutate(robust = ifelse(half == "Minor" & report == "novel", "Conditional", robust) ), "../TableS7_all.pair.gi.csv")



###WGD versus SSD.
ohnologs <- read.table("../datafiles/hsapiens.Pairs.Strict.2R.txt", header=T, sep="\t")
ohnologs <- ohnologs %>% rename(gene1 = Symbol1, gene2 = Symbol2)
ohnologs <- ohnologs %>% mutate(pair1 = paste(gene1, gene2, sep="."), pair2 = paste(gene2, gene1, sep="."))

all.pair.wgd  <- all.pair.gi.count %>% mutate(gi = ifelse(Y >= 1, "Negative", "others")) %>% mutate( "Duplication" = ifelse(pair %in% c(ohnologs$pair1, ohnologs$pair2), "WGD", "SSD")) %>% select(gi, Duplication) %>% table() %>% as.data.frame() %>% dcast(Duplication~gi, value.var="Freq") %>% mutate(sum = (Negative + others)) %>% mutate(ratio = Negative/sum) %>% mutate(sd = sqrt(ratio*(1-ratio)/sum))
all.pair.wgd$Duplication <- factor(all.pair.wgd$Duplication, levels=c("WGD", "SSD"))

all.pair.wgd.p.value <- fisher.test(all.pair.wgd %>% select(Negative, others),alternative = 'less')$p.value %>% format(digits = 2)
p.all.wgd.ssd <- all.pair.wgd %>% ggplot() + geom_bar(aes(x=Duplication, y=ratio, fill=Duplication), stat="identity") + ylim(c(0, 0.5))  + geom_errorbar(aes(x=Duplication, ymax=ratio + sd,ymin=ratio - sd), width=0.2) + theme_bw() + scale_fill_discrete( guide = guide_legend(title = "Type")) + xlab("Duplication type") + ylab ("Duplicate pair negative GI rate") + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position="none") + geom_text(aes(x=Duplication, y=value2, label=Freq2),size=5.8, data = all.pair.wgd %>% mutate(value2 = 0.5, Freq2=paste(Negative, sum, sep="/")) ) + geom_text(aes(x=Duplication, y=value2, label=label),size=5.8, data = all.pair.wgd %>% mutate(value2 = ratio + 0.03, label = format( ratio, digits = 1 )) ) + geom_path(aes(x=x, y=y), data= data.frame("x"=c(1, 1, 2, 2), "y"=c(0.38, 0.40, 0.40, 0.38)), lwd=0.7) + geom_text(aes(x=x, y=y, label=label), data=data.frame("x"=1.5, "y"=0.42, label = all.pair.wgd.p.value ), size=5.8)



###Reciprocal enrichment on public data.
enrichment.ratio <- function(dataset1, dataset2){
s <- rbind.data.frame(dataset1 %>% inner_join(dataset2 %>% select(gene1, gene2), by=c("gene1"="gene1", "gene2"="gene2")), dataset1 %>% inner_join(dataset2 %>% select(gene1, gene2), by=c("gene1"="gene2", "gene2"="gene1")) )  %>% unique() %>% select(gi) %>% table()
print(s)
}

enrichment.ratio(compare.gi, union.pair)
enrichment.ratio(Hart2, union.pair)
enrichment.ratio(nbt, union.pair)
enrichment.ratio(Parrish, union.pair)
enrichment.ratio(Thompson, union.pair)

#s0 <- data.frame("Study"=c("Our", "Hart", "nbt", "Parrish", "Thompson"), "Negative"=c(42,24,271,122,73), "others"=c(636,378,417,1030,645), "Type"="Expected" ) %>% mutate(ratio = Negative/(Negative + others)) %>% mutate(sd = sqrt(ratio*(1-ratio)/(Negative + others)))    ###This is the background level of several studies.
s0 <- data.frame("Study"=c("Yu", "Dede", "Parrish", "Thompson"), "Negative"=c(42,24,122,73), "others"=c(636,378,1030,645), "Type"="Expected" ) %>% mutate(ratio = Negative/(Negative + others)) %>% mutate(sd = sqrt(ratio*(1-ratio)/(Negative + others)))    ###This is the background level of several studies.

s3 <- data.frame("Study"=c("Yu", "Dede", "Parrish", "Thompson"), "Negative"=c(16,7,49,27), "others"=c(21,17,48,50), "Type"="Union" ) %>% mutate(ratio = Negative/(Negative + others)) %>% mutate(sd = sqrt(ratio*(1-ratio)/(Negative + others)))

s0[,1] <- factor(s0[,1], levels=c("Dede", "Parrish", "Thompson", "Yu") )
s3[,1] <- factor(s3[,1], levels=c("Dede", "Parrish", "Thompson", "Yu") )

p.enrichment <- ggplot(s3) + geom_bar(aes(x=Study, y=ratio, fill=Study), stat="identity", width=0.6) + geom_bar(aes(x=Study, y=ratio), stat="identity", lty="dashed", width=0.8, fill=NA, col="black", data=s0) + geom_errorbar(aes(x=Study, ymin=ratio - sd, ymax=ratio + sd), width=0.3) + geom_text(aes(x=Study, y=ratio+sd+0.05, label=paste(Negative,Negative + others,sep="/")), size=5.8, data = s3)  + geom_text(aes(x=Study, y=0.68, label=enrich), size=5.8, data = s3 %>% inner_join(s0, by="Study") %>% mutate(enrich = format(ratio.x/ratio.y, digits=2)) ) + geom_text(aes(x=x, y=y, label=label), size=5.8, data = data.frame("x"=0.6, "y"=0.68, "label"="OR"))+ theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), legend.position="none", plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold")) + xlab("Studies") + ylab("Negative GI rate of pairs \nfrom the association dataset") + scale_x_discrete(labels=c( "Dede \net al.", "Parrish \net al.", "Thompson \net al.", "Yu \net al."), position="top") + scale_fill_manual( values=c("#9D97C4", "#5585C5", "#C16852", "#FAD36A"))  # ggtitle("386 in silico predicted GI pairs from \n above two studies")


###Same Corum complexes
corum <- read.table("../datafiles/CORUM_allComplexes.txt",sep ="\t")
corum <- corum %>% inner_join(corum, by="V1")
names(corum) <- c("complex", "genename1", "genename2")
corum <- corum %>% filter(genename1 != genename2) %>% mutate(pairs = paste(genename1, genename2, sep="_"))

compare.gi$corum <- ifelse(paste(compare.gi$gene1, compare.gi$gene2, sep="_") %in% corum$pairs, "Same","Not.same")
Parrish$corum     <- ifelse(paste(Parrish$gene1, Parrish$gene2, sep="_") %in% corum$pairs, "Same","Not.same")
nbt$corum        <- ifelse(paste(nbt$gene1, nbt$gene2, sep="_") %in% corum$pairs, "Same","Not.same")
Thompson$corum        <- ifelse(paste(Thompson$gene1, Thompson$gene2, sep="_") %in% corum$pairs, "Same","Not.same")

complex.plot <- function(dataset){
s <- dataset %>% select(gi, corum) %>% table() %>% as.data.frame()
s <- s %>% inner_join( aggregate(s$Freq, by=list(s$corum), sum), by=c("corum" = "Group.1") ) %>% mutate(ratio = Freq/x) %>% mutate(position = ifelse(ratio > 0.5, ratio/2, 1-ratio/2))
s %>% dcast(gi~corum, value.var="Freq") %>% select(Same, Not.same) %>% fisher.test(alternative='greater') %>% print()  ###Test significance
k <- (fisher.test( s %>% dcast(gi~corum, value.var="Freq") %>% select(Same, Not.same),alternative='greater' ))$p.value %>% format(digits = 2)    ###P.value

p <- s %>% ggplot(aes(x=corum, y=ratio)) + geom_bar( aes(fill=gi), stat="identity") + geom_text( aes(x=corum, y=position, label=Freq), col="black", size=5.8 ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17),  plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold")) + xlab("Protein complex") + ylab("Percent of pairs") + scale_fill_manual(values=c("#F2CA66", "#3BA6D9"), labels=c("negative GI", "others"), guide=guide_legend( title = "GI class") ) + geom_text(aes(x=x, y=y, label=label), data = data.frame(x=1.5, y=1.03, label = paste("P=",k,sep='')), size=5.8)
return(p)
}

p.complex.Ours   <- complex.plot( compare.gi )
p.complex.Dede   <- complex.plot( Hart2 )
p.complex.Thomas <- complex.plot( nbt )
p.complex.Parrish <- complex.plot( Parrish )
p.complex.Thompson <- complex.plot( Thompson )

###Enrichment on small-sized protein family. From PANTHER Release 17.
PTHR_family <- read.csv( "../datafiles/PTHR17.0_human.csv", stringsAsFactors = F )
names(PTHR_family) <- c("genename", "subfamilyID", "familyID")

PTHR_family_size  <- PTHR_family %>% select(genename, familyID) %>% unique() %>% select(familyID) %>% table() %>% as.data.frame() %>% arrange( desc(Freq) ) %>% filter(Freq >= 2) %>% rename( familyID = ".", family_size = Freq)

PTHR_family_gene  <- PTHR_family %>% inner_join( PTHR_family_size, by="familyID" ) %>% arrange( family_size, familyID )

PTHR_family_all   <- PTHR_family_gene %>% select(genename, familyID, family_size) %>% inner_join(PTHR_family_gene %>% select(genename, familyID), by="familyID" ) %>% select(genename.x, genename.y, familyID, family_size) %>% filter( genename.x != genename.y ) %>% mutate(fam_type=ifelse(family_size<=3,'Small','Large'))


compare.gi.fam <- compare.gi %>% inner_join(PTHR_family_all,by=c('gene1'='genename.x','gene2'='genename.y'))
nbt <- nbt %>% inner_join(PTHR_family_all,by=c('gene1'='genename.x','gene2'='genename.y'))
Parrish.fam <- Parrish %>% inner_join(PTHR_family_all,by=c('gene1'='genename.x','gene2'='genename.y'))
Hart2.fam <- Hart2 %>% inner_join(PTHR_family_all,by=c('gene1'='genename.x','gene2'='genename.y'))
Thompson.fam <- Thompson %>% inner_join(PTHR_family_all,by=c('gene1'='genename.x','gene2'='genename.y'))


PTHR_family.plot <- function(dataset){
  p <- dataset %>% ggplot(aes(x=gi, y=family_size)) +geom_violin(aes(fill=gi),col='white')+ geom_boxplot( col='black', outlier.shape=NA, width=0.2 )  + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text (size=17), legend.position="none", plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold")) + xlab("GI class") + ylab("Gene family size") + scale_fill_manual(values=c("#F2CA66", "#3BA6D9"), labels=c("negative GI", "others"), guide=guide_legend( title = "GI class") )  + geom_signif(comparisons = list(c("Negative", "others")), textsize=5.8, y_pos=7) + scale_x_discrete(labels=c("Negative GI", "Non-nGI") ) + scale_y_continuous(trans = 'log2')
  return(p)
}

p.PTHR_family.Ours   <- PTHR_family.plot( compare.gi.fam )
p.PTHR_family.Dede   <- PTHR_family.plot( Hart2.fam )
p.PTHR_family.nbt <- PTHR_family.plot( nbt )
p.PTHR_family.Parrish <- PTHR_family.plot( Parrish.fam )
p.PTHR_family.Thompson <- PTHR_family.plot( Thompson.fam )

library(patchwork)
p.PTHR_family.Ours/p.PTHR_family.Dede/p.PTHR_family.Parrish/p.PTHR_family.Thompson


###Duplicate pair identities
ensembl <- read.table("../datafiles/Ensembl_Ver96_Paralogous_Identity.txt", header=T, sep="\t")

names(ensembl) <- c("query", "target", "query.identity", "target.identity")
ensembl <- ensembl[!(duplicated(ensembl[,1:2])), ]
ensembl <- ensembl %>% mutate(identity = (query.identity + target.identity)/2) %>% filter( identity > 20) %>% arrange(query, -identity)

compare.gi.identity  <- compare.gi  %>% inner_join(ensembl, by=c("gene1" = "query", "gene2" = "target"))
Hart2.identity       <- Hart2  %>% inner_join(ensembl, by=c("gene1" = "query", "gene2" = "target"))
nbt.identity         <- nbt    %>% inner_join(ensembl, by=c("gene1" = "query", "gene2" = "target"))
Parrish.identity     <- Parrish    %>% inner_join(ensembl, by=c("gene1" = "query", "gene2" = "target"))
Thompson.identity      <- Thompson    %>% inner_join(ensembl, by=c("gene1" = "query", "gene2" = "target"))

identity.plot <- function(dataset1){
p <- dataset1 %>% ggplot(aes(x=gi, y=identity)) + geom_boxplot( aes(col=gi), outlier.shape=NA, width=0.8 ) + geom_beeswarm(aes(col=gi), size=1)  + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text (size=17), legend.position="none", plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold")) + xlab("GI class") + ylab("Protein identity") + scale_color_manual(values=c("#F2CA66", "#3BA6D9"), labels=c("negative GI", "others"), guide=guide_legend( title = "GI class") )  + geom_signif(comparisons = list(c("Negative", "others")), textsize=5.8, y_pos=105) + scale_x_discrete(labels=c("negative GI", "others") ) 
return(p)
}

p.identity.Ours    <- identity.plot( compare.gi.identity )
p.identity.Dede    <- identity.plot( Hart2.identity )
p.identity.Thomas  <- identity.plot( nbt.identity )
p.identity.Parrish <- identity.plot( Parrish.identity )
p.identity.Thompson  <- identity.plot( Thompson.identity )
grid.arrange(p.identity.Ours, p.identity.Dede, p.identity.Thomas, p.identity.Parrish, p.identity.Thompson, layout_matrix = matrix( 1:5, nrow=1))


###Duplicate pair GO similarities
goSim <- read.csv("../datafiles/all_paralogs_GOsimilarity.csv", stringsAsFactors = F)

compare.gi.goSim  <- rbind.data.frame( compare.gi %>% inner_join(goSim, by=c("gene1" = "gene1", "gene2" = "gene2")), compare.gi %>% inner_join(goSim, by=c("gene1" = "gene2", "gene2" = "gene1")) )
Hart2.goSim       <- rbind.data.frame( Hart2  %>% inner_join(goSim, by=c("gene1" = "gene1", "gene2" = "gene2")), Hart2  %>% inner_join(goSim, by=c("gene1" = "gene2", "gene2" = "gene1")) )
nbt.goSim         <- rbind.data.frame( nbt    %>% inner_join(goSim, by=c("gene1" = "gene1", "gene2" = "gene2")), nbt    %>% inner_join(goSim, by=c("gene1" = "gene2", "gene2" = "gene1")) )
Parrish.goSim     <- rbind.data.frame( Parrish    %>% inner_join(goSim, by=c("gene1" = "gene1", "gene2" = "gene2")), Parrish    %>% inner_join(goSim, by=c("gene1" = "gene2", "gene2" = "gene1")) )
Thompson.goSim      <- rbind.data.frame( Thompson    %>% inner_join(goSim, by=c("gene1" = "gene1", "gene2" = "gene2")), Thompson    %>% inner_join(goSim, by=c("gene1" = "gene2", "gene2" = "gene1")) )


goSim.plot <- function(dataset1){
p <- dataset1 %>% ggplot(aes(x=gi, y=similarity.BMA)) + geom_boxplot( aes(col=gi), outlier.shape=NA, width=0.8 ) + geom_beeswarm(aes(col=gi), size=1)  + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text (size=17), legend.position="none", plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold")) + xlab("GI class") + ylab("GO similarity") + scale_color_manual(values=c("#F2CA66", "#3BA6D9"), labels=c("negative GI", "others"), guide=guide_legend( title = "GI class") )  + geom_signif(comparisons = list (c("Negative", "others")), textsize=5.8, y_pos=0.85) + scale_x_discrete(labels=c("negative GI", "others") )
return(p)
}

p.goSim.Ours   <- goSim.plot( compare.gi.goSim )
p.goSim.Dede   <- goSim.plot( Hart2.goSim )
p.goSim.Thomas <- goSim.plot( nbt.goSim )
p.goSim.Parrish <- goSim.plot( Parrish.goSim )
p.goSim.Thompson <- goSim.plot( Thompson.goSim )
grid.arrange(p.goSim.Ours, p.goSim.Dede, p.goSim.Thomas, p.goSim.Parrish, p.goSim.Thompson, layout_matrix = matrix( 1:5, nrow=1))


###Bioplex PPI partners. This data is used in later ANOVA analysis. Mind the pair name.
bioplex <- read.table("../datafiles/BioPlex_293T_Network_10K_Dec_2019.tsv", sep="\t", header=T, stringsAsFactors = F)

duplicate.pair <- rbind.data.frame(compare.gi %>% select(gene1, gene2), Hart2 %>% select(gene1, gene2), nbt %>% select(gene1, gene2), Parrish %>% select(gene1, gene2), Thompson %>% select(gene1, gene2)) %>% unique() %>% mutate(pair1 = paste(gene1, gene2, sep='.')) %>% mutate(pair2 = paste(gene2, gene1, sep='.'))
 
s <- bioplex %>% filter(SymbolA == "SRSF6" | SymbolB == "SRSF6") %>% select(SymbolA, SymbolB) %>% mutate("genename" = "SRSF6") %>% melt(id.var="genename") %>% filter(genename != value) %>% select(variable)
s <- s[as.numeric(0), ]
for(i in (c(duplicate.pair$gene1, duplicate.pair$gene2) %>% unique)){ k <- bioplex %>% filter(SymbolA == i | SymbolB == i) %>% select(SymbolA, SymbolB) %>% mutate("genename" = i) %>% melt(id.var="genename") %>% filter(genename != value)  %>% select(genename, value) %>% unique(); s <- rbind.data.frame(s, k)}

bioplex.size <- s %>% select(genename) %>% table() %>% as.data.frame()
names(bioplex.size) <- c("genename", "Size")

k <- rbind.data.frame(duplicate.pair %>% inner_join(s, by=c("gene1" = "genename")), duplicate.pair %>% inner_join(s, by=c("gene2" = "genename")) )
bioplex.pair1 <- k %>% select(pair1, value) %>% table() %>% as.data.frame() %>% filter(Freq > 0) %>% select(pair1, Freq) %>% table() %>% as.data.frame() %>% dcast(pair1~Freq, value.var="Freq.1") %>% rename( pair = pair1 )
bioplex.pair2 <- k %>% select(pair2, value) %>% table() %>% as.data.frame() %>% filter(Freq > 0) %>% select(pair2, Freq) %>% table() %>% as.data.frame() %>% dcast(pair2~Freq, value.var="Freq.1") %>% rename( pair = pair2 )
bioplex.pair  <- rbind.data.frame(bioplex.pair1, bioplex.pair2) %>% unique()
names(bioplex.pair)[2:3] <- c("Unique", "Shared")


bioplex.plot <- function(dataset1){
s  <- dataset1 %>% mutate(pair = paste(gene1, gene2, sep=".")) %>% inner_join(bioplex.pair, by=c("pair" = "pair")) %>% mutate(Ratio = Shared/(Shared + Unique)) %>% inner_join(bioplex.size, by=c("gene1" = "genename")) %>% inner_join(bioplex.size, by=c("gene2" = "genename"))

p  <- s %>% filter(Size.x >= 3, Size.y >= 3) %>% ggplot(aes(x=gi, y=Ratio)) + geom_boxplot( aes(col=gi), outlier.shape=NA, width=0.8 ) + geom_beeswarm(aes(col=gi), size=1) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text (size=17), legend.position="none", plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold")) + xlab("GI class") + ylab("PPI sharing propotion") + scale_color_manual(values=c("#F2CA66", "#3BA6D9"), labels=c("negative GI", "others"), guide=guide_legend( title = "GI class") )  + geom_signif(comparisons = list (c("Negative", "others")), textsize=5.8, y_pos=1.05) + scale_x_discrete(labels=c("negative GI", "others") )
return(p)
}

p.bioplex.Ours   <- bioplex.plot( compare.gi )
p.bioplex.Dede   <- bioplex.plot( Hart2 )
p.bioplex.Thomas <- bioplex.plot( nbt )
p.bioplex.Parrish  <- bioplex.plot( Parrish )
p.bioplex.Thompson <- bioplex.plot( Thompson )
grid.arrange(p.bioplex.Ours, p.bioplex.Dede, p.bioplex.Thomas, p.bioplex.Parrish, p.bioplex.Thompson, layout_matrix = matrix( 1:5, nrow=1))



### A comparison file of negative GI rate
s <- all.pair.gi %>% select(pair, A549.gi, HT29.gi, OVCAR8.gi, PC9.gi, HeLa.gi, A375.gi, Mewo.gi, RPE1.gi, HCT116.gi, GM12878.gi) %>% melt( id.vars="pair" ) %>% select(variable, value) %>% table() %>% as.data.frame() %>% dcast(variable~value, value.var="Freq") %>% mutate(Total = N + Y) %>% mutate(ratio = Y/Total) %>% mutate( sd = sqrt(ratio * (1-ratio)/Total))
s$CellLine <- sub(".gi", "", s$variable )
s$source <- c("Dede", "Dede", "Dede", "Parrish", "Parrish", "Thompson", "Thompson", "Thompson", "Ours", "Ours")

s$CellLine <- factor( s$CellLine, levels=s$CellLine )
s$source   <- factor( s$source, levels=c("Dede", "Parrish", "Thompson", "Ours") )

p.gi.rate <-  s %>% ggplot() + geom_bar(aes(x=CellLine, y=ratio, fill=source), stat="identity") + ylim(c(0, 0.4))  + geom_errorbar(aes(x=CellLine, ymax=ratio + sd,ymin=ratio - sd), width=0.2) + theme_bw() + scale_fill_manual( values=c("#9D97C4", "#5585C5", "#C16852", "#FAD36A"), guide = guide_legend(title = "Studies"), labels=c("Dede et al.", "Parrish et al.", "Thompson et al.", "Yu et al.") ) + xlab("Human cell lines") + ylab("The proportion of duplicate pairs\n with negative GI") + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17),  plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position=c(0.8,0.6)) + geom_text(aes(x=CellLine, y=value2, label=Freq2),size=5.8, data = s %>% mutate(value2 = 0.4, Freq2=paste(Y, Total, sep="/")) ) + geom_text(aes(x=CellLine, y=value2, label=label), size=5.8, data = s %>% mutate(value2 = ratio + 0.03, label = format( ratio, digits = 1 )) ) + geom_hline( yintercept=0.015, lty="dashed")



###Reciprocal enrichments. Odds ratios are calculated.
odds.ratio <- function(dataset1, dataset2){
s <- rbind.data.frame(dataset1 %>% inner_join(dataset2, by=c("gene1" = "gene1", "gene2" = "gene2")), dataset1 %>% inner_join(dataset2, by=c("gene1" = "gene2", "gene2" = "gene1"))) %>% select(gi.x, gi.y) %>% table() %>% as.matrix()
print(s)
print(paste("Odds ratio ", (s[1,1] * s[2,2]/(s[1,2] * s[2,1]) ) %>% format(digits=2), sep=": " ))
print(paste("P value ", fisher.test(s)$p.value %>% format(digits=2), sep=": " ))
}

odds.ratio(compare.gi, Hart2)
odds.ratio(compare.gi, nbt)
odds.ratio(compare.gi, Parrish)
odds.ratio(compare.gi, Thompson)
odds.ratio(Hart2, nbt)
odds.ratio(Hart2, Parrish)
odds.ratio(Hart2, Thompson)
odds.ratio(nbt, Parrish)
odds.ratio(nbt, Thompson)
odds.ratio(Parrish, Thompson)

odds.ratio(nbt %>% mutate(gi = ifelse(rpe1.negative.gi == "Negative", "Negative", "others")), Thompson %>% mutate(gi = ifelse(RPE1.gi == "Y", "Negative", "others")) )   #We specifically focused on RPE1 cell line data.


###Summary of feature enrichment
feature.table <- data.frame(data = c("Yu et al.", "Dede et al.", "Parrish et al.", "Thompson et al."), "complex" = c( "***", "N.S.", "***","***"),"Small family"=c("N.S.","**","N.S.","*"), "identity" = c( "**", "***", "***", "***"), "ppi" = c( "**", "**", "***", "***"), "goSim" = c( "***", "N.S.", "***", "***"), stringsAsFactors = F) %>% melt(id.var = "data") 
feature.table$data <- factor(feature.table$data, levels=c( "Dede et al.", "Parrish et al.", "Thompson et al.", "Yu et al.") %>% rev() )

p.feature <- feature.table %>% ggplot(aes(x=variable, y=data)) + geom_text(aes(label=value), size=5.8) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.title = element_text(size=20), axis.text = element_text(size=17), axis.ticks = element_line(size = 1), panel.grid = element_blank(), panel.border = element_rect(fill = NA, size=NA), plot.tag=element_text(size=28, face="bold")) + scale_x_discrete(labels=c("Complex","Family\nsize", "Identity", "PPI", "GO \nsimilarity"), position = "top") + xlab("Enrichment on features") + ylab("Studies")


OR.table <- data.frame(data = c( "Yu et al.", "Dede et al.", "Parrish et al.", "Thompson et al."), "Yu et al." = c(1.1, 1, 20.8, 28.9, 10.6), "Dede et al." = c(2.5, 20.8, 1, 18.5, 6.7), "Parrish et al." = c(1.7, 28.1, 18.5, 1, 7.4), "Thompson et al." = c(2.0, 10.6, 6.7, 7.4, 1), check.names = F ) %>% melt(id.var = "data")
OR.table$data <- factor(OR.table$data, levels=c("Dede et al.", "Parrish et al.", "Thompson et al.",  "Yu et al.") %>% rev() )
OR.table$variable <- factor(OR.table$variable, levels=c("Dede et al.", "Parrish et al.", "Thompson et al.", "Yu et al.") )

p.OR <- OR.table %>% filter(value != 1) %>% rename(`Odds ratio` = value) %>% ggplot(aes(x=variable, y=data)) + geom_tile(aes(fill=`Odds ratio`)) + geom_text(aes(label=`Odds ratio`), size=5.8) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text.x = element_text(size=17), axis.text.y = element_blank(), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_blank(),panel.grid = element_blank(), panel.border = element_rect(fill = NA, size=NA), legend.text=element_text(size=17), legend.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold")) +  scale_fill_viridis_c( option="inferno", alpha=0.6, direction = -1) + scale_x_discrete(labels=c("Dede \net al.", "Parrish \net al.", "Thompson\net al.", "Yu \net al."), position = "top") + xlab("Enrichment on other studies") + ylab("")


###Save above results
write.csv( rbind.data.frame(
nbt         %>% left_join(ensembl, by=c("gene1" = "query", "gene2" = "target")) %>% select(gene1, gene2, gi, corum, identity) %>% left_join( nbt.goSim %>% select(gene1, gene2, similarity.BMA), by=c("gene1", "gene2") )  %>% arrange(gi) %>% mutate( study = "G.P. et al.") ,
Hart2       %>% left_join(ensembl, by=c("gene1" = "query", "gene2" = "target")) %>% select(gene1, gene2, gi, corum, identity) %>% left_join( Hart2.goSim %>% select(gene1, gene2, similarity.BMA), by=c("gene1", "gene2") )  %>% arrange(gi) %>% mutate( study = "Dede et al.") ,
Parrish     %>% left_join(ensembl, by=c("gene1" = "query", "gene2" = "target")) %>% select(gene1, gene2, gi, corum, identity) %>% left_join( Parrish.goSim %>% select(gene1, gene2, similarity.BMA), by=c("gene1", "gene2") )  %>% arrange(gi) %>% mutate( study = "Parrish et al.") ,
Thompson    %>% left_join(ensembl, by=c("gene1" = "query", "gene2" = "target")) %>% select(gene1, gene2, gi, corum, identity) %>% left_join( Thompson.goSim %>% select(gene1, gene2, similarity.BMA), by=c("gene1", "gene2") )  %>% arrange(gi) %>% mutate( study = "Thompson et al.") ,
compare.gi  %>% left_join(ensembl, by=c("gene1" = "query", "gene2" = "target")) %>% select(gene1, gene2, gi, corum, identity) %>% left_join( compare.gi.goSim %>% select(gene1, gene2, similarity.BMA), by=c("gene1", "gene2") )  %>% arrange(gi) %>% mutate( study = "Yu et al.") ) %>% 
mutate(pair = paste(gene1, gene2, sep=".")) %>% mutate(report = ifelse(pair %in% union.pair.list, "Reported", "novel")) %>% left_join(bioplex.size, by= c("gene1" = "genename")) %>% left_join(bioplex.size, by=c("gene2" = "genename")) %>% left_join(bioplex.pair, by=c("pair" = "pair")) %>% mutate(Ratio = Shared/(Shared + Unique)) %>% mutate(Ratio = ifelse(Size.x >=3 & Size.y >=3, Ratio, "NA") ), "../TableS6_Studies.features.benchmark.csv" )



###Combine duplicate pairs.
robust.pair.gi.genes  <- robust.pair.gi %>% select(gene1, gene2, pair) %>% melt(id.var="pair") %>% select(value) %>% unlist() %>% as.character() %>% unique()
all.pair.gi.genes     <- all.pair.gi.count %>% filter(Y > 0) %>% inner_join(all.pair.gi, by="pair") %>% select(gene1, gene2, pair) %>% melt(id.var="pair") %>% select(value) %>% unlist() %>% as.character() %>% unique()
condition.pair.gi.genes  <- all.pair.gi.genes %>% setdiff(robust.pair.gi.genes)
other.pair.gi.genes   <- all.pair.gi.count %>% filter(Y == 0) %>% inner_join(all.pair.gi, by="pair") %>% select(gene1, gene2, pair) %>% melt(id.var="pair") %>% select(value) %>% unlist() %>% as.character() %>% unique()
other.pair.gi.genes   <- other.pair.gi.genes %>% setdiff(all.pair.gi.genes)

combine.gene.all <- rbind.data.frame( data.frame("genename" = robust.pair.gi.genes, "type" = "Robust"), data.frame("genename" = condition.pair.gi.genes, "type" = "Conditional"),  data.frame("genename" = other.pair.gi.genes, "type" = "Others") )
combine.gene.all$type   <- factor(combine.gene.all$type,  levels=c("Robust", "Conditional", "Others"))


robust.gi.pairs     <- robust.pair.gi %>% select(gene1, gene2, pair)
all.gi.pairs        <- all.pair.gi.count %>% filter(Y > 0) %>% inner_join(all.pair.gi, by="pair") %>% select(gene1, gene2, pair)
condition.gi.pairs  <- all.gi.pairs      %>% anti_join(robust.gi.pairs, by=c("pair"))
other.gi.pairs      <- all.pair.gi.count %>% filter(Y == 0) %>% inner_join(all.pair.gi, by="pair")%>% select(gene1, gene2, pair)

combine.pair.all <- rbind.data.frame(robust.gi.pairs %>% mutate("type" = "Robust"), condition.gi.pairs %>% mutate("type" = "Conditional"), other.gi.pairs %>% mutate("type" = "Others") )
combine.pair.all$type   <- factor(combine.pair.all$type,  levels=c("Robust", "Conditional", "Others"))
combine.pair.all$corum  <- ifelse(paste(combine.pair.all$gene1, combine.pair.all$gene2, sep="_") %in% corum$pairs, "Same", "Not.same")

a <- combine.pair.all %>% inner_join(PTHR_family_all,by=c('gene1'='genename.x','gene2'='genename.y')) 
gi_famsize <- a%>% select(type,family_size) %>%  ggplot(aes(x=type,y=family_size))+ geom_violin(aes(fill=type),col=NA) + geom_boxplot(fill="white",width=0.15,outlier.shape = NA,col='black') +theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17),  plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"))+ geom_signif(comparisons = list(c("Robust", "Conditional"), c("Conditional","Others"),c("Robust", "Others")), textsize=5.8, y_position=c(7.5, 7.5,8),test = 'wilcox.test',test.args = 'less') + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Non-nGI" ) ) + ylab('Family size')+xlab("Negative GI type")+scale_y_continuous(trans = 'log2',breaks = c(2,4,8,32,128,512)) +scale_x_discrete(labels=c('Core','Conditional','Non-nGI'))+theme(legend.position = 'none')

###Coverage of gene family. From Ensembl Ver96. From Ensembl Ver102, No such information is provided.
ensembl.family <- read.csv("../datafiles/Ensembl.Biomart.Family.V96.txt", header=T, sep="\t", stringsAsFactors = F)
chromosome <- c(1:22, "X", "Y")
ensembl.family <- ensembl.family %>% filter( Chromosome.scaffold.name %in% chromosome ) %>%  select(Gene.stable.ID, Gene.name, Ensembl.Protein.Family.ID.s. ) %>% na.omit() %>% filter(Ensembl.Protein.Family.ID.s. != "") %>% unique()
names(ensembl.family) <- c("geneID", "genename", "subfamilyID")

ensembl.family$familyID <- sub("_.*$", "", ensembl.family$subfamilyID )
#PTHR.family <- ensembl.family %>% select(familyID) %>% table() %>% as.data.frame() %>% filter(Freq >= 2)
PTHR.family <- ensembl.family %>% select(geneID, genename, familyID) %>% unique()  %>% select(familyID) %>% table() %>% as.data.frame() %>% filter(Freq >= 2, Freq <= 100)
names(PTHR.family)[1] <- "familyID"

all.pair.PTHR.family  <- all.pair.gi %>% select(pair, gene1, gene2) %>% inner_join( ensembl.family %>% select(genename, familyID), by=c("gene1" = "genename") ) %>% inner_join( ensembl.family %>% select(genename, familyID), by=c("gene2" = "genename") ) %>% filter( familyID.x == familyID.y) %>% unique()

a = PTHR.family %>% select( familyID ) %>% unique() %>% nrow()
b = all.pair.PTHR.family %>% select( familyID.x ) %>% unique() %>% nrow()
c = a - b
all.pair.PTHR.covered  <- data.frame( "Type"=c("All", "Covered", "Not covered"), "Frequency"=c(a, b, c) )
all.pair.PTHR.covered$Type <- factor( all.pair.PTHR.covered$Type, levels=c("Not covered", "Covered" ) ) 

p.PTHR.covred <- all.pair.PTHR.covered %>% filter(Type != "All") %>% ggplot(aes(x=factor(1), y=Frequency)) + geom_bar(aes(fill=Type), stat="identity") + geom_text( aes(x=factor(1), y=Frequency, label=Frequency), data=all.pair.PTHR.covered %>% filter(Type != "All"), size=5.8) + coord_polar("y", direction = 1) + theme_bw() + theme(title = element_text(size=20), axis.title = element_blank(), axis.text = element_text(size=17), axis.ticks = element_line(size=rel(1)), legend.title=element_text(size=20, face="bold"), legend.text=element_text(size=18), plot.tag=element_text(size=28, face="bold")) + ggtitle("PANTHER family")


###Coverage of duplicate pair.
a <- nrow(ensembl)/2
b <- combine.pair.all %>% inner_join( ensembl, by=c("gene1" = "query", "gene2" = "target")) %>% nrow()
c <- a - b
all.pair.Ensembl.covered  <- data.frame( "Type"=c("All", "Covered", "Not covered"), "Frequency"=c(a, b, c) )
all.pair.Ensembl.covered$Type <- factor( all.pair.Ensembl.covered$Type, levels=c("Not covered", "Covered" ) ) 

p.Ensembl.covred <- all.pair.Ensembl.covered %>% filter(Type != "All") %>% ggplot(aes(x=factor(1), y=Frequency)) + geom_bar(aes(fill=Type), stat="identity") + geom_text( aes(x=factor(1), y=Frequency, label=Frequency), data=all.pair.Ensembl.covered %>% filter(Type != "All"), size=5.8) + coord_polar("y", direction = 1) + theme_bw() + theme(title = element_text(size=20), axis.title = element_blank(), axis.text = element_text(size=17), axis.ticks = element_line(size=rel(1)), legend.title=element_text(size=20, face="bold"), legend.text=element_text(size=18), plot.tag=element_text(size=28, face="bold")) + ggtitle("All paralogous pairs") 



###A comparison of negative GI rate in Thompson et al.
Thompson <- read.csv("../datafiles/2021.NC.SData2.csv", stringsAsFactors = F)
Thompson.gi <- read.csv("../datafiles/2021.NC.SData5.csv", stringsAsFactors = F)

Thompson <- data.frame("pair" = Thompson[,2], "gene1" = sub("_.*", "", Thompson[,2]), "gene2" = sub(".*_", "", Thompson[,2]), "Type" = Thompson[,3] )
Thompson <- Thompson %>% mutate( Type = ifelse( Type == "Paralogous_gene_pair", "Duplicate", "Non-duplicate") )
Thompson <- Thompson %>% mutate(A375.gi = ifelse(pair %in% Thompson.gi[,1], "Y", "N") ) %>% mutate(Mewo.gi = ifelse(pair %in% Thompson.gi[,2], "Y", "N") ) %>% mutate(RPE1.gi = ifelse(pair %in% Thompson.gi[,3], "Y", "N") )

Thompson.frequency <- Thompson %>% select(pair, A375.gi, Mewo.gi, RPE1.gi) %>% melt( id.var = "pair") %>% select( pair, value ) %>% table() %>% as.data.frame() %>% dcast(pair~value, value.var="Freq") %>% arrange(desc(Y))

Thompson <- Thompson %>% inner_join( Thompson.frequency, by="pair" ) %>% arrange(desc(Y), Type) %>% mutate( sharing = ifelse(Y == 1, "1/3", "") ) %>% mutate( sharing = ifelse(Y > 1, "2/3 or 3/3", sharing) ) %>% select(gene1, gene2, pair, Type, Y, N, sharing, A375.gi, Mewo.gi, RPE1.gi) 
write.csv(Thompson, "../TableS7_inset_Thompson_data_final.csv")


Thompson.rate <- Thompson %>% select(pair, Type, A375.gi, Mewo.gi, RPE1.gi) %>% melt(id.var = c("Type", "pair")) %>% select(Type, variable, value) %>% table() %>% as.data.frame() %>% dcast(Type+variable~value, value.var = "Freq") %>% mutate(ratio = Y/(Y + N))
k <- t.test(Thompson.rate[Thompson.rate$Type == "Duplicate",]$ratio, Thompson.rate[Thompson.rate$Type == "Non-duplicate",]$ratio )$p.value %>% format(digits = 2)
s <- data.frame("source" = c("Duplicate", "Non-duplicate"), "value" = c( mean( Thompson.rate[Thompson.rate$Type == "Duplicate",]$ratio ), mean( Thompson.rate[Thompson.rate$Type == "Non-duplicate",]$ratio ) ), "se" = c( sd( Thompson.rate[Thompson.rate$Type == "Duplicate",]$ratio ), sd( Thompson.rate[Thompson.rate$Type == "Non-duplicate",]$ratio ) ), source = c("Thompson", "Thompson") )

p.gi.rate.Thompson <- s %>% ggplot() + geom_bar(aes(x=source, y=value, fill=source), stat="identity") + geom_errorbar( aes(x=source, ymin=value-se, ymax=value+se), width=0.3 ) + ylim(c(0, 0.4))  + theme_bw() + scale_fill_manual( values=c("#5585C5", "gray70"), guide = guide_legend(title = "Studies"), labels=c("Thompson et al.") ) + xlab("Gene pairs") + ylab ("The proportion of gene pairs\n with negative GI") + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.title=element_text(size=21, face="bold"), legend.text=element_text(size=17),  plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position=c(0.5,0.7)) + geom_text(aes(x=source, y=value2, label=label), size=5.8, data = s %>% mutate(value2 = value + 0.03, label = format( value, digits = 2 )) ) + geom_path(aes(x=x, y=y), data= data.frame("x"=c(1, 1, 2, 2), "y"=c(0.14, 0.15, 0.15, 0.14)), lwd=0.7) + geom_text(aes(x=x, y=y, label=label), data=data.frame(x=1.5, y=0.17, label=k), size=5.8)


###A comparison of negative GI sharing frequency in Thompson et al.
Thompson.sharing <- Thompson %>% filter(Y > 0) %>% mutate( Cell = ifelse(sharing == "1/3", "One", "Two") ) %>% select( Type, Cell) %>% table() %>% as.data.frame() %>% dcast(Type~Cell, value.var="Freq") %>% mutate( Sum = One + Two ) %>% melt(id.var = c("Type", "Sum")) %>% mutate(ratio = value/Sum)
c <- (Thompson.sharing %>% select( Type, variable, value) %>% dcast(Type ~ variable) %>% mutate( Type = NULL) %>% fisher.test(alternative = 'less') ) $p.value %>% format( digits = 2 )

p.gi.sharing.Thompson <- Thompson.sharing %>% ggplot() + geom_bar(aes(x=variable, y=ratio, fill=Type), stat="identity", position = "dodge")  + theme_bw() + scale_fill_manual( values=c( "#5585C5", "gray70"), guide = guide_legend(title = "Gene pairs") ) + ylab("Proportion") + xlab ("Sharing frequency \nin Thompson et al.") + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.title=element_text(size=21, face="bold"), legend.text=element_text(size=17),  plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position=c(0.7, 0.7) ) + scale_x_discrete(labels=c("1/3", "2/3 or 3/3")) + scale_y_continuous( breaks = (0:5)/5, labels= (0:5)/5) + geom_path(aes(x=x, y=y), data= data.frame("x"=c(1, 1, 2, 2), "y"=c(0.90, 0.92, 0.92, 0.90)), lwd=0.7) + geom_text(aes(x=x, y=y, label=label), data=data.frame(x=1.5, y=0.95, label=c), size=5.8)


###An analysis on Horbeck et al.(2018) data.
Cell_data <- read.csv( "../datafiles/2018.GeneticLandscape.Cell.Table_S5.csv" )
Cell_data_K562   <- Cell_data %>% filter(gene1 != gene2) %>% select(gene1, gene2, K562.GI  ) %>% na.omit() %>% mutate( gi = ifelse(K562.GI < -3, "Y", "N") )
Cell_data_Jurkat <- Cell_data %>% filter(gene1 != gene2) %>% select(gene1, gene2, Jurkat.GI) %>% na.omit() %>% mutate( gi = ifelse(Jurkat.GI < -3, "Y", "N") )
( Cell_data_K562   %>% filter( gi=="Y" ) %>% nrow() )/( Cell_data_K562   %>% nrow() )
( Cell_data_Jurkat %>% filter( gi=="Y" ) %>% nrow() )/( Cell_data_Jurkat %>% nrow() )

Cell_data_cells <- Cell_data_K562 %>% inner_join( Cell_data_Jurkat, by=c("gene1", "gene2") )
Cell_data_cells %>% select(gi.x, gi.y) %>% table()
( Cell_data_cells %>% filter(gi.x == "Y", gi.y == "Y") %>% nrow() )/( Cell_data_cells %>% filter(gi.x == "Y" | gi.y == "Y") %>% nrow() )


###Comparing sharing frequency between duplicate pairs in combined set and non-duplicate pairs in Horbeck et al. (2018).
fisher.test( matrix(c(50, 160, 82, 272 + 850), nrow=2)  )


###Compositary plots

pdf("../Figure3_all.pair.gi_Rate_2021Oct.pdf", width=18, height=18)
grid.arrange( p.gi.rate + labs(tag="A"), p.gi.rate.Thompson + labs(tag = "B"), p2 + labs(tag="C"), p.all.ratio + labs(tag="D") + theme(axis.text.x  = element_text(angle = 20)), p.gi.sharing.Thompson + labs(tag = "E"), p.all.report + labs(tag="F") + xlab("Sharing frequency \nacross cell lines") , ggplot() + theme_minimal(), layout_matrix = matrix(c(rep(1,6), rep(2,2), rep(3,3), rep(4,5), rep(5,3), rep(6,2), rep(7,3) ), nrow=3, byrow=T) )

dev.off()

pdf("../Figure2_Compare.gi.GI_Enrichment_2021Oct.pdf", width=20, height=10)
grid.arrange( p.enrichment+ labs(tag="A"), p.feature + labs(tag="B"), p.OR + labs(tag=""), p.Ensembl.covred + labs(tag="C"), p.PTHR.covred + labs(tag="D"), p1 + labs(tag="E") , layout_matrix = matrix(c( rep(1,6), rep(2,6), rep(3,6), rep(4,5), rep(5,5), rep(6,6), rep(7,2) ), nrow=2, byrow=T) )

dev.off()



pdf("../FigureS2_Compare.gi.features.addNC_2021Oct.pdf", width=29, height=23)
grid.arrange( 
p.complex.Dede + labs(tag=""),p.PTHR_family.Dede + labs(tag=""), p.identity.Dede + labs(tag=""), p.bioplex.Dede + labs(tag=""), p.goSim.Dede + labs(tag=""), 
p.complex.Parrish + labs(tag=""), p.PTHR_family.Parrish + labs(tag=""),p.identity.Parrish + labs(tag=""), p.bioplex.Parrish + labs(tag=""), p.goSim.Parrish + labs(tag=""), 
p.complex.Thompson + labs(tag=""), p.PTHR_family.Thompson + labs(tag=""), p.identity.Thompson + labs(tag=""), p.bioplex.Thompson + labs(tag=""), p.goSim.Thompson + labs(tag=""), 
p.complex.Ours + labs(tag=""),p.PTHR_family.Ours + labs(tag=""), p.identity.Ours + labs(tag=""), p.bioplex.Ours + labs(tag=""), p.goSim.Ours + labs(tag=""), layout_matrix = matrix(c(1:20), nrow=4, byrow=T) )
dev.off()

pdf('../FigureS3_Compare_cellline_familysize_2022Jun.pdf',width=8,height=6)
gi_famsize
dev.off()###Exploring function of neagative GI dupliate genes and robust negative GI duplicate pairs.

p.combine.identity <- combine.pair.all %>% inner_join(ensembl, by=c("gene1" = "query", "gene2" = "target") ) %>% ggplot(aes(x=type, y=identity)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(101, 105)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Non-nGI" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Protein identity") + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI")) + ylim(c(20, 110))

p.combine.bioplex  <- combine.pair.all %>% inner_join(bioplex.pair, by=c("pair" = "pair")) %>% mutate(Ratio = Shared/(Shared + Unique)) %>% inner_join(bioplex.size, by= c("gene1" = "genename")) %>% inner_join(bioplex.size, by=c("gene2" = "genename")) %>% filter(Size.x >=3, Size.y >=3) %>% ggplot(aes(x=type, y=Ratio)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(0.95, 1.04)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9")) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("PPI partner sharing ratio") + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI")) + ylim(c(0, 1.1))

p.combine.shared   <- combine.pair.all %>% inner_join(bioplex.pair, by=c("pair" = "pair")) %>% mutate(Ratio = Shared/(Shared + Unique)) %>% inner_join(bioplex.size, by= c("gene1" = "genename")) %>% inner_join(bioplex.size, by=c("gene2" = "genename")) %>% filter(Size.x >=3, Size.y >=3) %>% ggplot(aes(x=type, y=Shared)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(105, 115)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9")) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Number of shared PPI partners") + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))

p.combine.PPI      <- combine.pair.all %>% inner_join(bioplex.size, by=c("gene1"="genename")) %>% inner_join(bioplex.size, by=c("gene2"="genename")) %>% mutate(Size.sum = (Size.x + Size.y)) %>% ggplot(aes(x=type, y=Size.sum)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(2.58, 2.73)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9")) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Number of total PPI partners") + scale_y_log10() + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))

k1 <- (combine.pair.all %>% select(type, corum) %>% table() %>% as.data.frame() %>% dcast(type ~ corum, value.var="Freq")  %>% filter(type %in% c("Robust", "Conditional") ) %>% mutate(type = NULL) %>% as.matrix() %>% fisher.test(alternative = 'less'))$p.value %>% format(digits=2)
k2 <- (combine.pair.all %>% select(type, corum) %>% table() %>% as.data.frame() %>% dcast(type ~ corum, value.var="Freq")  %>% filter(type %in% c("Robust", "Others") ) %>% mutate(type = NULL) %>% as.matrix() %>% fisher.test(alternative='less'))$p.value %>% format(digits=2)
p.combine.complex  <- combine.pair.all %>% select(type, corum) %>% table() %>% as.data.frame() %>% dcast(type ~ corum, value.var="Freq") %>% mutate(ratio = Same/(Not.same + Same)) %>% mutate(sd = sqrt(ratio*(1 - ratio)/(Same + Not.same)) ) %>% ggplot(aes(x=type,  y=ratio)) + geom_bar(aes(fill = type), stat="identity") + geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), width=0.3) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9")) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Probability in same complex") + geom_text(aes(x=x, y=y, label=label), data=data.frame(x=1.5, y=0.32, label=k1), size=5.8) + geom_path(aes(x=x, y=y), data= data.frame("x"=c(1, 1, 2, 2), "y"=c(0.30, 0.31, 0.31, 0.30)) ) + geom_text(aes(x=x, y=y, label=label), data=data.frame(x=2.0, y=0.35, label=k2), size=5.8) + geom_path(aes(x=x, y=y), data= data.frame("x"=c(1, 1, 3, 3), "y"=c(0.33, 0.34, 0.34, 0.33)) )



###Total GO annotation number 
Go_sum <- read.csv("../datafiles/all_paralogous_GOsum.csv", stringsAsFactors = F)

p.GO.BP <- combine.pair.all %>% inner_join(Go_sum %>% mutate(type = NULL), by=c("gene1" = "gene1", "gene2" = "gene2") ) %>% mutate( BP_sum = log10(BP_sum + 1), MF_sum = log10(MF_sum + 1), CC_sum = log10(CC_sum + 1), ) %>%  ggplot(aes(x=type, y=BP_sum)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(2.2, 2.5)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Number of GO terms\n in biological processes")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI")) + scale_y_continuous(breaks=c(0, 1, 2), labels=c(0, 10, 100))


p.GO.CC <- combine.pair.all %>% inner_join(Go_sum %>% mutate(type = NULL), by=c("gene1" = "gene1", "gene2" = "gene2") ) %>% mutate( BP_sum = log10(BP_sum + 1), MF_sum = log10(MF_sum + 1), CC_sum = log10(CC_sum + 1), ) %>%  ggplot(aes(x=type, y=CC_sum)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(2.2, 2.5)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Number of GO terms\n in cellular components") + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI")) + scale_y_continuous(breaks=c(0, 1, 2), labels=c(0, 10, 100))


###Total expression level. Prenatal development
###Transcriptome Tau index. Prenatal development

gene.expression.prenatal  <- readRDS( "../datafiles/gene.expression.prenatal.rds" )
gene.tau.prenatal  <- readRDS( "../datafiles/gene.tau.prenatal.rds")
gene.tau.prenatal.tissue  <- readRDS( "../datafiles/gene.tau.prenatal.tissue.rds" )
 
p.combine.expression.level <- combine.pair.all %>% inner_join( gene.expression.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.expression.prenatal, by=c("gene2" = "genename") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2) %>% mutate(TPM = log2(TPM + 1)) %>% ggplot(aes(x=type, y=TPM)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(14, 16)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression level")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))

p.combine.expression.tau <- combine.pair.all %>% inner_join( gene.tau.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.tau.prenatal, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% ggplot(aes(x=type, y=tau)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(0.9, 1.0)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression specificity")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))  + scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), labels = c(0.2, 0.4, 0.6, 0.8, 1.0))

p.combine.expression.tau.tissue <- combine.pair.all %>% inner_join( gene.tau.prenatal.tissue, by=c("gene1" = "genename") ) %>% inner_join( gene.tau.prenatal.tissue, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% ggplot(aes(x=type, y=tau)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(1.0, 1.1)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression specificity \n(organ level)")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))



###Chicken duplicate pairs. Protein Identity

ensembl.chicken <- read.table("../datafiles/Chicken_Ensembl_Ver102_Paralogous_Identity.txt", header=T, sep="\t", stringsAsFactors=F)
ensembl.chicken <- ensembl.chicken[,-c(3:4)] %>% na.omit()
names(ensembl.chicken) <- c("query", "target", "query.identity", "target.identity")
ensembl.chicken <- ensembl.chicken[!(duplicated(ensembl.chicken[,1:2])), ]
ensembl.chicken <- ensembl.chicken %>% mutate(identity = (query.identity + target.identity)/2) %>% filter( identity > 20) %>% arrange(query, -identity)

hcop.chicken <- read.table("../datafiles/human_chicken_hcop_fifteen_column.parsed.txt", sep="\t", stringsAsFactors=F)
names( hcop.chicken ) <- c("human", "chicken")
combine.pair.all.chicken <- combine.pair.all %>% inner_join( hcop.chicken, by=c("gene1" = "human") ) %>% inner_join( hcop.chicken, by=c("gene2" = "human") )

p.combine.identity.chicken <- combine.pair.all.chicken %>% inner_join(ensembl.chicken, by=c("chicken.x" = "query", "chicken.y" = "target") ) %>% ggplot(aes(x=type, y=identity)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(101, 105)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Protein identity")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI")) + ylim(c(20, 110))


###Chicken duplicate pairs. Expression level

chicken_gene.expression.prenatal <- readRDS( "../datafiles/chicken.gene.expression.prenatal.rds" )
chicken_gene.tau.prenatal  <- readRDS( "../datafiles/chicken.gene.tau.prenatal.rds")
chicken_gene.tau.prenatal.tissue  <- readRDS( "../datafiles/chicken.gene.tau.prenatal.tissue.rds")

p.combine.expression.level.chicken <- combine.pair.all.chicken %>% inner_join( chicken_gene.expression.prenatal, by=c("chicken.x" = "geneID") ) %>% inner_join( chicken_gene.expression.prenatal, by=c("chicken.y" = "geneID") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2) %>% mutate(TPM = log2(TPM + 1)) %>% ggplot(aes(x=type, y=TPM)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(14, 16)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression level")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))

p.combine.expression.tau.chicken  <- combine.pair.all.chicken %>% inner_join( chicken_gene.tau.prenatal, by=c("chicken.x" = "geneID") ) %>% inner_join( chicken_gene.tau.prenatal, by=c("chicken.y" = "geneID") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% ggplot(aes(x=type, y=tau)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(0.9, 1.0)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression specificity")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))  + scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), labels = c(0.2, 0.4, 0.6, 0.8, 1.0))

p.combine.expression.tau.tissue.chicken  <- combine.pair.all.chicken %>% inner_join( chicken_gene.tau.prenatal.tissue, by=c("chicken.x" = "geneID") ) %>% inner_join( chicken_gene.tau.prenatal.tissue, by=c("chicken.y" = "geneID") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% ggplot(aes(x=type, y=tau)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(1.0, 1.1)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression specificity \n(organ level)")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))



###Gentree age dating. Gene branch. Lower branch of a pair
gentree <- read.csv("../datafiles/hg19_ver73_age.tsv", sep="\t", stringsAsFactors = T)
taxon <- c( "Euteleostomi \nor before", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Boreoeutheria", "Euarchontoglires", "Simiiformes", "Catarrhini", "Hominoidea", "Hominidae", "Homo sapiens &\n Pan troglodytes", "Homo sapiens" )

combine.pair.all.branch <- combine.pair.all %>% inner_join(genename, by=c("gene1" = "V2")) %>% inner_join(gentree, by=c("V1" = "gene")) %>% mutate(V1 = NULL) %>% inner_join(genename, by=c("gene2" = "V2")) %>% inner_join(gentree, by=c("V1" = "gene")) %>% mutate(V1 = NULL) %>% mutate(branch = ifelse(branch.x < branch.y, branch.y, branch.x) ) %>% select( type, branch) %>% table() %>% as.data.frame()

k1 <- (combine.pair.all.branch %>% filter(branch %in% c(0, 1, 2, 3))  %>% dcast(type ~ branch, value.var="Freq") %>% filter(type %in% c("Robust", "Conditional") ) %>% mutate(type = NULL) %>% as.matrix() %>% fisher.test(alternative='greater'))$p.value %>% format(digits = 2)
k2 <- (combine.pair.all.branch %>% filter(branch %in% c(0, 1, 2, 3))  %>% dcast(type ~ branch, value.var="Freq") %>% filter(type %in% c("Robust", "Others") ) %>% mutate(type = NULL) %>% as.matrix() %>% fisher.test(alternative='greater'))$p.value %>% format(digits = 2)

p.combine.branch  <-  combine.pair.all.branch %>% inner_join(aggregate(combine.pair.all.branch$Freq, by=list(combine.pair.all.branch$type), sum) %>% rename(sum = x), by=c("type" = "Group.1") ) %>% mutate(ratio = Freq/sum) %>% ggplot(aes()) + geom_bar(aes(x=branch, y=ratio, fill=type), stat="identity", position = "dodge") + theme_bw() + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), guide = guide_legend( title = "Negative GI type"), labels= c("Core", "Conditional", "Non-nGI") ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text.x=element_text(size=17, angle=30), axis.text.y=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = c(0.8, 0.7) ) + xlab("Gene origination branch") + ylab("Proportion") + geom_text(aes(x=x, y=y, label=label), data=data.frame(x=5.5, y=c(0.87, 0.80), label=paste("P value =", c(k1, k2)), sep=' ' ), size=5.8) + scale_x_discrete(labels = taxon)


branch.name <- data.frame("branch" = 0:13, "Nodes" = taxon)
branch.name$Nodes <- sub("\n", "", branch.name$Nodes )


###When we only focused on WGD pairs.
combine.pair.all.WGD <- combine.pair.all %>% mutate( "Duplication" = ifelse(pair %in% c(ohnologs$pair1, ohnologs$pair2), "WGD", "SSD")) %>% select(type, Duplication) %>% table() %>% as.data.frame() %>% dcast(type~Duplication, value.var="Freq") %>% mutate(ratio = WGD/(WGD + SSD)) %>% mutate( sd = sqrt( ratio*(1-ratio)/(SSD + WGD) ) )

k1 <- (combine.pair.all.WGD %>% filter(type %in% c("Robust", "Conditional") ) %>% select(SSD, WGD) %>% as.matrix() %>% fisher.test(alternative='less'))$p.value %>% format(digits=2)
k2 <- (combine.pair.all.WGD %>% filter(type %in% c("Robust", "Others") ) %>% select(SSD, WGD) %>% as.matrix() %>% fisher.test(alternative='less'))$p.value %>% format(digits=2)

p.combine.ohnologs  <- combine.pair.all.WGD %>% ggplot(aes(x=type,  y=ratio)) + geom_bar(aes(fill = type), stat="identity") + geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), width=0.3) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9")) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Proportion of \nWGD pairs") +  geom_text(aes(x=type, y=y, label=label), data=combine.pair.all.WGD %>% mutate(y=0.98, label = paste(WGD, WGD + SSD, sep='/') ), size=5.8)  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))

p.combine.ohnologs.identity <- combine.pair.all %>% filter( pair %in% c(ohnologs$pair1, ohnologs$pair2 )) %>% inner_join(ensembl, by=c("gene1" = "query", "gene2" = "target") ) %>% ggplot(aes(x=type, y=identity)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(95, 105)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Non-nGI" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Protein identity (WGD pairs)")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI")) + ylim(c(20, 110))



###Supplementary. Human post-natal transcriptome data
gene.expression.postnatal  <- readRDS( "../datafiles/gene.expression.postnatal.rds" )
gene.tau.postnatal  <- readRDS( "../datafiles/gene.tau.postnatal.rds")
gene.tau.postnatal.tissue  <- readRDS( "../datafiles/gene.tau.postnatal.tissue.rds")


p.combine.expression.level.postnatal <- combine.pair.all %>% inner_join( gene.expression.postnatal, by=c("gene1" = "genename") ) %>% inner_join( gene.expression.postnatal, by=c("gene2" = "genename") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2) %>% mutate(TPM = log2(TPM + 1)) %>% ggplot(aes(x=type, y=TPM)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(14, 16)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression level \n(human postnatal)")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))

p.combine.expression.tau.postnatal   <- combine.pair.all %>% inner_join( gene.tau.postnatal, by=c("gene1" = "genename") ) %>% inner_join( gene.tau.postnatal, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% ggplot(aes(x=type, y=tau)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(0.9, 1.0)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression specificity \n(human postnatal)")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI")) + scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), labels = c(0.2, 0.4, 0.6, 0.8, 1.0))

p.combine.expression.tau.postnatal.tissue <- combine.pair.all %>% inner_join( gene.tau.postnatal.tissue, by=c("gene1" = "genename") ) %>% inner_join( gene.tau.postnatal.tissue, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% ggplot(aes(x=type, y=tau)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(1.0, 1.1)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression specificity \n(organ level)")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))



###Chicken duplicate pairs. Expression level. All stages
chicken_gene.expression.postnatal <- readRDS( "../datafiles/chicken.gene.expression.postnatal.rds" )
chicken_gene.tau.postnatal  <- readRDS( "../datafiles/chicken.gene.tau.postnatal.rds")
chicken_gene.tau.postnatal.tissue  <- readRDS( "../datafiles/chicken.gene.tau.postnatal.tissue.rds")


p.combine.expression.level.chicken.postnatal <- combine.pair.all.chicken %>% inner_join( chicken_gene.expression.postnatal, by=c("chicken.x" = "geneID") ) %>% inner_join( chicken_gene.expression.postnatal, by=c("chicken.y" = "geneID") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2) %>% mutate(TPM = log2(TPM + 1)) %>% ggplot(aes(x=type, y=TPM)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(14, 16)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression level \nin chicken genome (post-natal)")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))

p.combine.expression.tau.chicken.postnatal  <- combine.pair.all.chicken %>% inner_join( chicken_gene.tau.postnatal, by=c("chicken.x" = "geneID") ) %>% inner_join( chicken_gene.tau.postnatal, by=c("chicken.y" = "geneID") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% ggplot(aes(x=type, y=tau)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(0.9, 1.0)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression specificity \nin chicken genome")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI")) + scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), labels = c(0.2, 0.4, 0.6, 0.8, 1.0))

p.combine.expression.tau.tissue.chicken.postnatal  <- combine.pair.all.chicken %>% inner_join( chicken_gene.tau.postnatal.tissue, by=c("chicken.x" = "geneID") ) %>% inner_join( chicken_gene.tau.postnatal.tissue, by=c("chicken.y" = "geneID") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% ggplot(aes(x=type, y=tau)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Robust", "Conditional"), c("Robust", "Others")), textsize=5.8, y_position=c(1.0, 1.1)) + scale_fill_manual(values=c("#7372B7", "#F2CA66", "#3BA6D9"), labels=c("Core", "Conditional", "Others" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression specificity \nin chicken genome (organ level)")  + scale_x_discrete(labels=c("Core", "Conditional", "Non-nGI"))



###Controlling Pleiotropy
library('MatchIt')

data <- combine.pair.all %>% inner_join( gene.expression.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.expression.prenatal, by=c("gene2" = "genename") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2) %>% mutate(TPM = log2(TPM + 1)) %>% inner_join( gene.tau.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.tau.prenatal, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% inner_join(Go_sum %>% mutate(type = NULL), by=c("gene1" = "gene1", "gene2" = "gene2") ) %>% mutate( BP_sum = BP_sum, MF_sum = MF_sum + 1, CC_sum = CC_sum + 1 )

data_raw         <- data %>% mutate(type=ifelse(type=='Robust', 'Y', 'N'))
data_raw$type    <- factor(data_raw$type,levels=c('N','Y'))

data_match_out    <- matchit(type~tau+BP_sum, data=data_raw, method='nearest', distance='glm')
data_matched      <- match.data( data_match_out )
data_matched$type <- factor(data_matched$type,levels=c('Y','N'))

p_matchsample     <- data_matched  %>%  ggplot(aes(x=type, y=TPM)) + geom_violin(aes(fill = type)) + geom_boxplot(outlier.shape=NA, width=0.2) + geom_signif(comparisons = list(c("Y", "N")),  textsize=5.8, y_position=15) + scale_fill_manual(values=c("#7372B7", "#AAAAAA"), labels=c("Core", "Non-core" ) ) + theme_bw() + theme(title=element_text(size=20, face="bold"), axis.text=element_text(size=17), legend.text=element_text(size=17), plot.title=element_text(size=20, face="bold"), plot.tag=element_text(size=28, face="bold"), legend.position = "none" ) + xlab("Negative GI type") + ylab("Expression level in  prenatal samples\n(pleiotropy controlled)") + scale_x_discrete(labels=c("Core", "Non-core"))



###Compositary plot
pdf("../Figure5_Robust.gi.Features.identity_2021Oct.pdf", width=18, height=18)
grid.arrange(  p.combine.branch + labs(tag="A"),  p.GO.BP + labs(tag="B"), p.combine.identity + labs(tag="C"), p.combine.expression.tau  + labs(tag="C"), p.combine.expression.level + labs(tag="C"), p.combine.identity.chicken + labs(tag="D"), p.combine.expression.tau.chicken  + labs(tag="D"), p.combine.expression.level.chicken  + labs(tag="D"), layout_matrix = matrix(c( rep(1,2), 2:8 ), nrow=3, byrow=T) )
dev.off()

pdf("../FigureS5_Robust.gi.Features.identity.supplementary6_2021Oct.pdf", width=20, height=18)
grid.arrange(  p.combine.ohnologs + labs(tag="A"),  p.GO.CC + labs(tag="C"), p_matchsample + labs(tag='D'), 
p.combine.expression.tau.tissue + labs(tag="B"), p.combine.expression.tau.postnatal + labs(tag="B"),  p.combine.expression.tau.postnatal.tissue + labs(tag="B"),  p.combine.expression.level.postnatal + labs(tag="B"),  p.combine.expression.tau.tissue.chicken + labs(tag="B"), p.combine.expression.tau.chicken.postnatal + labs(tag="B"), p.combine.expression.tau.tissue.chicken.postnatal + labs(tag="B"),   p.combine.expression.level.chicken.postnatal + labs(tag="B"), layout_matrix = matrix( rep(1:11,c(4,4,4,3,3,3,3,3,3,3,3)), nrow=3, byrow=T) )
dev.off()



###Saving all data.
write.csv( 
all.pair.gi.count %>% inner_join(combine.pair.all, by=c("pair")) %>% left_join(ensembl, by=c("gene1" = "query", "gene2" = "target") ) %>% select(gene1, gene2, pair, type, identity, corum) %>% left_join( genename %>% inner_join( gentree, by=c("V1" = "gene")) %>% select(V2, branch), by=c("gene1" = "V2")) %>% left_join( genename %>% inner_join( gentree, by=c("V1" = "gene")) %>% select(V2, branch), by=c("gene2" = "V2")) %>% mutate(branch = ifelse(branch.x < branch.y, branch.y, branch.x) ) %>% left_join(branch.name, by="branch") %>% mutate( "Duplication" = ifelse(pair %in% c(ohnologs$pair1, ohnologs$pair2), "WGD", "SSD")) %>% left_join(Go_sum %>% mutate(type = NULL), by=c("gene1" = "gene1", "gene2" = "gene2") ) %>% 

left_join( combine.pair.all %>% left_join( gene.expression.prenatal, by=c("gene1" = "genename") ) %>% left_join( gene.expression.prenatal, by=c("gene2" = "genename") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2), by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 
left_join( combine.pair.all %>% left_join( gene.tau.prenatal, by=c("gene1" = "genename") ) %>% left_join( gene.tau.prenatal, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2),  by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 
left_join( combine.pair.all %>% left_join( gene.tau.prenatal.tissue, by=c("gene1" = "genename") ) %>% left_join( gene.tau.prenatal.tissue, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2), by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 
left_join( combine.pair.all %>% left_join( gene.expression.postnatal, by=c("gene1" = "genename") ) %>% left_join( gene.expression.postnatal, by=c("gene2" = "genename") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2), by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 
left_join( combine.pair.all %>% left_join( gene.tau.postnatal, by=c("gene1" = "genename") ) %>% left_join( gene.tau.postnatal, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2),  by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 
left_join( combine.pair.all %>% left_join( gene.tau.postnatal.tissue, by=c("gene1" = "genename") ) %>% left_join( gene.tau.postnatal.tissue, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2), by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 

left_join( combine.pair.all.chicken %>% left_join( chicken_gene.expression.prenatal, by=c("chicken.x" = "geneID") ) %>% left_join( chicken_gene.expression.prenatal, by=c("chicken.y" = "geneID") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2), by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 
left_join( combine.pair.all.chicken %>% left_join( chicken_gene.tau.prenatal, by=c("chicken.x" = "geneID") ) %>% left_join( chicken_gene.tau.prenatal, by=c("chicken.y" = "geneID") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) , by=c("gene1", "gene2", "pair", "type", "corum") ) %>%
left_join( combine.pair.all.chicken %>% left_join( chicken_gene.tau.prenatal.tissue, by=c("chicken.x" = "geneID") ) %>% left_join( chicken_gene.tau.prenatal.tissue, by=c("chicken.y" = "geneID") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) , by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 

left_join( combine.pair.all.chicken %>% left_join( chicken_gene.expression.postnatal, by=c("chicken.x" = "geneID") ) %>% left_join( chicken_gene.expression.postnatal, by=c("chicken.y" = "geneID") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2), by=c("gene1", "gene2", "pair", "type", "corum") ) %>% 
left_join( combine.pair.all.chicken %>% left_join( chicken_gene.tau.postnatal, by=c("chicken.x" = "geneID") ) %>% left_join( chicken_gene.tau.postnatal, by=c("chicken.y" = "geneID") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) , by=c("gene1", "gene2", "pair", "type", "corum") ) %>%
left_join( combine.pair.all.chicken %>% left_join( chicken_gene.tau.postnatal.tissue, by=c("chicken.x" = "geneID") ) %>% left_join( chicken_gene.tau.postnatal.tissue, by=c("chicken.y" = "geneID") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) , by=c("gene1", "gene2", "pair", "type", "corum") )

, "../TableS10_combine.pair.gi.features.summary.csv" )



###Testing the independency of pleiotropy and expression level with multiple linear regression analysis.
combine.pair.all$type <- factor( combine.pair.all$type, levels=c("Others", "Conditional", "Robust") )

( lm(TPM~type + tau + BP_sum, data =combine.pair.all %>% inner_join( gene.expression.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.expression.prenatal, by=c("gene2" = "genename") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2) %>% mutate(TPM = log2(TPM + 1)) %>% inner_join( gene.tau.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.tau.prenatal, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% inner_join(Go_sum %>% mutate(type = NULL), by=c("gene1" = "gene1", "gene2" = "gene2") ) %>% mutate( BP_sum = log10(BP_sum + 1), MF_sum = log10(MF_sum + 1), CC_sum = log10(CC_sum + 1) ) ) %>% summary() )$coefficients

( lm(tau~type + TPM + BP_sum, data =combine.pair.all %>% inner_join( gene.expression.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.expression.prenatal, by=c("gene2" = "genename") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2) %>% mutate(TPM = log2(TPM + 1)) %>% inner_join( gene.tau.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.tau.prenatal, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% inner_join(Go_sum %>% mutate(type = NULL), by=c("gene1" = "gene1", "gene2" = "gene2") ) %>% mutate( BP_sum = log10(BP_sum + 1), MF_sum = log10(MF_sum + 1), CC_sum = log10(CC_sum + 1) ) ) %>% summary() )$coefficients

( lm(BP_sum~type + TPM + tau, data =combine.pair.all %>% inner_join( gene.expression.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.expression.prenatal, by=c("gene2" = "genename") ) %>% mutate(TPM = TPM.x + TPM.y) %>% mutate(TPM = TPM/2) %>% mutate(TPM = log2(TPM + 1)) %>% inner_join( gene.tau.prenatal, by=c("gene1" = "genename") ) %>% inner_join( gene.tau.prenatal, by=c("gene2" = "genename") ) %>% mutate(tau = tau.x + tau.y) %>% mutate(tau = tau/2) %>% inner_join(Go_sum %>% mutate(type = NULL), by=c("gene1" = "gene1", "gene2" = "gene2") ) %>% mutate( BP_sum = log10(BP_sum + 1), MF_sum = log10(MF_sum + 1), CC_sum = log10(CC_sum + 1) ) ) %>% summary() )$coefficients



