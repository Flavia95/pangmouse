library (VennDiagram)
library(grid)
library(scales)
library(tidyverse)

#vennposVCFcompare,another way to obtained a venn for variable sites from VCFcomparepos.txt
myd = read.csv("vcfcomparepos.txt",sep = "",header = F)
my_data = myd %>% select(V1)
final_df <- as.data.frame(t(my_data))
colnames(final_df) = c ("Genomic","Inters","Pang")
ListGeno <- select(final_df, Genomic) 
ListPang <- select(final_df, Pang) 
ListInt <- select(final_df, Inters)
png('vennpositions.png')
venn.plot <- draw.pairwise.venn(ListPang+ListInt, ListGeno+ListInt, ListInt, category=c('GATK', 'VG'),fill = c("light blue", "#bbffbb"))
dev.off()


#Piechart_vtpeek

myd = read.delim("types_variants_geno.txt",  sep="\t", header = T)
df = data.frame(myd)
c = ggplot(df, aes (x="", y = Count, fill = factor(Type))) + 
geom_col(position = 'stack', width = 1) +geom_text(aes(label = paste(round(Count / sum(Count) * 100, 1), "%"), x = 1.3), position = position_stack(vjust = 0.6), size = 2.5 , fontface="bold") + theme_classic() + theme(plot.title = element_text(hjust=0.7), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + labs(fill = "Type", x = NULL, y = NULL, title = "Variants distribution on pangenome of 39 strains after smoothxg") + coord_polar("y")
ggsave("./distributionofvariants_geno.png", plot= c, device="png", width = 20, height = 15, units = "cm", dpi = 300)


myd = read.delim("types_variants_pang.txt",  sep="\t", header = T)
df = data.frame(myd)
c = ggplot(df, aes (x="", y = Count, fill = factor(Type))) +
geom_col(position = 'stack', width = 1) +geom_text(aes(label = paste(round(Count / sum(Count) * 100, 1), "%"), x = 1.3), position = position_stack(vjust = 0.6), size = 2.5 , fontface="bold") + theme_classic() + theme(plot.title = element_text(hjust=0.7), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + labs(fill = "Type", x = NULL, y = NULL, title = "Variants distribution on pangenome of 39 strains after smoothxg") + coord_polar("y")
ggsave("./distributionofvariants_pang.png", plot= c, device="png", width = 20, height = 15, units = "cm", dpi = 300)



#vennvariantsbcftools

myd = read.csv("overlapvariants.tsv", header = T, sep = "\t")
colnames(myd) <- c("SN", "ID","Type", "Count")                                                     
variants = c("number of records:")
Rec <- filter (myd,Type %in% variants)
my_data = Rec %>% select(Count)
final_df <- as.data.frame(t(my_data))
colnames(final_df) = c("Geno","Pang","Int")
Geno <- select(final_df, Geno) 
Pang <- select(final_df, Pang) 
Int <- select(final_df, Int)
png('vennvariants.png')
venn.plot <- draw.pairwise.venn(Pang+Int, Geno+Int, Int, category=c('VG', 'GATK'),fill = c("light blue", "#bbffbb"))
dev.off()



#Concordance and Discordance Classes
myd = read.csv("overlapvariants.tsv", header = T, sep = "\t")
colnames(myd) <- c("SN", "ID","Type", "Count")                                                     
variants = c("number of SNPs:", "number of MNPs:", "number of indels:", "number of others:")
OnlyGeno <- filter (myd, ID == 0 & Type %in% variants)
mydG<- spread(OnlyGeno, Type, Count)  
colnames(mydG)<- c("SN", "ID","Indels", "MNPs","Others","SNPs") 
mydG$ID  <- "Only Genomic"                                                                           

OnlyPang <- filter (myd, ID == 1 & Type %in% variants)
mydP<- spread(OnlyPang, Type, Count)
colnames(mydP)<- c("SN", "ID","Indels", "MNPs","Others","SNPs")
mydP$ID<-" Only Pangenomic"                                                                                                                                                       

Both <- filter (myd, ID == 2 & Type %in% variants)
mydB<- spread(Both, Type, Count)
colnames(mydB)<- c("SN", "ID","Indels", "MNPs","Others","SNPs")
mydB$ID  <- "Both"                                                                            


tab <- rbind (mydB,mydG,mydP)
NewTab = tab %>% gather(Type,Variants, starts_with('SNPs'), starts_with('MNPs'),starts_with('Indels'), starts_with('Others'))
NewTab$Variants = as.numeric(NewTab$Variants)                  
NewTab$Type <- factor(NewTab$Type, levels = c("Mnps", "SNPs", "Indels","Others"))
myPlot = ggplot(NewTab, aes(ID, Variants/sum(Variants)*100, fill= Type)) + geom_bar (stat = "identity", position = "dodge") + theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5)) + ylab("Percentage") + xlab("Class") + ggtitle("Discordance and Concordance Classes") + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_discrete(name = "Variant Type") 
g= myPlot + scale_fill_discrete(labels=c("SNPs", "INDELs", "Others", "MNPs")) + theme_bw()
ggsave("./ConcordanceandDiscordanceclasses.png", plot= g, device="png", width = 20, height = 15, units = "cm", dpi = 300)

#Distribution Indels

myd = read.csv("geno.indel.hist", sep="\t", header=T)
h <- (subset(myd, LENGTH!=0))
h$Methods  <- "Genomic"
myd2 = read.csv("pang.indel.hist", sep="\t", header=T)
h1 <- (subset(myd2, LENGTH!=0))
h1$Methods  <- "Pangenomic"
tab <- rbind (h,h1)
p = ggplot(tab)+
  geom_col(aes(x = LENGTH, y = COUNT, fill=Methods)) +xlim(-200,200) +ylim(0,10000) + facet_grid(. ~ Methods) +theme_bw() + theme(legend.position = "none")
 ggsave("./Lenindelspangvsgen.png", plot= p, device="png", width = 20, height = 15, units = "cm", dpi = 300) 
  
#Zoom indels
lenindels<-ggplot(tab, aes(LENGTH, COUNT , fill=Methods) ) + geom_col  (position='dodge') + xlim(-1000, 1000) + scale_y_log10() + theme_bw()
ggsave("./Lenindelszoompangvsgen.png", plot= lenindels, device="png", width = 20, height = 15, units = "cm", dpi = 300)

#Bonus plot

#1. Quality
myd = read.delim("qual_sites_gen.lqual", sep = "\t", skip = 1)
colnames(myd) = c("chr", "pos", "qual")
myd1 = read.delim("qual_sites_pan.lqual", sep = "\t", skip = 1)
colnames(myd1) = c("chr", "pos", "qual")
myd$Methods  <- "Genomic"
myd1$Methods  <- "Pangenomic"
h = rbind(myd,myd1)
p <- ggplot(h, aes(qual, fill = Methods), position = "identity") + geom_density(alpha = 0.3) + xlim(20,1000) + theme_light()
ggsave("./qualitypangvsgen.png", plot= lenindels, device="png", width = 20, height = 15, units = "cm", dpi = 300)

#2. MAF
myd = read.delim("freq_geno.frq", sep = "\t", skip = 1)
colnames(myd) = c("chr", "pos", "nalleleles","nchr","a1","a2")
myd$maf <- myd %>% select(a1,a2) %>% apply(1,function(z)min(z))
myd1 = read.delim("freq_pan.frq", sep = "\t", skip = 1)
colnames(myd1) = c("chr", "pos", "nalleleles","nchr","a1","a2")
myd1$maf <- myd1 %>% select(a1,a2) %>% apply(1,function(z)min(z))
myd$Methods  <- "Genomic"
myd1$Methods  <- "Pangenomic"
h = rbind(myd,myd1)
a <- ggplot(h, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0.002, 0.8)
x= a+ facet_grid(. ~ Methods) + theme_light()
ggsave("./MAFpangvsgen.png", plot= lenindels, device="png", width = 20, height = 15, units = "cm", dpi = 300)

#3. Ibread
ind_het1 <- read_delim("het_gen.het", delim = "\t",col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
a <- ggplot(ind_het1, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()
ggsave("./inbread.png", plot= lenindels, device="png", width = 20, height = 15, units = "cm", dpi = 300)
