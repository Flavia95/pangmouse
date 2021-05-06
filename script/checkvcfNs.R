#After run the script
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' genomic.vcf.gz > checkposNgeno.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' pangenomic.vcf.gz > checkposNpang.vcf

#For check numbers of positions that have N as ALT and REF allele

myd = read.csv("pangenomic.vcf", sep="\t", header=F)
colnames(myd) = c("CHROM","POS","REF","ALT")
p<- myd %>% select (POS,REF, ALT)
x = filter(p, grepl('N', ALT))
d = data.frame(cbind(dim(p),dim(x)))
dt <- d %>% slice(1)
colnames(dt)<- c("length", "Ns")
dt$Methods = "Pangenomic"

myd1 = read.csv("genomic.vcf", sep="\t", header=F)
colnames(myd1) = c("CHROM","POS","REF","ALT")
p1<- myd1 %>% select (POS,REF, ALT)
x1 = filter(p1, grepl('N', ALT))
d1 = data.frame(cbind(dim(p1),dim(x1)))
dt1 <- d1 %>% slice(1)
colnames(dt1)<- c("length", "Ns")
dt1$Methods = "Genomic"
Nasalt = rbind(dt,dt1)

write.table (Nasalt, quote = F, sep = " \t",row.names = F, "NumbersofNasALTallele.txt") #same for REF allele
