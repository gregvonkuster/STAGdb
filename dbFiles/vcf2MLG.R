setwd("F:/PSU/Hybridization Project/SNPs/database")

#Required R-packages for multi-locus genotype calling
library(vcfR)
library(poppr)
library(adegenet)
library(vegan)
library(ape)

#Read in VCF file with array SNVs
vcf <- read.vcfR('baitsSNV.recode.vcf')

#Missing GT in samples submitted
library(RColorBrewer)
palette(brewer.pal(n=12, name = 'Set3'))

gt <- extract.gt(vcf, element = "GT", as.numeric=FALSE)
myMiss <- apply(gt, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- (myMiss/nrow(vcf))*100

write.table(myMiss,file="missingData.txt")

#Convert VCF file into a genind needed for Poppr package
gind <- vcfR2genind(vcf)

#add population or species information to the genind pop slot
poptab<-read.table("F:/PSU/Hybridization Project/SNPs/database/112017_databaseDevelopment/popInfo.txt",
                   check.names=FALSE, header=T, na.strings = c("", "NA"))

gind@pop <- as.factor(poptab$region)

obj2<-as.genclone(gind)
obj2

#calculate the bitwise distance between individuals
xdis<-bitwise.dist(obj2)

#multilocus genotypes (threshold of 1%)
mlg.filter(obj2, distance= xdis) <- 0.01 #threshold of 0.01
m<-mlg.table(obj2, background=TRUE, color=TRUE)
nmll(obj2)

#create table of MLGs
id<-mlg.id(obj2)
df <- data.frame(matrix((id), byrow=T))
write.table(df,file="mlg-id.txt")

#rarifaction curve
H.obj <- mlg.table(obj2, plot = TRUE)
pdf ("geno_rareCurve.pdf", width=10, height=7)
rarecurve(H.obj, ylab="Number of expected MLGs", sample=min(rowSums(H.obj)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")
dev.off()

#genotype accumulation curve, sample = number of loci randomly selected for to make the curve
pdf ("geno_accumulationCurve.pdf", width=10, height=7)
gac <- genotype_curve(gind, sample = 5, quiet = TRUE)
gac
dev.off()

##Create a phylogeny of samples based on distance matrices
cols <- c("skyblue2","#C38D9E", '#E8A87C',"darkcyan","#e27d60")
set.seed(999)
pdf ("NJphylogeny.pdf", width=10, height=7)
theTree <- obj2 %>%
  aboot(dist = provesti.dist, sample = 100, tree = "nj",
        cutoff = 50, quiet = TRUE) %>%
  ladderize() # organize branches by clade
plot.phylo(theTree, tip.color = cols[obj2$pop],label.offset = 0.0125,cex=0.7, 
           font=2, lwd=4)
add.scale.bar(length = 0.05, cex=0.65) # add a scale bar showing 5% difference.
nodelabels(theTree$node.label, cex=.5, adj = c(1.5, -0.1), frame = "n", font=3,
           xpd = TRUE)
dev.off()

