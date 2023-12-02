#Reading and preparation of count files as the output of htseq-count
getwd()
setwd("C:/Users/User/Documents/HCT116 vs SW48/")
files<- list.files("C:/Users/User/Documents/HCT116 vs SW48/", ".*txt")
htcounts<-lapply(files, read.delim , header=F, comment.char="-")
htcounts<-do.call(cbind, htcounts)
dim(htcounts)
head(htcounts)
row.names(htcounts)<-htcounts [,1]
head(htcounts)
htcounts<-htcounts[, -seq(1, ncol(htcounts), 2)]
head(htcounts)
colnames(htcounts)=sub(".txt", "", files)
head(htcounts)
dim(htcounts)

#Normalization and differential expression analysis with DESeq2
library("DESeq2")
gr=factor(c(rep("HCT",4),rep("SW48",3)))
colData<-data.frame(group=gr, type="paired-end")   
dds=DESeqDataSetFromMatrix(countData = htcounts, colData = colData, design = ~ group)
normalizeddds=DESeq(dds)


#final results of Normalized counts
res=results(normalizeddds)
normalizedcounts=counts(normalizeddds, normalized=T)

#Merge final results with raw counts
resdata= merge(as.data.frame(res), as.data.frame(rawcounts), by="row.names",sort=FALSE)
head(resdata)
write.table(resdata, file = "Final Results.txt", dec = ".", sep = ";", quote = FALSE)


#MAPLOT 
pdf("MAplot.pdf")
plotMA(normalizeddds, ylim=c(-15,15),xlab="Mean of normalized counts")
dev.off()

#Heatmap
install.packages("gplots")
library(gplots)
vsd=varianceStabilizingTransformation(normalizeddds, blind = T)=res[order(res$padj),]
res=res[order(res$padj),]
head(res)
select= head(row.names(res), 1000)
my.pallete=colorRampPalette(c("blue", "white", "red"))(n=1000)
class(vsd)
vsd=assay(vsd)
class(vsd)
pdf("Heatmap.pdf")
heatmap.2(vsd[select,],col = my.pallete,scale = "row", key = "T", keysize = 1, symkey = "T", density.info = "none", trace = "none", cexCol = 1, labRow = T, main = "" )
dev.off()


