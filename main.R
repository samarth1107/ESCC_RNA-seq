#####################
# Samarth Chauhan
# 2018410
# CSB
#####################



#####################
getwd()
# change the working directory
setwd("E:/Project/RNA-seq")
getwd()
#####################




#####################
# TO do filtering and DE analysis on HTseq data
library(DESeq2)
#####################




#####################
#To take files from HTseq as count matix and covert to DESeq data type
#Take all files name
CM_Files <- grep("Tumour", 
                 list.files("Count_matrix/"),
                 value=TRUE)

CM_Samples <- c("Non_Tumour_1","Non_Tumour_2","Non_Tumour_3","Tumour_1","Tumour_2","Tumour_3")

#Take all samples here Tumour and Non-Tumour
CM_Condition <- sub("(.*Tumour).*",
                         "\\1",
                         CM_Files)

#Make table which is in the form of DE Seq data 
Table <- data.frame(sampleName = CM_Files,
                    fileName = CM_Files,
                    condition = CM_Condition)
Table$condition <- factor(Table$condition)

#Make DESeqDataSet data
DESeq_Data <- DESeqDataSetFromHTSeqCount(sampleTable = Table,
                                       directory = "Count_matrix/",
                                       design= ~ condition)
#####################





#####################
#Differenctial gene expressoion
DESeq <- DESeq(DESeq_Data)
DESeq$type<-c('single-read','single-read','single-read','single-read','single-read','single-read')
DESeq$type<-factor(c('single-read'))


DESeq_result <- results(DESeq,
                        pAdjustMethod = "BH",
                        alpha = 0.1)

####################
#convert gene names
library("org.Hs.eg.db")
DESeq_result$hgnc_symbol <- mapIds(org.Hs.eg.db,
                                   keys=gsub("\\..*","",row.names(DESeq_result)),
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")
DESeq_result$entrezid <- mapIds(org.Hs.eg.db,
                                keys=gsub("\\..*","",row.names(DESeq_result)),
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first")
library('edgeR')
DESeq_result$logCPM <- rowMeans(log(cpm(assay(DESeq_Data))))
DESeq_result$SD <- rowSds(log(cpm(assay(DESeq_Data))))


sum(DESeq_result$pvalue < 0.05, na.rm=TRUE)
sum(DESeq_result$padj < 0.1, na.rm=TRUE)
summary(DESeq_result)

#MA plot for dispersion of data
plotMA(DESeq_result, ylim=c(-10,10))

#Dispersion  estimation
plotDispEsts(DESeq)

plot(DESeq_result$logCPM, DESeq_result$SD, pch=".", main="Distribution of data", ylab="sd", xlab="Average logCPM")

#Standard deviationn
norm_DESeq <- normTransform(DESeq)
variance_stable_DESeq <- varianceStabilizingTransformation(DESeq,
                                                           blind = TRUE,
                                                           fitType = 'parametric')
log_DESeq <- rlog(DESeq)


#PCA
library("gplots")
library("RColorBrewer")
plotPCA(log_DESeq, 
        intgroup=c("condition", "type"))

sampleDists <- dist(t(assay(log_DESeq)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(log_DESeq$dex, log_DESeq$cell, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )



#Distance between sample to sample
library("RColorBrewer")
library("pheatmap")
sample_distance <- dist(t(assay(log_DESeq)))
sample_distance_matrix <- as.matrix(sample_distance)
rownames(sample_distance_matrix) <- paste(log_DESeq$condition, 
                                          log_DESeq$type, 
                                          sep="-")
colnames(sample_distance_matrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sample_distance_matrix,
         clustering_distance_rows=sample_distance,
         clustering_distance_cols=sample_distance,
         col=colors,
         legend_labels = CM_Condition,
         labels_row = CM_Samples,
         labels_col = CM_Samples)





#histogram with p value
hist(DESeq_result$pvalue, 
     breaks = 0:50/50, 
     main = "p-value histogram",
     xlab="p-value",
     ylab="frequency",
     plot = TRUE)

hist(DESeq_result$pvalue[DESeq_result$baseMean > 1], 
     breaks = 0:50/50, 
     main = "p-value histogram",
     xlab="p-value",
     ylab="frequency",
     plot = TRUE)

use <- DESeq_result$pvalue < 0.05
h1 <- hist(DESeq_result$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(DESeq_result$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="red", `pass`="green")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


#MA plot
plotMA(DESeq_result,
       xlab = "mean of count (normalized)",
       main = "Differential vs Expression")




#######
library("vsn")
meanSdPlot(assay(norm_DESeq))
meanSdPlot(assay(variance_stable_DESeq))




# create bins using the quantile function
qs <- c( 0, quantile( DESeq_result$baseMean[DESeq_result$baseMean > 0], 0:7/7 ) )
# "cut" the genes into the bins
bins <- cut( DESeq_result$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of p values less than .01 for each bin
ratios <- tapply( DESeq_result$pvalue, bins, function(p) mean( p < .05, na.rm=TRUE ) )
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")


metadata(DESeq_result)$filterThreshold


plot(metadata(DESeq_result)$filterNumRej,type="b",
     xlab="quantiles of 'baseMean'",
     ylab="number of rejections")



#####################




#####################
#gene selection 
filtered_genes <- as.data.frame(DESeq_result[order(DESeq_result$padj),])

filtered_genes <- filtered_genes[!(filtered_genes$baseMean==0),]
filtered_genes <- filtered_genes[!is.na(filtered_genes$pvalue),]
filtered_genes <- filtered_genes[!is.na(filtered_genes$padj),]
filtered_genes <- filtered_genes[!is.na(filtered_genes$hgnc_symbol),]
filtered_genes <- filtered_genes[!(filtered_genes$logCPM==-Inf),]

filtered_genes <- filtered_genes[(filtered_genes$pvalue<0.05),]
filtered_genes <- filtered_genes[(filtered_genes$padj<0.1),]
filtered_genes <- filtered_genes[(filtered_genes$logCPM>2),]
filtered_genes <- filtered_genes[(filtered_genes$log2FoldChange>4 | filtered_genes$log2FoldChange< -4),]

write.table(c(filtered_genes$hgnc_symbol), 
            file = "names.txt", sep = "\n",
            row.names = FALSE)
#####################



#####################
#gene ontology
library(GOfuncR)
library(Homo.sapiens)


gene_ids = c(filtered_genes$hgnc_symbol)
input_hyper = data.frame(gene_ids, is_candidate=1)
res_hyper = go_enrich(input_hyper, n_randset=100)

res_hyper = go_enrich(input_hyper, godir="E:/Project/RNA-seq/Ontology/go_weekly-termdb-tables", )

top_gos_hyper = res_hyper[[1]][1:5, 'node_id']
plot_anno_scores(res_hyper, top_gos_hyper)
