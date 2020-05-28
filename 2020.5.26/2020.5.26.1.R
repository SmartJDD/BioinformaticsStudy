# 处理拷贝数片段信息
cnv <- as.matrix(read.csv('../第二次课/data/CNV_results/all_lesions.conf_90.txt',
                          sep = '\t', header = T))
cnv <- cnv[1:(dim(cnv)[1]/2),]
infor <- matrix(0, dim(cnv)[1], 3)
for (i in 1:dim(cnv)[1]) {
  split.info <- unlist(strsplit(cnv[i,5], '(', fixed = T))
  split.info1 <- unlist(strsplit(split.info[1], ':', fixed = T))
  split.info2 <- unlist(strsplit(split.info1[2], '-', fixed = T))
  infor[i,1] <- split.info1[1]
  infor[i,2] <- split.info2[1]
  infor[i,3] <- split.info2[2]
}
colnames(infor) <- c('chomr', 'start', 'end')
infor <- cbind(infor, cnv[,1:2])

write.table(infor, 'cnv_infor.txt',
            sep = "\t", col.names = T,row.names = F,quote = F)

# 基因比对
gene.info <- as.matrix(read.csv("../第二次课/data/protein_coding_gene_hg38.txt",
                                sep = "\t",header = F))

result<-c()
for(i in 1:dim(gene.info)[1]) {
  for(j in 1:dim(infor)[1]) {
    if (gene.info[i,3] == infor[j,1] & 
        (as.numeric(infor[j,2]) <= as.numeric(gene.info[i,5]) &
         as.numeric(infor[j,3]) >= as.numeric(gene.info[i,5]) |
         as.numeric(infor[j,3]) >= as.numeric(gene.info[i,6]) &
         as.numeric(infor[j,2]) <= as.numeric(gene.info[i,6])))
    {
      result_mid <- c(gene.info[i,2], infor[j,])
      result_mid <- t(result_mid)
      result <- rbind(result,result_mid)
    }
  }
}

result <-result[!duplicated(result[,1]),]
write.table(result,"cnv_gene.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)
Amp.gene <- result[grep("Amplification",result[,5]),]
Del.gene <- result[grep("Deletion",result[,5]),]

# 差异基因识别--limma
library(limma)
gene.exp <- as.matrix(read.csv("../第二次课/data/gene_exp_matrix.txt",
                               sep = "\t",header = T))
# 构建分组矩阵
sample <- colnames(gene.exp)
normal <- grep('.*11$', sample)
cancer <- grep('.*01$', sample)
tmp1 <- matrix(0, ncol = 1, nrow = length(sample))
tmp2 <- matrix(0, ncol = 1, nrow = length(sample))
tmp1[normal] <- 1
tmp2[cancer] <- 1
design.mat <- cbind(tmp1, tmp2)
colnames(design.mat) <- c('normal', 'cancer')
rownames(design.mat) <- sample

# 构建差异比较矩阵
contrast.mat<-makeContrasts(paste0(c('normal','cancer'),collapse = "-"), levels = design.mat)

fit <- lmFit(gene.exp, design.mat)
fit2 <- contrasts.fit(fit, contrast.mat)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
DEGene <- rownames(nrDEG)[which(nrDEG$adj.P.Val < 0.05)]
UP.gene <- rownames(nrDEG)[which((nrDEG$adj.P.Val < 0.05) & (nrDEG$logFC > 0))]
DOWN.gene <- rownames(nrDEG)[which((nrDEG$adj.P.Val < 0.05) & (nrDEG$logFC < 0))]

# 识别driver gene
Driver.gene <- c(intersect(UP.gene, Amp.gene), intersect(DOWN.gene, Del.gene))

# 功能富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
genelist <- bitr(DEGene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
genelist <- genelist$ENTREZID
go <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont = 'ALL',
               pAdjustMethod = 'BH', pvalueCutoff = 0.05,
               qvalueCutoff = 0.1, keyType = 'ENTREZID')
kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg',
                   pAdjustMethod = 'BH', pvalueCutoff = 0.05,
                   minGSSize = 10, maxGSSize = 500,
                   qvalueCutoff = 0.2,
                   use_internal_data = FALSE)
result_go<-go@result
result_kegg<-kegg@result
barplot(go, showCategory = 20, drop = T)
dotplot(go, showCategory = 20)
