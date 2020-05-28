gene.exp <- as.matrix(read.csv('data/mRNA_exp_FPKM.txt',
                               sep = '\t', header = T, row.names =1))
gene.exp <- log(gene.exp + 0.0001)

# 提取正常样本和疾病样本的表达
sample <- colnames(gene.exp)
normal.exp <- gene.exp[, grep('.*11$', sample)]
tumor.exp <- gene.exp[, grep('.*01$', sample)]

# 离散化
cutoff.left <- apply(normal.exp, 1, function(x){
  mean(x) - sd(x)
})
cutoff.right <- apply(normal.exp, 1, function(x){
  mean(x) + sd(x)
})
cutoff <- data.frame(cutoff.left, cutoff.right)

disease.exp <- matrix(0, nrow = nrow(tumor.exp), ncol = ncol(tumor.exp))
colnames(disease.exp) <- colnames(tumor.exp)
rownames(disease.exp) <- rownames(tumor.exp)

for (i in 1:nrow(tumor.exp)) {
  disease.exp[i, which(tumor.exp[i,] < cutoff$cutoff.left[i])] <- -1
  disease.exp[i, which(tumor.exp[i,] > cutoff$cutoff.right[i])] <- -1
}

# 导入CN数据
disease.cnv <- as.matrix(read.csv('data/CNV_results/all_thresholded.by_genes.txt',
                                  sep = '\t', header = T, row.names = 1))
disease.cnv <- disease.cnv[, -(1:2)]
colnames(disease.cnv) <- substr(colnames(disease.cnv), 1, 15)

comm.row <- intersect(rownames(disease.exp), rownames(disease.cnv))
comm.col <- intersect(colnames(disease.exp), colnames(disease.cnv))
data.exp <- disease.exp[comm.row, comm.col]
data.cnv <- disease.cnv[comm.row, comm.col]
data.cnv <- apply(data.cnv, 2, as.numeric)
rownames(data.cnv) <- comm.row

for (i in 1:nrow(data.cnv)) {
  data.cnv[i, which(data.cnv[i,] == 2)] <- 1
  data.cnv[i, which(data.cnv[i,] == -2)] <- -1
}

# 计算DES得分
gene.DES <- matrix(0, nrow(data.exp), 2)
gene.DES[,1] <- rownames(data.exp)
colnames(gene.DES) <- c('Gene_Symbol', 'DES')

for (i in 1:nrow(data.exp)) {
  amp.index <- which(data.cnv[i,] == 1)
  del.index <- which(data.cnv[i,] == -1)
  amp.num <- length(which(data.exp[i, amp.index] == 1))
  del.num <- length(which(data.exp[i, del.index] == -1))
  DES.num <- amp.num + del.num
  cnv.num <- length(amp.index) + length(del.index)
  gene.DES[i,2] <- DES.num/cnv.num
}

# 输出
write.table(gene.DES, 'DES.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)
