gene.exp <- as.matrix(read.csv('data/mRNA_exp_FPKM.txt',
                               sep = '\t', header = T, row.names =1))
gene.exp <- log(gene.exp + 0.0001)
gene.cnv <- as.matrix(read.csv('data/CNV_results/all_data_by_genes.txt',
                               sep = '\t', header = T, row.names =1))
gene.cnv <- gene.cnv[, -(1:2)]
colnames(gene.cnv) <- substr(colnames(gene.cnv), 1, 15)

# 获取重叠的数据
comm.row <- intersect(rownames(gene.exp), rownames(gene.cnv))
comm.col <- intersect(colnames(gene.exp), colnames(gene.cnv))
sen.cnv <- gene.cnv[comm.row, comm.col]
sen.exp <- gene.exp[comm.row, comm.col]
sen.cnv <- apply(sen.cnv, 2, as.numeric)
rownames(sen.cnv) <- comm.row
rownames(sen.exp) <- comm.row

# 筛选基因计算MS以及DSS
r <- nrow(sen.exp)
gene.DSS <- matrix(0, r, 2)
colnames(gene.DSS) <- c('Gene_Symbol', 'DSS')
gene.DSS[,1] <- rownames(sen.exp)
n <- 6

for (i in 1:r) {
  fit <- loess(sen.exp[i,] ~ sen.cnv[i,], span = 2/3, degree = 1, family = 'symmetric')
  x <- seq(min(sen.cnv[i,]), max(sen.cnv[i,]), length.out = n)
  px <- predict(fit, newdata = x)
  z <- outer(px, px, '-')
  z <- z[lower.tri(z)]
  z <- sign(z)
  MS <- 2 * sum(z) / (n * (n - 1))
  lma <- lm(sen.exp[i,] ~ sen.cnv[i,])
  slope <- as.numeric(lma$coefficients[2])
  DSS <- MS * abs(slope)
  gene.DSS[i,2] <- DSS
}

# 输出
write.table(gene.DSS, 'DSS.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)
