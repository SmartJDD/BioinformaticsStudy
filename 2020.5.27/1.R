rm(list = ls())

library(splines)
library(survival)

survival.info <- as.matrix(read.csv('../第三次课/data/survival_info.txt',
                                    sep = '\t', header = T, row.names = 1))
gene.exp <- as.matrix(read.csv('../第三次课/data/gene_exp.txt',
                               sep = '\t', header = T))

# 建立生存对象
train.time <- as.numeric(survival.info[,1])
train.status <- as.numeric(survival.info[,2])
y<-Surv(train.time, train.status)

driver.gene <- c('ZNF286A', 'LMOD1', 'NRP2', 'APOLD1')
p <- matrix(0, nrow = 2, ncol = 4)
colnames(p) <- driver.gene
rownames(p) <- c('cox', 'log-rank')

for (i in driver.gene) {
  # 单因素cox回归
  # 衡量driver gene的表达对生存的影响
  x <- as.numeric(as.character(gene.exp[i,]))
  train <- list(train.time,train.status,x)
  q <- coxph((y ~ x), train)
  p[1,i] <- summary(q)$coefficients[,5]
  
  # KM生存曲线
  # 用于比较不同组之间存活时间的差异
  # 利用显著与预后相关的基因表达中值，将病人分成两组
  train.group <- matrix(0,2,ncol(gene.exp))
  train.group[1,] <- as.numeric(as.character(gene.exp[i,]))
  colnames(train.group) <- colnames(gene.exp)
  cutoff <- median(train.group[1,])
  train.group[2, which(train.group[1,] <= cutoff)] <- 1
  train.group[2, which(train.group[1,] > cutoff)] <- 2
  
  # 建立生存对象
  train <- list(train.time, train.status, train.group[2,])
  y <- Surv(train.time, train.status)
  dif <- survfit(y ~ train.group[2,], train)
  
  # 绘制生存曲线
  png(filename = paste(i, '.png', sep = ''))
  plot(dif)
  plot(dif, xlab="Time(days)",  ylab="Overall survival", mark=3,mark.time=T,col=c('red','blue'))
  dev.off()
  
  # log-rank检验
  r <- survdiff(y ~ train.group[2,], train) 
  p[2,i] <- 1 - pchisq(r$chisq, length(r$n) - 1)
}

# 多因素cox分析
clinical.info <- as.data.frame(read.csv('../第三次课/data/clinical_data.txt',
                                    sep = '\t', header = T, row.names = 1))
#建立生存对象
y <- Surv(as.numeric(clinical.info$days_to_death), clinical.info$status)
#注意这里的生存时间需要是数值型的数据
#多因素cox回归，性别、年龄对生存的影响
cox <- coxph(y ~ sex + age, data = clinical.info)
summary(cox)
