library(GenomicFeatures)
library(rtracklayer)
library(AnnotationHub)
BiocManager::install("AnnotationHub")
getwd()
setwd("D:\\data\\CTEPH组织转录组学\\publicdata\\GSE84538")
hub <- AnnotationHub()
query(hub, c("Homo sapiens", "GRCh38", "84"))
# 加载GTF文件
gtf_path <- "Homo_sapiens.GRCh38.84.gtf"
gtf <- rtracklayer::import(gtf_path)

# 提取基因长度
txdb <- makeTxDbFromGFF(gtf_path, format="gtf")
genes <- genes(txdb)
# 获取长度信息
geneLengths <- width(genes)
modulegen<-comdata[genes,]
#计算FPKM
library(edgeR)
gene<-read.csv("gene.csv",header = T)
genes<-gene$gene
rownames(data1)<-data1$ID
rownames(fpkm)<-fpkm$X
data1<-data1[,-1]
fpkm<-fpkm[,-1]
fpkm<-read.csv("fpkmdata.csv",header = T)
gene<-rownames(count)
rownames(count)<-count$X
count<-count[,-1]
# 假设counts是你的读数计数矩阵，行为基因，列为样本
# 假设geneLengths是一个与counts的行对应的向量，包含了相应的基因长度

# 首先，确保geneLengths是以千为单位的基因长度（kb）
geneLengthsKb <- geneLengths / 1000

# 计算每个样本的总映射读数（百万为单位）
totalMappedReads <- colSums(count) / 1e6

# 初始化FPKM矩阵
fpkm <- count

# 对每个样本进行循环，计算FPKM值
for (i in 1:ncol(count)) {
  fpkm[, i] <- (count[, i] / geneLengthsKb) / totalMappedReads[i]
}

# 现在fpkm变量包含了每个样本的FPKM值
print(fpkm)
write.csv(fpkm,file = "fpkmdata.csv")
#ID转换
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
gene<-rownames(count)

mappedSymbols <-AnnotationDbi::select(org.Hs.eg.db, 
                                     keys = gene, 
                                     columns = c("SYMBOL"), 
                                    keytype = "ENSEMBL")

print(mappedSymbols)
#匹配
# 检查并移除任何可能的重复基因ID，因为一个基因ID可能对应多个基因符号
uniqueSymbols <- mappedSymbols[!duplicated(mappedSymbols$ENSEMBL), ]
# 检查并移除任何可能的重复基因symbol，保留第一个映射的symbol
uniqueSymbols <- uniqueSymbols %>%
  distinct(SYMBOL, .keep_all = TRUE)
# 去除基因符号为NA的行
cleanSymbols <- uniqueSymbols[!is.na(uniqueSymbols$SYMBOL), ]
# 使用cleanSymbols更新表达矩阵的行名
# 首先，确保表达矩阵中的基因ID与cleanSymbols中的基因ID完全匹配
# 这一步将会移除表达矩阵中那些无法在cleanSymbols中找到匹配基因符号的行（基因）
exprMatrixClean <- count[rownames(count) %in% cleanSymbols$ENSEMBL, ]
# 更新行名为基因符号
rownames(exprMatrixClean ) <- cleanSymbols$SYMBOL[match(rownames(exprMatrixClean), cleanSymbols$ENSEMBL)]
count1<-exprMatrixClean
#清除低表达数据
cutoff = 0.5
count1<- data.frame(count1[which(apply(count1, 1, function(x){length(which
                                                                      (x!= 0))/length(x)}) >= cutoff),])
fpkm1<-data.frame(fpkm1[which(apply(fpkm1, 1, function(x){length(which
                                                                 (x!= 0))/length(x)}) >= cutoff),])
#读入自己表达矩阵
# 确保exprMatrix2的行顺序与exprMatrix1相匹配
colnames(count1)<-c(paste("Control",c(6:9),sep = ""),paste("CTEPH",C(6:8),sep = ""))
count2 <- count2[match(rownames(count1), rownames(count2)),]
# 使用cbind()函数合并矩阵
combinedMatrix <- cbind(count1, count2)
combinedMatrix<-na.omit(combinedMatrix)


fpkm2 <- fpkm2[match(rownames(fpkm1), rownames(fpkm2)),]
colnames(fpkm1)<-c(paste("Control",c(6:9),sep = ""),paste("CTEPH",C(6:8),sep = ""))
# 使用cbind()函数合并矩阵
all_fpkm <- cbind(fpkm1, fpkm2)
all_fpkm<-na.omit(all_fpkm)
log_fpkm <- log2(all_fpkm + 1)
log_fpkm<-na.omit(log_fpkm)
# 进行批次效应校正
batch_corrected_fpkm <- ComBat(dat = log_fpkm, batch = batchVector, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
min_positive_value <- min(batch_corrected_fpkm[batch_corrected_fpkm > 0])
batch_corrected_fpkm[batch_corrected_fpkm < 0] <- min_positive_value
# 创建批次向量
# 假设exprMatrix1有m列，exprMatrix2有n列
batchVector <- c(rep("Batch1", ncol(count1)), rep("Batch2", ncol(count2)))
library(limma)
library(sva)
library(edgeR)
combinedMatrix<-na.omit(combinedMatrix)

# 创建DGEList对象
dge <- DGEList(counts = combinedMatrix)

# 使用TMM归一化（也有助于后续批次效应处理）
dge <- calcNormFactors(dge)

# 获取归一化后的CPM值
normMatrix <- cpm(dge, normalized.lib.sizes=TRUE)
#log转化
normMatrix<-log2(normMatrix + 1)
comdata<-batchCorrectedMatrix
comfpkm<-batch_corrected_fpkm
write.csv(comdata,file="comcount.csv")
write.csv(comfpkm,file="comfpkm.csv")
batchCorrectedMatrix <- ComBat(dat=normMatrix, batch=batchVector, mod=NULL, par.prior=T, prior.plots=FALSE)
min_positive_value <- min(batchCorrectedMatrix[batchCorrectedMatrix > 0])
batchCorrectedMatrix[batchCorrectedMatrix < 0] <- min_positive_value
# 对整合前的数据进行PCA
pcaBefore <- prcomp(t(log2(normMatrix+1)), scale. = TRUE)

# 对整合并去除批次效应后的数据进行PCA
pcaAfter <- prcomp(t(log2(batchCorrectedMatrix+1)), scale. = TRUE)
library(ggplot2)
comdata<-batchCorrectedMatrix
pca<-prcomp(t(data1), scale. = TRUE)
# 整合前的PCA结果绘制
dfBefore <- data.frame(PC1 = pcaBefore$x[,1], PC2 = pcaBefore$x[,2], Batch = batchVector)
pBefore <- ggplot(dfBefore, aes(x = PC1, y = PC2, color = as.factor(Batch))) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  labs(x = "PC1", y = "PC2", title = "PCA Before Batch Correction") +
  coord_fixed()

# 去除批次效应后的PCA结果绘制
library(ggthemes)  # 额外的美化主题和比例尺
library(ggpubr)
color_values <- c(Control = "#377eb8", CTEPH = "#e41a1c")
fill_values <- c("Batch1" = "#A6CEE3", "Batch2" = "#FB9A99")  # 填充颜色稍微透明
# 定义Nature风格的颜色
color_palette <- c("#377eb8", "#e41a1c")
dfBefore <- data.frame(
  Dim1 = pca$x[,1],
  Dim2 = pca$x[,2],
  Condition =c(rep("Control",3), rep("CTEPH", 5)))

dfAfter <- data.frame(
  Dim1 = pcaAfter$x[,1],
  Dim2 = pcaAfter$x[,2],
  Condition =c(rep("Batch1", ncol(count1)), rep("Batch2", ncol(count2))))
group_list <- factor(c(rep("Control",9),rep("CTEPH",8)),levels = c("Control","CTEPH"))  ## 随机定义分组
draw_pca(comdata,
         group_list,
         color = c("#377eb8","#e41a1c"),
         addEllipses =T,
         style = 'ggplot2'#3D创建
)+ggprism::theme_prism()
p <- ggplot(dfBefore, aes(x = Dim1, y = Dim2, color = Condition, fill = Condition)) +
  geom_point(shape = 21, size = 3) +  # 实心点
  scale_color_manual(values = c("Control" = "#377eb8", "CTEPH" = "#e41a1c")) +
  scale_fill_manual(values = c("Control" = "#377eb8", "CTEPH" = "#e41a1c")) +
  stat_ellipse(geom = "polygon", alpha = 0.2, linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "Dim.1", y = "Dim.2", title = "PCA Plot", color = "Condition", fill = "Condition") +
  coord_fixed()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
print(pafter)
colnames(batchCorrectedMatrix)
comfpkm<-comfpkm[,c(8:12,1:4,13:17,5:7)]
comdata<-comdata[,c(8:12,1:4,13:17,5:7)]
colnames(comdata)
#差异分析
condition = factor(c(rep("Ctrl",9),rep("PE",8)),levels = c("Ctrl", "PE"))
design <- model.matrix(~0+factor(condition))
names<-colnames(comdata)
rownames(design)<-names
colnames(design) <- c("Control","CTEPH")
fit <- lmFit(comdata,design)
cont.matrix<-makeContrasts(CTEPH-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
DEGs2<-na.omit(allDiff)
diff2<- DEGs2[with(DEGs2, (abs(DEGs2$logFC)> 1 & adj.P.Val< 0.05 )), ]
colnames(diff2)[2]<-"logFC"
#6、数据清洗####
#1、上下调 
tb.2 <- diff2 %>%
  mutate(Disease = ifelse(logFC > 0, "Up","Down"))%>%
  mutate(Gene.symbol=rownames(diff2))
table(tb.2$Disease)
write.csv(tb.2,file = "genediff.csv")
#write table
#3. 分别提取上下调的基因
CTEPHup <- tb.2$symbol[tb.2$Disease == "Up"]
CTEPHdown <- tb.2$symbol[tb.2$Disease == "Down"]
####二、火山图####
##1、加载包 ----
library(tidyverse)
library(ggrepel)#给基因添加名字注释，让名字注释不重叠
colnames(diff1)[3]<-"logFC"
##2、准备火山图数据 ####
allDiff2<-diff1
DEG_all2 <- allDiff %>% 
  rownames_to_column("symbol") %>% #行名变成一列
  dplyr::select(symbol,logFC=logFC, Pvalue =adj.P.Val) %>% #选择三列，select可以在选取的同时重命名
  mutate(direction = factor(ifelse(Pvalue < 0.05 & abs(logFC) >1,#添加direction一列
                                   ifelse(logFC > 0, "Up",
                                          "Down"),
                                   "NS"),
                            levels=c('Up','Down','NS')
  )
  )#这个代码嵌套太多，比较复杂，注意理清楚逻辑
#### 绘制火山图 ----
library(ggplot2)
library(tidyverse)
library(ggrepel)#给基因添加名字注释，让名字注释不重叠
summary(DEG_all$logFC)
DEG_all<-na.omit(DEG_all)
ggplot(data = DEG_all, aes(x = logFC, y = -log10(Pvalue), colour = direction)) + #数据映射
  geom_point(size=3,alpha = 1) +    #散点图，alpha就是点的透明度
  theme_bw() +#设定主题 
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +         #不想要direction这个title
  theme(text=element_text(size=10,  family="sans"))+
  ylab(expression(-log[10]("P Value"))) +     #expression的作用就是让log10的10下标
  xlab(expression(log[2]("Fold Change"))) +
  xlim(- 4, 4) +
  ylim(0, 8) +
  geom_vline(xintercept = c(-0.585, 0.585), #加垂直线，在-1和+1之间划线
             lty = 2,
             col = "black",
             lwd = 0.6) +
  geom_hline(yintercept = -log10(0.05),
             lty = 2,
             col = "black",
             lwd = 0.6) +theme(text = element_text(size = 10))+ #标记点
  scale_color_manual(values = c("#f00000", "#1d2671", "#808080")) + #手动调颜色
  geom_text_repel(data = DEG_all %>% filter(Pvalue < 0.05, abs(logFC) >10),#加注释，筛选出差异非常显著的基因
                  aes(label = symbol),#便签就是基因名 geom_text_repel：adds text directly to the plot.
                  size =6,
                  segment.color = "black", #连接线的颜色，就是名字和点之间的线
                  show.legend = FALSE)+theme(legend.text = element_text(size = 15))+theme(axis.text = element_text(size = 15))#坐标轴

##绘制热图####
#### 准备热图数据 ----
Up<-DEG_all[which(DEG_all$direction=="Up"),]
Down<-DEG_all[which(DEG_all$direction=="Down"),]
#取出上下调的表达矩阵
Updata <-comdata[Up$symbol,]
Downdata <- comdata[Down$symbol,]
# 绘制热图 
#### 准备热图数据 ---
#1、按照表达差异倍数取
genes_sig1 <-tb.2[order(tb.2$logFC,decreasing = T),]%>% #按照表达矩阵从
  dplyr::slice(1:20) %>% #选择前50个p值最小的基因
  pull(symbol)#拉出Symbol作为一个向量

genes_sig2 <- Down[order(Down$Pvalue,decreasing = T),]%>% #按照表达矩阵从
  dplyr::slice(1:20) %>% #选择前50个p值最小的基因
  pull(symbol)#拉出Symbol作为一个向量
genes_sig2 <- Down %>% arrange(logFC) %>% #按照Pvalue从小到大进行排序
  dplyr::slice(1:25) %>% #选择前50个p值最小的基因
  pull(symbol)
#2、按照P值
genes_sig2 <- tb.2 %>% arrange(P.Value) %>% #按照Pvalue从小到大进行排序
  dplyr::slice(1:40) %>% #选择前50个p值最小的基因
  pull(Gene.symbol)#拉出Symbol作为一个向量
genes_sig2 <- tb.2 %>% arrange(P.Value) %>% #按照Pvalue从小到大进行排序
  dplyr::slice(1:10) %>% #选择前50个p值最小的基因
  pull(symbol)
#提取表达矩阵
DATA<-tb.2[genes,]
rownames(exp)<-exp$symbol
genes_sig<-tb.2$symbol
genes_sig <-  rbind(genes_sig1,genes_sig2)
heatmap_df4 <- as.data.frame(data)#筛选出这50个基因的表达矩阵
heatmap_df4 <- as.data.frame(comdata[genes_sig2, ])#筛选出这50个基因的表达矩阵
colnames(heatmap_df1)
heat_df2<-heatmap_df4[,c(1:15,28:64)]
#### 绘制热图 ----
library(graphics)
heat_df1<-scale(t(heatmap_df4),scale = T,center = T)
head(heat_df1)
heat_df1<-as.data.frame(scale(heatmap_df4))
heat_df1<-as.data.frame(scale(modulegen))
modulegen<-t(modulegen)
#构建注释信息
#构建列注释信息
annotation_row = data.frame(
  Group = factor(rep(c("Control", "CTEPH"), c(9, 8)),
                 levels =c("Control", "CTEPH"))
)
rownames(annotation_row) = rownames(modulegen)
#构建行注释
# 构建行注释信息
annotation_col = data.frame(
  GeneClass = factor(rep(c("Up", "Down"), c(10, 10)))
)
rownames(annotation_col) = genes
head(annotation_row)
# 自定注释信息的颜色列表
#ffe000','blue'
gene<-c()
ann_colors = list(
  Group = c(Control = "#2F7FC1",CTEPH="#D8383A"),
  GeneClass = c(Up = "#7570B3", Down = "#66A61E")
)
pdf(file="heatmap1.pdf")
mycols <- colorRamp2(breaks = c(-2, 0, 2),
                     colors = colorRampPalette(c("#0000EF", "white", "red")))
heat1<-ann[c("KDR","VEGFA"),]
library(pheatmap)
pheatmap(modulegen, 
         scale = "column",#每个基因的水平scale一下，只对行进行转换
         border_color = "black",
         cellwidth = 16,
         cellheight = 8,
         cluster_rows = F,
         cluster_cols =F,
         #main = "Example heatmap",#标题
         annotation_row = annotation_row,#注释
         annotation_colors = ann_colors,
         annotation_legend = T,#注释图
         #treeheight_row = 10,
         #treeheight_col = 12,
         show_colnames =T,
         show_rownames = F,
         annotation_names_row = T, 
         annotation_names_col = T,
         density.info = "none",
         #cutree_rows = 2,#聚类分割
         #gaps_col = c(9),
         color = colorRampPalette(c("#1AA3FF","white","#FF6B67"))(100),
         fontsize_row = 12)
graph2pdf(file = 'heatmap1.pdf',width = 10,height =7)
write.csv(ann,file="allexp.csv")
####读入和读出数据####
gene<-tb.2$symbol
heatmap_df4<-ann1[gene,]
write.csv(heatmap_df4,file="finaldata.csv")
write.table(heatmap_df4,file = "Top25DEGs.txt",
            sep = "\t",row.names = T,col.names = T,
            quote = F)

data <- read.table(file = "mRNA/rawdata.csv",
                   sep = "\t",
                   header = T,
                   check.names = F,
                   quote = "")

#4.写出到excel中
library(openxlsx)
deg<-ann[rownames(diffSig),]
heatmap_df4 <- heatmap_df4 %>% as.data.frame() %>%
  rownames_to_column("symbol")
ann2<-as.matrix(ann1)
ann2$symbol<-as.vector(rownames(ann2))
Downdata$symbol<-rownames(Downdata)
out.data2 <- list(G1data = ann1,
                  G2data=Updata,
                  G3data=Downdata,
                  G4data=exp1
                  
)
save
outXlsx <- "ALLdata.xlsx"
write.xlsx(out.data2,file =outXlsx)

#ssGSEA分析
清空变量
remove(list = ls())
comfpkm<-batch_corrected_fpkm
comfpkm<-comfpkm[,c(8:12,1:4,13:17,5:7)]
# 读取标准化后的表达矩阵
input <- read.table("data/normalize.txt", sep="\t",
                    header=T, row.names = 1, check.names=F)#这是一个经过标准化后的芯片表达矩阵
exp <- as.matrix(input)

# 加载所需R包

if(!"XML"%in%installed.packages()) 
  install.packages("XML", type = "binary")
library(GSEABase)
library(GSVA)

# 设置gmt文件路径，这次我们用MsigDb（分子标签数据库）数据集
getwd()
msigdb_GMTs <- "h.all.v2023.2.Hs.symbols.gmt"
msigdb <- "h.all.v2023.2.Hs.symbols.gmt"#定义相关基因集

# 获取gmt文件内容
all_msigdb <- read.gmt(file.path(msigdb))
geneset2 <- getGmt(file.path(msigdb)) #作为GSVA的功能数据集 
# 进行GSVA分析
class(fpkm2)
fpkm2<-as.matrix(log2(fpkm2+1))
data1<-as.matrix(log2(data1+1))
data1<-na.omit(data1)
gsva_data <- gsva(data1, geneset2, 
                  method="ssgsea",
                  mx.diff=FALSE, verbose=FALSE, 
                  parallel.sz=1)#线程数

# 创建临床表型
Type=c(rep("Control",3),rep("CTEPH",5))#和表达矩阵样本队对应
clin<-as.data.frame(cbind(as.vector(colnames(data1)),Type))
#合并信息
gsva_data <-t(gsva_data)
datafpkm1<-cbind(clin,gsva_data)

datafpkm1=na.omit(datafpkm1)
gene_columns <- names(datafpkm1)[-c(1,2)]
datafpkm1<-datafpkm1[,-1]
# 使用aggregate函数进行分组和求均值操作
resultfpkm1 <- do.call(data.frame, aggregate(. ~ Type, data = datafpkm1, FUN = function(x) mean(x, na.rm = TRUE)))
#基因集转换
strings <- gene_columns
# 创建一个空的向量保存处理后的结果
result1 <- vector("list", length(strings))
# 对每个字符串进行处理
for (i in 1:length(strings)) {
  splitted <- strsplit(strings[i], "_", fixed = TRUE)[[1]]
  rest_part <- paste(splitted[-1], collapse = "_")
  rest_part_lower <- tolower(rest_part)
  rest_part_title <- paste(toupper(substring(rest_part_lower, 1, 1)), substring(rest_part_lower, 2), sep = "")
  result1[[i]] <- rest_part_title
}
# 输出处理后的结果
print(result1)
result_string <- paste(result1, collapse = " ")
gene_names <- unlist(strsplit(result_string, " ")) 
# 输出转换后的字符串
print(result_string)
resultfpkm1<-resultfpkm1[,-1]
datafpkm1<-datafpkm1[,-1]
colnames(datafpkm1)<-gene_names
rownames(resultfpkm1)<-c("Control","CTEPH")
colnames(resultfpkm1)<-gene_names
resultfpkm1<-t(resultfpkm1)
#构建表达矩阵
#构建注释信息
#构建列注释信息
annotation_col = data.frame(
  Group = factor(rep(c("Control", "CTEPH"), c(3, 5)),
                 levels = c("Control","CTEPH"))
)
datafpkm1<-t(datafpkm1)
rownames(annotation_col) = colnames(datafpkm1)
#构建行注释
# 构建行注释信息
annotation_row = data.frame(
  GeneClass = factor(rep(c("Up", "Down"), c(25, 25)))
)
rownames(annotation_row) = rownames(resultfpkm1)
head(annotation_row)
# 自定注释信息的颜色列表
#ffe000','blue'
ann_colors = list(
  Group = c(Control = "blue", CTEPH = "red"),
  GeneClass = c(Up = "#7570B3", Down = "#66A61E")
)
library(pheatmap)
pheatmap(datafpkm1, 
         scale = "row",#每个基因的水平scale一下，只对行进行转换
         border_color = "black",
         cellwidth = 8,
         cellheight = 12,
         cluster_rows = T,
         cluster_cols =F,
         #main = "Example heatmap",#标题
         annotation_col = annotation_col,#注释
         annotation_colors = ann_colors,
         annotation_legend = T,#注释图
         #treeheight_row = 10,
         #treeheight_col = 12,
         show_colnames =F,
         show_rownames = T,
         annotation_names_row = T, 
         annotation_names_col = T,
         density.info = "none",
         cutree_rows = 3,#聚类分割
         gaps_col = c(3),
         color = colorRampPalette(c("#1AA3FF","white","#FF6B67"))(100),
         fontsize_row = 12)
save(datafpkm,file="SSGSEA.csv")
### 火山图
annotation_col = data.frame(
  Group = factor(rep(c("Control", "CTEPH"), c(1, 1)),
                 levels = c("Control","CTEPH"))
)
data<-t(data)
rownames(annotation_col) = colnames(resultfpkm)
#构建行注释
# 构建行注释信息
annotation_row = data.frame(
  GeneClass = factor(rep(c("Up", "Down"), c(25, 25)))
)
rownames(annotation_row) = rownames(resultfpkm)
head(annotation_row)
# 自定注释信息的颜色列表
#ffe000','blue'
ann_colors = list(
  Group = c(Control = "blue", CTEPH = "red"),
  GeneClass = c(Up = "#7570B3", Down = "#66A61E")
)
library(pheatmap)
pheatmap(resultfpkm, 
         scale = "row",#每个基因的水平scale一下，只对行进行转换
         border_color = "black",
         cellwidth = 16,
         cellheight = 12,
         cluster_rows = T,
         cluster_cols =F,
         #main = "Example heatmap",#标题
         annotation_col = annotation_col,#注释
         annotation_colors = ann_colors,
         annotation_legend = T,#注释图
         #treeheight_row = 10,
         #treeheight_col = 12,
         show_colnames =F,
         show_rownames = T,
         annotation_names_row = T, 
         annotation_names_col = T,
         density.info = "none",
         cutree_rows = 2,#聚类分割
         gaps_col = c(1),
         color =  colorRampPalette(c("#1AA3FF","white","#FF6B67"))(100),
         fontsize_row = 12)

















library(pheatmap)
GSVA<-t(gsva_data)
pheatmap(GSVA, #热图的数据
         cluster_rows = T,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         ##annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         ##scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 10, cellheight = 8,# 格子比例
         fontsize = 10)
#选择特定通路，小提琴图，组间比较
data1$HALLMARK_INFLAMMATORY_RESPONSE
data2<-data1[,c("V1","Type","HALLMARK_INFLAMMATORY_RESPONSE")]
colnames(data2)=c("id","Type","Expression")
rt<-data2
#设置比较租
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
my_comparisons=list(c("Control","CTEPH"))
library(ggplot2)
#绘制boxplot
scores_data <- data.frame(VEGFA = rt$Expression,
                          group = rt$Type)
col1<-cols
cols<-c(Control="blue",CTEPH="red")
MG1<-ggplot(scores_data,
            aes(x=group, y= VEGFA,fill=group))+
  geom_violin(aes(fill=group), scale = "area", trim=F, size=0.5)+
  geom_violin(scale = "area", trim=F, size=0.5,color="white",cex=1,alpha=1)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="white")+
  #facet_grid(.~cell_types_2_groups)+                     
  xlab("")+
  ylab("Hypoxia score")+
  scale_fill_manual(values=cols)+
  theme_classic(base_size=14)+
  ggsignif::geom_signif(comparisons = list(c("Control", "CTEPH")), # 第2组与第3组的比较
                        map_signif_level = T)+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))


#----------------------------clusterProfiler------------------------------------------
setwd("E:/生物信息学/GSEA分析")
library(org.Hs.eg.db)
library(clusterProfiler)
diff<-allDiff
diff <- read.csv("diff.csv", header = T, row.names = 1)

genesymbol <- rownames(diff)
entrezID <- bitr(genesymbol,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")#genesymbol转化为ID

genelist <- diff$logFC
names(genelist) <- rownames(diff)
genelist <- genelist[names(genelist) %in% entrezID[,1]]
names(genelist) <- entrezID[match(names(genelist),entrezID[,1]),2]
genelist <- sort(genelist,decreasing = T)

R.utils::setOption( "clusterProfiler.download.method",'auto' )

#GSEA分析：
KEGG_gesa <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)
#提取结果，将结果中的ID转化为genesymbol
KEGG_ges_result = setReadable(KEGG_gesa,
                              OrgDb = "org.Hs.eg.db",
                              keyType = "ENTREZID")
KEGG_ges_result <- KEGG_ges_result@result


GO_kk_entrez <- gseGO(geneList     = genelist,
                      ont          = "BP",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = org.Hs.eg.db,#人类org.Hs.eg.db 鼠org.Mm.eg.db
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.25)   #实际为padj阈值可调整


GO_ges_result = setReadable(GO_kk_entrez ,
                              OrgDb = "org.Hs.eg.db",
                              keyType = "ENTREZID")
GO_ges_result <- GO_ges_result@result



write.csv(GO_ges_result,file="gsea.GO.csv")



#可视化
library(ggplot2)
library(enrichplot)
gseaplot2(KEGG_gesa,
          geneSetID = 2,
          color = "#93C15A",
          rel_heights = c(1.5, 0.3, 1),
          subplots = 1:3,
          pvalue_table = F,
          ES_geom = "line",
          title=KEGG_gesa$Description[2])



gseaplot2(KEGG_gesa,
          geneSetID = c(1,3,5),
          color = c("#93C15A","red","purple"),
          rel_heights = c(1.5, 0.3, 1),
          subplots = 1:2,
          pvalue_table = F,
          ES_geom = "line")



#----------------------------fgsea------------------------------------------
setwd("E:/生物信息学/GSEA分析")
diff <- read.csv("diff.csv", header = T, row.names = 1)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)


#human <- msigdbr(species = "Homo sapiens")
#table(human$gs_subcat)

## 预定义基因集
geneset_KEGG = msigdbr(species = "Homo sapiens",#mouse:Mus musculus
                       category = "C2", 
                       subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#去除前缀KEGG_
geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)#将大写换为小写
geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#将_转化为空格
library(Hmisc)
geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)#首字母大写
GSEA_geneset_KEGG <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)


## 根据logfc降序排列基因
df <- diff[order(diff$avg_log2FC,decreasing = T),]
ranks <- df$avg_log2FC
names(ranks) <- rownames(df)

## GSEA分析
GSEA_df <- fgsea(pathways = GSEA_geneset_KEGG, 
                 stats = ranks,
                 minSize=10,
                 maxSize=500,
                 eps=0.0)
library(data.table)
fwrite(GSEA_df, file="GSEA_df.txt", sep="\t", sep2=c("", " ", ""))


#可视化
library(ggplot2)
plotEnrichment(GSEA_geneset_KEGG[["Alzheimers disease"]],
               ranks,ticksSize = 0.5) + labs(title="Alzheimers disease")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))




#可视化
ToPlot  <- sapply(GSEA_df$pathway, function(Pathway){
  pathway <- GSEA_geneset_KEGG[[Pathway]]
  stats <- ranks
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^1)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0), pathway = Pathway)
  return(list(toPlot))
})

gseplotdata <- ToPlot$`Alzheimers disease`

p1 <- ggplot(gseplotdata, aes(x = x, y = y)) + 
  geom_line(size = 1, show.legend = FALSE,color="darkgreen") + 
  geom_hline(yintercept = 0, linetype=2, lwd = 1) +
  ylab("Enrichment Score") +
  xlab('')+
  labs(title="Alzheimers disease")+
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 12,colour = 'black'),
        legend.position="none",
        plot.title = element_text(size=15,hjust =0.1, face = 'bold'),
        axis.ticks.x = element_blank(),
        plot.margin=margin(t=.2, r = .2, b= -0.5, l=.2, unit="cm"))



test <- gseplotdata[1:floor(nrow(gseplotdata)/length(1)),]
test$xend = test$x+1
p2 <- ggplot(test)+ 
  geom_rect(aes(xmin = x,xmax = xend , ymin = 0 , ymax = 1, fill=x), 
            color="black", size=0.5)+
  theme_bw() +
  theme(legend.position = "none",
        plot.margin = margin(t=0, b=0,unit="cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand=c(0,0))



p3 <- ggplot(test)+ 
  geom_rect(aes(xmin = x-50,xmax = xend+300, 
                ymin = 0 , ymax = 1, fill=x), lwd=4)+
  scale_fill_gradientn(colours = c("#035BFD", "#397EFC", "#5B94FB","white",
                                   "#F77A7C", "#F45557", "#FB0407"))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  
  theme(legend.position = "none",
        plot.margin = margin(t=0, b=0,unit="cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank()) +
  scale_y_continuous(expand=c(0,0))



rel_heights <- c(1, .15, .1)
plotlist <- list(p1, p2, p3)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(size = 12, face = "bold"))


library(cowplot)
plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)


###############################
ggplot(gseplotdata, aes(x = x, y = y)) + 
  geom_point(fill="darkgreen", shape=21, size=5) + 
  geom_hline(yintercept = 0, linetype=2, lwd = 1) +
  ylab("Enrichment Score") +
  xlab('')+
  labs(title="Alzheimers disease")+
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 15,colour = 'black'),
        legend.position="none",
        plot.title = element_text(size=15,hjust =0.1, face = 'bold'),
        axis.ticks.x = element_blank())




