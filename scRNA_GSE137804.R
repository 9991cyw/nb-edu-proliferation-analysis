#参考教程：
#https://www.jianshu.com/p/159df5c4ff21

library(dbplyr) #保证这个包优先加载
library(harmony)
library(tidyverse)
library(patchwork)
library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(DoubletFinder)
library(qs)
setwd("~/scRNA/R_Project/Schwann_project/Neuroblastoma/")
#GSE253903,GSE179159,GSE155512,GSE131778动脉粥样硬化单细胞分析
#设置好本项目专用随机种子Seed：240701292
set.seed(240701292)

###################----一、 单细胞数据下载、读取及以及质控----##################
#手动将下载的10X文件按样本分文件夹，每个文件夹内部包含三个文件
#"barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz"

### 1. GEO手动下载数据集并处理 ### 
untar(tarfile = "./GSE137804_RAW.tar", exdir = "./GSE137804/")

###读取GSE137804样本中的神经母细胞瘤数据
setwd("~/scRNA/R_Project/Schwann_project/Neuroblastoma/GSE137804/")
assays<-dir("~/scRNA/R_Project/Schwann_project/Neuroblastoma/GSE137804/")
count_files <- grep(".xls.gz$", assays, value = TRUE)
# 假设 project 和 list 都是你已定义的列表
project <- str_sub(count_files,1,10)
scelist<-list()
for (i in 1:length(count_files)) {
  count_data<-fread(count_files[[1]])
  count_data <- count_data[!duplicated(count_data$Symbol), ]
  count_data <- column_to_rownames(count_data,var = "Symbol")
  count_data <- count_data[, -which(names(count_data) == "Gene_ID")]
  count_data <- as.matrix(count_data)
  
  # 将数据转换为稀疏矩阵格式
  counts <- as(count_data, "sparseMatrix")
  
  # 创建 Seurat 对象，并设置 orig.ident 为 project[i]
  scelist[[i]] <- CreateSeuratObject(counts = counts, project = project[[i]],
                                     min.cells = 3, min.features = 200)
  
  # 给细胞 barcode 加个前缀，防止合并后 barcode 重名
  scelist[[i]] <- RenameCells(scelist[[i]], add.cell.id = project[[i]])
  
  # 计算线粒体基因比例
  scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^MT-")
  
  # 计算核糖体基因比例
  scelist[[i]][["percent.rb"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^RP[SL]")
  
  # 计算红细胞基因比例
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(scelist[[i]]))
  scelist[[i]][["percent.HB"]] <- PercentageFeatureSet(scelist[[i]], features = HB.genes)
}
###数据15-22的表达矩阵和之前不一致，更改读取
for (i in 15:22) {
  count_data<-fread(count_files[[i]])
  count_data <- count_data[!duplicated(count_data$V1), ]
  count_data <- column_to_rownames(count_data,var = "V1")
  count_data <- as.matrix(count_data)
  
  # 将数据转换为稀疏矩阵格式
  counts <- as(count_data, "sparseMatrix")
  
  # 创建 Seurat 对象，并设置 orig.ident 为 project[i]
  scelist[[i]] <- CreateSeuratObject(counts = counts, project = project[[i]],
                                     min.cells = 3, min.features = 200)
  
  # 给细胞 barcode 加个前缀，防止合并后 barcode 重名
  scelist[[i]] <- RenameCells(scelist[[i]], add.cell.id = project[[i]])
  
  # 计算线粒体基因比例
  scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^MT-")
  
  # 计算核糖体基因比例
  scelist[[i]][["percent.rb"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^RP[SL]")
  
  # 计算红细胞基因比例
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(scelist[[i]]))
  scelist[[i]][["percent.HB"]] <- PercentageFeatureSet(scelist[[i]], features = HB.genes)
}
# 打印每个Seurat对象的信息
for (name in names(scelist)) {
  cat("File:", name, "\n")
  print(scelist[[name]])
  cat("\n")
}
for (i in 1:22) {
  names(scelist)[[i]]<-project[[i]]
}
###去除非神经母细胞瘤的数据，这六个数据包括2个胚胎样本及4份胎儿肾上腺样本
scelist<-scelist[-c(11:16)]
names(scelist)

dir.create("~/scRNA/R_Project/Schwann_project/Integrate/GSE137804")
setwd("~/scRNA/R_Project/Schwann_project/Integrate/GSE137804/")
system.time(qsave(scelist,file = "scRNAlist_NB_only.qs"))
rm(list = ls())
gc()
##使用Merge函数将scRNAlist合并为一个Seurat对象
scRNAlist<-qread("./scRNAlist_NB_only.qs")

scRNA <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA
table(scRNA$orig.ident)
### 2. 数据质控 ###
##绘制质控小提琴图
#设置可能用到的主题
theme_set2 = theme(axis.line.x = element_blank())
#设置绘图元素
plot.features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb","percent.HB")
group = "orig.ident"
#质控前小提琴图
plots=list()
for (i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA,group.by = group,pt.size = 0,raster=FALSE,
                       features = plot.features[i]) + theme_set2 + NoLegend()
}
violin<-wrap_plots(plots = plots,nrow=2)
dir.create("QC_NB_GSE137804")
ggsave("QC_NB_GSE137804/vlnPlot_before_qc.pdf",plot = violin, width = 16, height = 8)

#设置质控标准
#部分数据太夸张，看不出来需要个性化展示
hist(scRNA@meta.data[["nCount_RNA"]][scRNA@meta.data[["nCount_RNA"]] < 40000],breaks = c(0,5000,10000,15000,20000,25000,30000,35000,40000))
#UMI(Count)设置为2W是没问题的
minGene = 500
maxGene = 6000
maxUMI = 20000
pctMT = 15
pctHB = 1
##数据质控并绘制小提琴图
scRNA <- subset(scRNA,subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & nFeature_RNA < maxGene &
                  percent.mt < pctMT & percent.HB < pctHB)
###质控之前： 30256 features across 179283 sample
###质控之后： 30256 features across 166992 samples
table(scRNA$orig.ident)
#这一步去掉了低质量样本GSM4710703_PDAC_TISSUE_14
plots = list()
for (i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA,group.by = group,pt.size = 0,raster=FALSE,
                       features = plot.features[i]) + theme_set2 + NoLegend()
}
violin<-wrap_plots(plots = plots,nrow=2)
ggsave("QC_NB_GSE137804/vlnPlot_after_qc.pdf",plot = violin, width = 16, height = 8)

#### 查看批次效应
#标准化、归一化、PCA
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA <- RunPCA(scRNA,seed.use = 240701292)
ElbowPlot(scRNA, ndims = 50)
#选择拐点对应的PC数
pc.num<-1:30
scRNA <- FindNeighbors(scRNA,dims = pc.num) %>% FindClusters()
scRNA <- scRNA %>% RunUMAP(dims = pc.num,seed.use = 240701292)
scRNA$Sample <- scRNA$orig.ident
DimPlot(scRNA,label=T,raster=FALSE,seed = 240701292)
DimPlot(scRNA,label=T,raster=FALSE,seed = 240701292,group.by = "Sample")
qsave(scRNA,file = "~/scRNA/R_Project/Schwann_project/Integrate/GSE137804/scRNA_unHarmonyed_AS.qs")
rm(list = ls())
gc()
### 3.--去除双细胞 ###
#拆分为单个seurat子对象
#1. 按样本进行拆分
scRNA <- qread("~/scRNA/R_Project/Schwann_project/Integrate/GSE137804/scRNA_unHarmonyed_AS.qs")
DefaultAssay(scRNA)<-"RNA"
scRNA.list <- SplitObject(scRNA, split.by = "Sample")
dim.usage<-30 #根据PCA肘图
##DoubletFinder去双胞的标准流程，封装成一个函数
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(0.075*ncol(data)) #github官网提供数值
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  return(data)
}

#循环对每个样本进行双细胞标注
for (i in 1:length(scRNA.list)) {
  print(i)
  scRNA.list[[i]] <- NormalizeData(scRNA.list[[i]], verbose = FALSE)
  scRNA.list[[i]] <- FindVariableFeatures(scRNA.list[[i]], selection.method = "vst", 
                                          nfeatures = 2000, verbose = FALSE)
  scRNA.list[[i]] <- ScaleData(scRNA.list[[i]])
  scRNA.list[[i]] <- RunPCA(scRNA.list[[i]],seed.use = 240701292)
  scRNA.list[[i]] <- RunUMAP(scRNA.list[[i]], dims = 1:30,seed.use = 240701292)
  scRNA.list[[i]] <- Find_doublet(scRNA.list[[i]])
}
scRNA <- merge(scRNA.list[[1]],scRNA.list[2:length(scRNA.list)])
scRNA.sg<-scRNA[,scRNA@meta.data$doublet_info %in% "Singlet"]
table(scRNA@meta.data$doublet_info)
##去除双细胞之前：32335 features across 323502 samples#这个是包含GSE26155的
#Doublet Singlet 
# 24259  299243  
#不包括GSE26155的
#Doublet Singlet 
#12524  154468
#保存数据集
qsave(scRNA.sg,file = "~/scRNA/R_Project/Schwann_project/Integrate/GSE137804/scRNA_Dedoub_UnInte_AS.qs")
rm(list = ls())
gc()



################---- 二、用Harmony对多样本数据集进行整合----####################
### 1. 单纯的单细胞数据整合推荐使用harmony ###

scRNA.sg <- qread("~/scRNA/R_Project/Schwann_project/Integrate/GSE137804/scRNA_Dedoub_UnInte_AS.qs")
names(scRNA.sg@meta.data)
cellinfo<-subset(scRNA.sg@meta.data,select = c("orig.ident","Sample","percent.mt","percent.rb","percent.HB"))
dim(cellinfo)
scRNA.sg <- CreateSeuratObject(scRNA.sg@assays$RNA@counts,meta.data = cellinfo)
#使用SCT标准化,运行后默认矩阵为SCT，改的话需要defaultAssay()<-"RNA"
scRNA.sg <- NormalizeData(scRNA.sg) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.sg <- RunPCA(scRNA.sg,seed.use = 240701292)
ElbowPlot(scRNA.sg, ndims = 50)
###PCA
#scRNA.sg <- RunPCA(scRNA.sg,npcs = 50,seed.use = 240701292,assay = "SCT")
## 整合方法：1. 单个样本间进行整合（推荐，效果更好）
scRNA.sg <- RunHarmony(scRNA.sg,group.by.vars = "Sample",#assay.use = "SCT",
                       max.iter.harmony = 10)
# lambda参数默认是1，唯一一个决定了Harmony的整合力度的参数，值越小，整合力度越大，取0.5-2
pc.num=1:30
###找cluster
scRNA.sg <- FindNeighbors(scRNA.sg,dims = pc.num,reduction = "harmony"#,assay = "SCT"
)#本次不用SCT
scRNA.sg <- FindClusters(scRNA.sg,resolution = seq(0,2,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.sg)
sel.clust = "RNA_snn_res.1"
scRNA.sg <- SetIdent(scRNA.sg, value = sel.clust)
table(scRNA.sg@active.ident)
scRNA.sg = RunUMAP(scRNA.sg,reduction = "harmony",dims = pc.num,random.seed = 240701292)
scRNA.sg$seurat_clusters <- scRNA.sg$RNA_snn_res.1
DimPlot(scRNA.sg,group.by = "Sample",reduction = "umap",raster=FALSE)
DimPlot(scRNA.sg,label = T,raster=FALSE)
##合并部分临床信息
#meta<-read.table('../metadata/GSE137804SraRunTable.txt',sep = ',',header = T)
library(readxl)
# 读取meta信息
NB_meta <- read_excel("../metadata/GSE137804meta.xlsx")
NB_meta <- NB_meta[,c(1,6,7)]
# 将Seurat对象的元数据提取出来
seurat_metadata <- scRNA.sg@meta.data
seurat_metadata <- rownames_to_column(seurat_metadata,var = "barcodes")
# 合并元数据，确保匹配样本
combined_metadata <- merge(seurat_metadata, NB_meta, by.x = "Sample", by.y = "Sample", all.x = TRUE)
combined_metadata[1:5,30:37]
combined_metadata<-column_to_rownames(combined_metadata,var = "barcodes")
combined_metadata<-combined_metadata[,c(2,1,3:33,35,36,34)]
combined_metadata<-combined_metadata[rownames(scRNA.sg@meta.data),]
# 将合并后的元数据重新赋值给Seurat对象
scRNA.sg@meta.data <- combined_metadata
# 验证元数据是否成功添加
head(scRNA.sg@meta.data)
# 保存修改后的Seurat对象
qsave(scRNA.sg,file = "~/scRNA/R_Project/Schwann_project/Integrate/GSE137804/scRNA_harmonyed_AS.qs")
rm(list = ls())
gc()


###################---- 三、对细胞类型进行逐级手工注释----######################
library(patchwork)

scRNA.sg<-qread("~/scRNA/R_Project/Schwann_project/Integrate/GSE137804/scRNA_harmonyed_AS.qs")
DimPlot(scRNA.sg,group.by = "RNA_snn_res.1",label = T,raster = FALSE)
### 1.先看看大类###
Markergene<-c("PTPRC",#免疫细胞
               "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
               "NRP1","CLEC4C","IRF7","IRF8",#DC
               "CD3D","CD3E","CD52","CD7",#T细胞
               "BANK1","MS4A1","CD79A","CD79B",#B细胞
               "MZB1","IGHG3","JCHAIN",#浆细胞
               "CMC1","GNLY","GZMB","NKG7",#NK细胞
               "DCN","COL1A1","COL3A1","TAGLN",#成纤维
               "CDH5","VWF","ENG","CLDN5",#内皮
               "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
               "SOX10","S100B","PLP1","EGR2"#施旺细胞
)
SCpubr::do_DotPlot(scRNA.sg,group.by = "RNA_snn_res.1",
                   assay = "RNA",scale = T,
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
DimPlot(scRNA.sg,group.by = "RNA_snn_res.1",label = T,raster=FALSE)+FeaturePlot(scRNA.sg,features = c("SOX10","EGR2"),raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC","CHGA","DCN","VWF","CEACAM5","PECAM1","EGR2"),cols = c("lightgrey","darkred"),raster=FALSE,
            pt.size = 0.1#, # Adjust point size for visibility
            #blend = TRUE
) 
table(scRNA.sg$RNA_snn_res.1)
##2.先将大类注释的基本类型注释好，手动赋值
table(scRNA.sg$seurat_clusters)
scRNA.sg$seurat_clusters<-scRNA.sg$RNA_snn_res.1

clusterAnn = as.data.frame(matrix(nrow=36,ncol=2),col.names=NULL ,row.names=NULL ) 
colnames(clusterAnn)<-c("ClusterID","Annotation_Level1")
clusterAnn$ClusterID<-c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
                        "21","22","23","24","25","26","27","28","29","30","31","32","33","34","35")
clusterAnn[clusterAnn$ClusterID %in% c(0),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(1),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(2),2]='Immunocytes'  
clusterAnn[clusterAnn$ClusterID %in% c(3),2]='Immunocytes'  
clusterAnn[clusterAnn$ClusterID %in% c(4),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(5),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(6),2]='Immunocytes'  
clusterAnn[clusterAnn$ClusterID %in% c(7),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(8),2]='NEs'
clusterAnn[clusterAnn$ClusterID %in% c(9),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(10),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(11),2]='NEs'
clusterAnn[clusterAnn$ClusterID %in% c(12),2]='NEs' 
clusterAnn[clusterAnn$ClusterID %in% c(13),2]='Immunocytes'
clusterAnn[clusterAnn$ClusterID %in% c(14),2]='Immunocytes'  
clusterAnn[clusterAnn$ClusterID %in% c(15),2]='Schwann_cells'  
clusterAnn[clusterAnn$ClusterID %in% c(16),2]='Immunocytes'
clusterAnn[clusterAnn$ClusterID %in% c(17),2]='NEs' 
clusterAnn[clusterAnn$ClusterID %in% c(18),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(19),2]='Fibroblasts'  
clusterAnn[clusterAnn$ClusterID %in% c(20),2]='Undefined'  
clusterAnn[clusterAnn$ClusterID %in% c(21),2]='NEs'
clusterAnn[clusterAnn$ClusterID %in% c(22),2]='Immunocytes' 
clusterAnn[clusterAnn$ClusterID %in% c(23),2]='Schwann_cells' 
clusterAnn[clusterAnn$ClusterID %in% c(24),2]='NEs' 
clusterAnn[clusterAnn$ClusterID %in% c(25),2]='Immunocytes'  
clusterAnn[clusterAnn$ClusterID %in% c(26),2]='Undefined'  
clusterAnn[clusterAnn$ClusterID %in% c(27),2]='NEs'  
clusterAnn[clusterAnn$ClusterID %in% c(28),2]='NEs'
clusterAnn[clusterAnn$ClusterID %in% c(29),2]='Immunocytes' 
clusterAnn[clusterAnn$ClusterID %in% c(30),2]='Immunocytes' 
clusterAnn[clusterAnn$ClusterID %in% c(31),2]='Undefined' 
clusterAnn[clusterAnn$ClusterID %in% c(32),2]='Immunocytes' 
clusterAnn[clusterAnn$ClusterID %in% c(33),2]='NEs' 
clusterAnn[clusterAnn$ClusterID %in% c(34),2]='Schwann_cells' 
clusterAnn[clusterAnn$ClusterID %in% c(35),2]='Fibroblasts' 

names(scRNA.sg@meta.data)

#一次性获得大量区分度颜色
library(RColorBrewer)
library(viridis)
# Generate colors
set1_colors <- brewer.pal(9, "Set1")
set2_colors <- brewer.pal(8, "Set2")
set3_colors <- brewer.pal(12, "Set3")
viridis_colors <- viridis(14) # 14 colors from the viridis palette
# Combine the color sets
my43color <- c(set1_colors, set2_colors, set3_colors, viridis_colors)
#用for循环对初级注释赋值
for(i in 1:nrow(clusterAnn)){
  scRNA.sg@meta.data[which(scRNA.sg@meta.data$RNA_snn_res.1== clusterAnn$ClusterID[i]),'Annotation_Level1'] <- clusterAnn$Annotation_Level1[i]}
table(scRNA.sg$Annotation_Level1)
#
library(SCpubr)
SCpubr::do_DimPlot(scRNA.sg,reduction = "umap",group.by = "Annotation_Level1",
                   plot.title = "Umap",
                   legend.position = "top",
                   legend.ncol = 6,
                   dims = c(2, 1),label = T#,colors.use = my13SCpubr_colors
)
DimPlot(scRNA.sg,group.by = "Sample")+DimPlot(scRNA.sg,group.by = "Annotation_Level1",label = T,raster=FALSE)
table(scRNA.sg$Sample_ID,scRNA.sg$Annotation_Level1)
table(scRNA.sg$Annotation_Level1)
### 3. 将不确定的细胞进行进一步细分，已获得准确的细胞注释结果###
# 使用 Seurat 的 subset 函数来子集化
scRNA.un<-scRNA.sg[,scRNA.sg$Annotation_Level1 %in% "Undefined"]
#scRNA.un <- subset(scRNA.sg, subset = Annotation_Level1 == "Undefined")

names(scRNA.un@meta.data)
scRNA.un@meta.data<-scRNA.un@meta.data[,-c(8:28)]
table(scRNA.un$Annotation_Level1)
#每次提取新的细胞亚群，都要走标准流程
scRNA.un <- NormalizeData(scRNA.un) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.un <- RunPCA(scRNA.un,npcs = 50,seed.use = 240701292)
scRNA.un <- RunHarmony(scRNA.un,group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA.un,ndims = 50)
pc.num=1:20
scRNA.un <- FindNeighbors(scRNA.un,dims = pc.num,reduction = "harmony")
scRNA.un <- FindClusters(scRNA.un, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.un)
sel.clust = "RNA_snn_res.1"
scRNA.un <- SetIdent(scRNA.un, value = sel.clust)
scRNA.un = RunUMAP(scRNA.un,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.un,group.by = "Sample")+DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T)
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "NRP1","CLEC4C","IRF7","IRF8",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "CMC1","GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
              )

SCpubr::do_DotPlot(scRNA.un,group.by = "RNA_snn_res.1",
                   assay = "RNA",
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
scRNA.un$Annotation_Correction<-ifelse(scRNA.un$RNA_snn_res.1 %in% c("0","3","4","5","6","7","10","11","12"),"Immunocytes","NEs")


table(scRNA.un$Annotation_Correction)
DimPlot(scRNA.un,group.by = "Annotation_Correction")

#纠正注释结果
meta.un<-as.data.frame(scRNA.un$Annotation_Correction)
colnames(meta.un)<-"Annotation"
meta.un$Barcode<-rownames(meta.un)
meta.unun<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Level1 %in% "Undefined")]$Annotation_Level1)
colnames(meta.unun)<-"Annotation"
meta.unun$Barcode<-rownames(meta.unun)
meta<-rbind(meta.unun,meta.un) 
meta<-meta[colnames(scRNA.sg),]
scRNA.sg$Annotation_Level1<-meta$Annotation
table(scRNA.sg$Annotation_Level1)
DimPlot(scRNA.sg,group.by = "Annotation_Level1",raster=FALSE)+FeaturePlot(scRNA.sg,c("PTPRC","PHOX2A"),split.by = "Annotation_Level1",raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC"),raster=FALSE)
table(scRNA.sg$Annotation_Level1)
rm(scRNA.un)



#然而发现NEs还是混合了施旺细胞，分别提出来进一步区分

### 4.进一步对混合细胞进行进一步区分 ###

scRNA.un<-scRNA.sg[,scRNA.sg$Annotation_Level1 %in% "Fibroblasts"]
table(scRNA.un$Annotation_Level1)
names(scRNA.un@meta.data)
scRNA.un@meta.data<-scRNA.un@meta.data[,-c(8:28)]
scRNA.un <- NormalizeData(scRNA.un) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.un <- RunPCA(scRNA.un,npcs = 50,seed.use = 240701292)
scRNA.un <- RunHarmony(scRNA.un,group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA.un,ndims = 50)
pc.num=1:20
scRNA.un <- FindNeighbors(scRNA.un,dims = pc.num,reduction = "harmony")
scRNA.un <- FindClusters(scRNA.un, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.un)
sel.clust = "RNA_snn_res.1"
scRNA.un <- SetIdent(scRNA.un, value = sel.clust)
scRNA.un = RunUMAP(scRNA.un,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.un,group.by = "Sample")+DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T)
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "NRP1","CLEC4C","IRF7","IRF8",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
              )
SCpubr::do_DotPlot(scRNA.un,group.by = "RNA_snn_res.1",
                   assay = "RNA",scale = T,
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")

scRNA.un$Annotation_Correction<-ifelse(scRNA.un$RNA_snn_res.1 %in% c("7"),"Schwann_cells",
                                       ifelse(scRNA.un$RNA_snn_res.1 %in% c("11"),"Immunocytes",
                                              ifelse(scRNA.un$RNA_snn_res.1 %in% c("9"),"NEs","Fibroblasts")))

table(scRNA.un$RNA_snn_res.1)
table(scRNA.un$Annotation_Correction)
DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T,raster=FALSE)+DimPlot(scRNA.un,group.by = "Annotation_Correction",raster=FALSE)
#纠正注释结果
meta.un<-as.data.frame(scRNA.un$Annotation_Correction)
colnames(meta.un)<-"Annotation"
meta.un$Barcode<-rownames(meta.un)
meta.unun<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Level1 %in% "Fibroblasts")]$Annotation_Level1)
colnames(meta.unun)<-"Annotation"
meta.unun$Barcode<-rownames(meta.unun)
meta<-rbind(meta.unun,meta.un) 
meta<-meta[colnames(scRNA.sg),]
scRNA.sg$Annotation_Level1<-meta$Annotation
table(scRNA.sg$Annotation_Level1)
DimPlot(scRNA.sg,group.by = "Annotation_Level1",raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC","PHOX2A"),split.by = "Annotation_Level1",raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC"),raster=FALSE)
table(scRNA.sg$Annotation_Level1)
rm(scRNA.un)




###将施旺细胞中的细胞分离
scRNA.un<-scRNA.sg[,scRNA.sg$Annotation_Level1 %in% "Schwann_cells"]
table(scRNA.un$Annotation_Level1)
names(scRNA.un@meta.data)
scRNA.un@meta.data<-scRNA.un@meta.data[,-c(8:28)]
scRNA.un <- NormalizeData(scRNA.un) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.un <- RunPCA(scRNA.un,npcs = 50,seed.use = 240701292)
scRNA.un <- RunHarmony(scRNA.un,group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA.un,ndims = 50)
pc.num=1:20
scRNA.un <- FindNeighbors(scRNA.un,dims = pc.num,reduction = "harmony")
scRNA.un <- FindClusters(scRNA.un, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.un)
sel.clust = "RNA_snn_res.1"
scRNA.un <- SetIdent(scRNA.un, value = sel.clust)
scRNA.un = RunUMAP(scRNA.un,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.un,group.by = "Sample")+DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T)
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "NRP1","CLEC4C","IRF7","IRF8",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
)
SCpubr::do_DotPlot(scRNA.un,group.by = "RNA_snn_res.1",
                   assay = "RNA",scale = T,
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")

scRNA.un$Annotation_Correction<-ifelse(scRNA.un$RNA_snn_res.1 %in% c("5","11"),"Endothelials",
                                       ifelse(scRNA.un$RNA_snn_res.1 %in% c("14"),"Immunocytes",
                                                     ifelse(scRNA.un$RNA_snn_res.1 %in% c("4","9","13"),"SCs","NEs")))

table(scRNA.un$RNA_snn_res.1)
table(scRNA.un$Annotation_Correction)
DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T,raster=FALSE)+DimPlot(scRNA.un,group.by = "Annotation_Correction",raster=FALSE)
#纠正注释结果
meta.un<-as.data.frame(scRNA.un$Annotation_Correction)
colnames(meta.un)<-"Annotation"
meta.un$Barcode<-rownames(meta.un)
meta.unun<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Level1 %in% "Schwann_cells")]$Annotation_Level1)
colnames(meta.unun)<-"Annotation"
meta.unun$Barcode<-rownames(meta.unun)
meta<-rbind(meta.unun,meta.un) 
meta<-meta[colnames(scRNA.sg),]
scRNA.sg$Annotation_Level1<-meta$Annotation
table(scRNA.sg$Annotation_Level1)
DimPlot(scRNA.sg,group.by = "Annotation_Level1",label = T,raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC","PHOX2A"),split.by = "Annotation_Level1",raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC"),raster=FALSE)
table(scRNA.sg$Annotation_Level1)
rm(scRNA.un)


###将混杂的内皮细胞分离
scRNA.un<-scRNA.sg[,scRNA.sg$Annotation_Level1 %in% "Endothelials"]
table(scRNA.un$Annotation_Level1)
names(scRNA.un@meta.data)
scRNA.un@meta.data<-scRNA.un@meta.data[,-c(8:28)]
scRNA.un <- NormalizeData(scRNA.un) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.un <- RunPCA(scRNA.un,npcs = 50,seed.use = 240701292)
scRNA.un <- RunHarmony(scRNA.un,group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA.un,ndims = 50)
pc.num=1:20
scRNA.un <- FindNeighbors(scRNA.un,dims = pc.num,reduction = "harmony")
scRNA.un <- FindClusters(scRNA.un, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.un)
sel.clust = "RNA_snn_res.1"
scRNA.un <- SetIdent(scRNA.un, value = sel.clust)
scRNA.un = RunUMAP(scRNA.un,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.un,group.by = "Sample")+DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T)
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "NRP1","CLEC4C","IRF7","IRF8",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
              )
SCpubr::do_DotPlot(scRNA.un,group.by = "RNA_snn_res.1",
                   assay = "RNA",scale = T,
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")

scRNA.un$Annotation_Correction<-ifelse(scRNA.un$RNA_snn_res.1 %in% c("7"),"NEs","Endothelials")

table(scRNA.un$RNA_snn_res.1)
table(scRNA.un$Annotation_Correction)
DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T,raster=FALSE)+DimPlot(scRNA.un,group.by = "Annotation_Correction",raster=FALSE)
#纠正注释结果
meta.un<-as.data.frame(scRNA.un$Annotation_Correction)
colnames(meta.un)<-"Annotation"
meta.un$Barcode<-rownames(meta.un)
meta.unun<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Level1 %in% "Endothelials")]$Annotation_Level1)
colnames(meta.unun)<-"Annotation"
meta.unun$Barcode<-rownames(meta.unun)
meta<-rbind(meta.unun,meta.un) 
meta<-meta[colnames(scRNA.sg),]
scRNA.sg$Annotation_Level1<-meta$Annotation
table(scRNA.sg$Annotation_Level1)
DimPlot(scRNA.sg,group.by = "Annotation_Level1",label = T,raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC","PHOX2A"),split.by = "Annotation_Level1",raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC"),raster=FALSE)
table(scRNA.sg$Annotation_Level1)
rm(scRNA.un)



### 4.进一步对免疫细胞进行精细注释 ###
#提取免疫细胞，走一套标准流程
scRNA.imm<-scRNA.sg[,scRNA.sg$Annotation_Level1 %in% "Immunocytes"]
names(scRNA.imm@meta.data)
scRNA.imm@meta.data<-scRNA.imm@meta.data[,-c(8:28)]
scRNA.imm <- NormalizeData(scRNA.imm) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.imm <- RunPCA(scRNA.imm,npcs = 50,seed.use = 240701292)
scRNA.imm <- RunHarmony(scRNA.imm,group.by.vars = "orig.ident",max.iter.harmony = 10)
ElbowPlot(scRNA.imm,ndims = 50)
pc.num=1:20
scRNA.imm <- FindNeighbors(scRNA.imm,dims = pc.num,reduction = "harmony")
scRNA.imm <- FindClusters(scRNA.imm, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.imm)
scRNA.imm = RunUMAP(scRNA.imm,reduction = "harmony",dims = pc.num,seed.use = 240701292)
DimPlot(scRNA.imm,group.by = "Sample",raster=FALSE)
sel.clust = "RNA_snn_res.1"
scRNA.imm <- SetIdent(scRNA.imm, value = sel.clust)
table(scRNA.imm@active.ident)
#收集到的marker
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "IL3RA","NRP1",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
)
SCpubr::do_DotPlot(scRNA.imm,group.by = "RNA_snn_res.1",
                   assay = "RNA",scale = T,
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
DimPlot(scRNA.imm,group.by = "RNA_snn_res.1",label = T,raster=FALSE)+FeaturePlot(scRNA.imm,c("CD14","C1QA","FCGR3A","CD3E"),raster=FALSE)
table(scRNA.imm$RNA_snn_res.1)
#创建一个clusterAnn_imm空对象
clusterAnn_imm = as.data.frame(matrix(nrow=25,ncol=2),col.names=NULL ,row.names=NULL ) 
colnames(clusterAnn_imm)<-c("ClusterID","Annotation_Level2")
clusterAnn_imm$ClusterID<-c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
                            "19","20","21","22","23","24")
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(0),2]='T_cells'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(1),2]='NEs'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(2),2]='NEs'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(3),2]='Myeloid_cells'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(4),2]='T_cells'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(5),2]='Undefined'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(6),2]='NEs'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(7),2]='B_cells'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(8),2]='T_cells'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(9),2]='T_cells'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(10),2]='Undefined'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(11),2]='Myeloid_cells'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(12),2]='Myeloid_cells'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(13),2]='T_cells' 
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(14),2]='Myeloid_cells'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(15),2]='Undefined'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(16),2]='NEs'  
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(17),2]='Undefined'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(18),2]='T_cells'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(19),2]='Undefined'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(20),2]='T_cells'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(21),2]='Undefined'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(22),2]='Undefined'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(23),2]='NEs'
clusterAnn_imm[clusterAnn_imm$ClusterID %in% c(24),2]='SCs'

#赋值
for(i in 1:nrow(clusterAnn_imm)){
  scRNA.imm@meta.data[which(scRNA.imm@meta.data$RNA_snn_res.1== clusterAnn_imm$ClusterID[i]),'Annotation_Level2'] <- clusterAnn_imm$Annotation_Level2[i]}
table(scRNA.imm$Annotation_Level2)
#检查一下
DimPlot(scRNA.imm,group.by = "RNA_snn_res.1",label = T,raster = FALSE)+
  DimPlot(scRNA.imm,group.by = "Annotation_Level2",label = T,raster = FALSE)+FeaturePlot(scRNA.imm,c("CD14","C1QA","FCGR3A","CD3E","S100A8","NKG7"))

###将Undifined提取出来细分
scRNA.un<-scRNA.imm[,scRNA.imm$Annotation_Level2 %in% "Undefined"]
table(scRNA.un$Annotation_Level2)
names(scRNA.un@meta.data)

scRNA.un <- NormalizeData(scRNA.un) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.un <- RunPCA(scRNA.un,npcs = 50,seed.use = 240701292)
scRNA.un <- RunHarmony(scRNA.un,group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA.un,ndims = 50)
pc.num=1:20
scRNA.un <- FindNeighbors(scRNA.un,dims = pc.num,reduction = "harmony")
scRNA.un <- FindClusters(scRNA.un, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.un)
sel.clust = "RNA_snn_res.1"
scRNA.un <- SetIdent(scRNA.un, value = sel.clust)
scRNA.un = RunUMAP(scRNA.un,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.un,group.by = "Sample")+DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T)
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "IL3RA","NRP1",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
              )
SCpubr::do_DotPlot(scRNA.un,group.by = "RNA_snn_res.1",
                   assay = "RNA",scale = T,
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")

scRNA.un$Annotation_Correction<-ifelse(scRNA.un$RNA_snn_res.1 %in% c("0","5","10","11"),"Myeloid_cells",
                                       ifelse(scRNA.un$RNA_snn_res.1 %in% c("1","9"),"NK_cells",
                                              ifelse(scRNA.un$RNA_snn_res.1 %in% c("2","7","8","12"),"T_cells",
                                                     ifelse(scRNA.un$RNA_snn_res.1 %in% c("3","14"),"NEs",
                                                            ifelse(scRNA.un$RNA_snn_res.1 %in% c("4"),"DCs","B_cells")))))

table(scRNA.un$RNA_snn_res.1)
table(scRNA.un$Annotation_Correction)
DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T,raster=FALSE)+DimPlot(scRNA.un,group.by = "Annotation_Correction",raster=FALSE)
#纠正注释结果
meta.un<-as.data.frame(scRNA.un$Annotation_Correction)
colnames(meta.un)<-"Annotation"
meta.un$Barcode<-rownames(meta.un)
meta.unun<-as.data.frame(scRNA.imm[,!(scRNA.imm$Annotation_Level2 %in% "Undefined")]$Annotation_Level2)
colnames(meta.unun)<-"Annotation"
meta.unun$Barcode<-rownames(meta.unun)
meta<-rbind(meta.unun,meta.un) 
meta<-meta[colnames(scRNA.imm),]
scRNA.imm$Annotation_Level2<-meta$Annotation
table(scRNA.imm$Annotation_Level2)
DimPlot(scRNA.imm,group.by = "Annotation_Level2",label = T,raster=FALSE)
FeaturePlot(scRNscRNA.immA.imm,c("PTPRC","PHOX2A"),split.by = "Annotation_Level2",raster=FALSE)
FeaturePlot(scRNA.imm,c("PTPRC"),raster=FALSE)
table(scRNA.imm$Annotation_Level2)
rm(scRNA.un)

####将信息返回到原数据
meta.un<-as.data.frame(scRNA.imm$Annotation_Level2)
colnames(meta.un)<-"Annotation"
meta.un$Barcode<-rownames(meta.un)
table(scRNA.sg$Annotation_Level1)
meta.unun<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Level1 %in% "Immunocytes")]$Annotation_Level1)
colnames(meta.unun)<-"Annotation"
meta.unun$Barcode<-rownames(meta.unun)
meta<-rbind(meta.unun,meta.un) 
meta<-meta[colnames(scRNA.sg),]
scRNA.sg$Annotation_Level2<-meta$Annotation
table(scRNA.sg$Annotation_Level2)
table(scRNA.sg$Annotation_Level1,scRNA.sg$Annotation_Level2)
table(scRNA.sg$Sample_ID,scRNA.sg$Annotation_Level2)
DimPlot(scRNA.sg,group.by = "Annotation_Level2",label = T,raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC","PHOX2A"),split.by = "Annotation_Level2",raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC"),raster=FALSE)
table(scRNA.sg$Annotation_Level1)



### 对髓系细胞进行注释 ###本项目研究施旺细胞暂时先不注释髓系细胞
#替换掉之前的scRNA.un，节省内存，它已经没用了
scRNA.un<-scRNA.imm[,scRNA.imm$Annotation_Level2 %in% "Myeloid_cells"]
table(scRNA.un$Annotation_Level2)
names(scRNA.un@meta.data)

scRNA.un <- NormalizeData(scRNA.un) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.un <- RunPCA(scRNA.un,npcs = 50,seed.use = 240701292)
scRNA.un <- RunHarmony(scRNA.un,group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA.un,ndims = 50)
pc.num=1:20
scRNA.un <- FindNeighbors(scRNA.un,dims = pc.num,reduction = "harmony")
scRNA.un <- FindClusters(scRNA.un, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.un)
sel.clust = "RNA_snn_res.1"
scRNA.un <- SetIdent(scRNA.un, value = sel.clust)
scRNA.un = RunUMAP(scRNA.un,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.un,group.by = "Sample")+DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T)
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","TYROBP",#Myeloid_cells
              "CD68","CD163","MRC1","C1QA","C1QB","C1QC","MSR1",#巨噬细胞 区分单核/巨噬/DC系参考PMID33420488
              "CD14","FCGR3A",#单核细胞高表达
              "IL3RA","CLEC4C","IRF7","IRF8",#DC
              "S100A8","S100A9","MMP9","CSF3R","FCGR3B",#"Neutrophil"
              "KIT","IL1RL1",#肥大细胞
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
)
SCpubr::do_DotPlot(scRNA.un,group.by = "RNA_snn_res.1",
                   assay = "RNA",scale = T,
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")


###发现施万细胞、成纤维细胞、内皮细胞有些混杂 分一下
scRNA.un<-scRNA.sg[,scRNA.sg$Annotation_Level2 %in% c("SCs","Endothelials","Fibroblasts")]
table(scRNA.un$Annotation_Level2)
names(scRNA.un@meta.data)
scRNA.un@meta.data<-scRNA.un@meta.data[,-c(8:28)]
scRNA.un <- NormalizeData(scRNA.un) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.un <- RunPCA(scRNA.un,npcs = 50,seed.use = 240701292)
scRNA.un <- RunHarmony(scRNA.un,group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA.un,ndims = 50)
pc.num=1:20
scRNA.un <- FindNeighbors(scRNA.un,dims = pc.num,reduction = "harmony")
scRNA.un <- FindClusters(scRNA.un, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.un)
sel.clust = "RNA_snn_res.1"
scRNA.un <- SetIdent(scRNA.un, value = sel.clust)
scRNA.un = RunUMAP(scRNA.un,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.un,group.by = "Sample")+DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T)
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "NRP1","CLEC4C","IRF7","IRF8",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
)
SCpubr::do_DotPlot(scRNA.un,group.by = "RNA_snn_res.1",
                   assay = "RNA",scale = T,
                   colors.use = c("#83BF46", "#DF7CB1"),
                   features = Markergene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")

scRNA.un$Annotation_Correction<-ifelse(scRNA.un$RNA_snn_res.1 %in% c("0","1","2","3","6","8","10","15","16","17","18"),"Fibroblasts",
                                       ifelse(scRNA.un$RNA_snn_res.1 %in% c("4","7","9","11"),"SCs",
                                              ifelse(scRNA.un$RNA_snn_res.1 %in% c("13"),"T_cells","Endothelials")))

table(scRNA.un$RNA_snn_res.1)
table(scRNA.un$Annotation_Correction)
DimPlot(scRNA.un,group.by = "RNA_snn_res.1",label = T,raster=FALSE)+DimPlot(scRNA.un,group.by = "Annotation_Correction",raster=FALSE)
#纠正注释结果
meta.un<-as.data.frame(scRNA.un$Annotation_Correction)
colnames(meta.un)<-"Annotation"
meta.un$Barcode<-rownames(meta.un)
meta.unun<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Level2 %in% c("SCs","Endothelials","Fibroblasts"))]$Annotation_Level2)
colnames(meta.unun)<-"Annotation"
meta.unun$Barcode<-rownames(meta.unun)
meta<-rbind(meta.unun,meta.un) 
meta<-meta[colnames(scRNA.sg),]
scRNA.sg$Annotation_Level2<-meta$Annotation
table(scRNA.sg$Annotation_Level2)
DimPlot(scRNA.sg,group.by = "Annotation_Level2",label = T,raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC","PHOX2A"),split.by = "Annotation_Level2",raster=FALSE)
FeaturePlot(scRNA.sg,c("PTPRC"),raster=FALSE)
table(scRNA.sg$Annotation_Level1)
rm(scRNA.un)






qsave(scRNA.sg,file = "~/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/scRNA_Annotationed.qs")
scRNA.sg<-qread("~/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/scRNA_Annotationed.qs")
rm(list = ls())
gc()


#################----四、InferCNV推断上皮细胞中的恶性细胞----###################
### 1. 提出肿瘤细胞和参考细胞  ###
#（免疫细胞不适合作为参考对照细胞，且数量太多耗内存）
scRNA.cnv<-qread("~/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/scRNA_Annotationed_AS.qs")

DimPlot(scRNA.cnv,group.by = "Annotation_Level2",label = T,raster=FALSE)
#scRNA.cnv<-scRNA.sg[,!(scRNA.sg$Annotation_Level1 %in% "Immunocytes")]
#走一遍标准流程，把上皮细胞分为多个cluster，以区分其中潜在的恶性细胞
scRNA.cnv <- NormalizeData(scRNA.cnv) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.cnv <- RunPCA(scRNA.cnv,npcs = 50,seed.use = 240701292)
scRNA.cnv <- RunHarmony(scRNA.cnv,group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA.cnv,ndims = 50)
pc.num=1:30
scRNA.cnv <- FindNeighbors(scRNA.cnv,dims = pc.num,reduction = "harmony")
scRNA.cnv <- FindClusters(scRNA.cnv, resolution = 1,random.seed = 240701292)
library(clustree)

sel.clust = "RNA_snn_res.1"
scRNA.cnv <- SetIdent(scRNA.cnv, value = sel.clust)
scRNA.cnv = RunUMAP(scRNA.cnv,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.cnv,group.by = "Annotation_Level6",raster = FALSE)+DimPlot(scRNA.cnv,group.by = "RNA_snn_res.1",raster = FALSE)
scRNA.cnv$Annotation_cnv<-ifelse(scRNA.cnv$Annotation_Level6 %in% "Neuroendocrine_cells",paste0("Neuroendocrine_C",scRNA.cnv$RNA_snn_res.1),scRNA.cnv$Annotation_Level6)
table(scRNA.cnv$Annotation_cnv)

### 2. InferCNV所需要的输入文件准备  ###
#准备count矩阵、细胞注释文件、基因位置|
library(infercnv)
### 2.1 制作GTF文件,
counts<-GetAssayData(object = scRNA.cnv,assay = "RNA",slot = "counts")
library(AnnoProbe)
geneInfor=annoGene(rownames(counts),"SYMBOL")
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
table(geneInfor$chr)
### 2.2 制作表达矩阵
counts<-as.data.frame(scRNA.cnv@assays$RNA@counts)
counts[1:4,1:4]
colnames(counts)<-gsub('-',replacement = '_',colnames(counts))
counts=counts[rownames(counts) %in% geneInfor[,1],]
counts=counts[match( geneInfor[,1], rownames(counts) ),] 
dim(counts)
expFile='~/scRNA/R_Project/Schwann_project/Integrate/GSE137804_GSE147766/InferCNV/expFile_counts.txt'
write.table(counts,file = expFile,sep = '\t',quote = F)

### 2.3 制作基因order表格
head(geneInfor)
geneFile='~/scRNA/R_Project/Schwann_project/Integrate/GSE137804_GSE147766/InferCNV/geneFile.txt'
#把单字的染色体编号改成统一制式
geneInfor$chr[which(geneInfor$chr=='chr1')]='chr01'
geneInfor$chr[which(geneInfor$chr=='chr2')]='chr02'
geneInfor$chr[which(geneInfor$chr=='chr3')]='chr03'
geneInfor$chr[which(geneInfor$chr=='chr4')]='chr04'
geneInfor$chr[which(geneInfor$chr=='chr5')]='chr05'
geneInfor$chr[which(geneInfor$chr=='chr6')]='chr06'
geneInfor$chr[which(geneInfor$chr=='chr7')]='chr07'
geneInfor$chr[which(geneInfor$chr=='chr8')]='chr08'
geneInfor$chr[which(geneInfor$chr=='chr9')]='chr09'
geneInfor<-geneInfor[with(geneInfor,order(chr)),]
head(geneInfor)
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)

### 2.4 制作注释文件
table(scRNA.cnv@meta.data$Annotation_Level6)
table(scRNA.cnv@meta.data$Annotation_cnv)
groupinfo<-scRNA.cnv@meta.data
groupinfo<-data.frame(rownames(groupinfo),groupinfo$Annotation_Level6,groupinfo$Annotation_cnv)
table(groupinfo$groupinfo.Annotation_Level6)
table(groupinfo$groupinfo.Annotation_cnv)
sum(table(groupinfo$groupinfo.Annotation_Level6))
sum(table(groupinfo$groupinfo.Annotation_cnv))
groupFiles='~/scRNA/R_Project/Schwann_project/Integrate/GSE137804_GSE147766/InferCNV/groupFiles.txt'
head(groupinfo)
groupinfo$rownames.groupinfo.<-gsub('-',replacement = '_',groupinfo$rownames.groupinfo.)
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
options(future.globals.maxSize=9891289600)

### 2.5 开始构建CNV对象 ###
#注意：一定要按照上面步骤构建好文件，不可以用变量，而且不按步骤来一堆问题
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,#表达矩阵
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names= c("Fibroblasts","Endothelial_Cells",
                                                       "Schwann_cells"))

## 2.6 切换到InferCNV_Run.R跑出结果
saveRDS(infercnv_obj,file = '~/scRNA/R_Project/Schwann_project/Integrate/GSE137804_GSE147766/InferCNV/infercnv_ref.RDS')
##开始跑，这一步先试运行一下，没问题就在linux上用Rscript xxx.R &> nohup&挂后台跑，因为耗时长
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir='~/R_project/PDAC/Article/tmp/InferCNV/PDAC/result/', 
                              min_cells_per_gene = 5,
                              cluster_by_groups=T, 
                              denoise=TRUE,
                              num_threads =30, 
                              HMM=TRUE,
                              no_prelim_plot = T,
                              output_format = "pdf")
#PS：使用Rscript跑一定要在R脚本最后加一句save(xxxx)，不然白跑
saveRDS(scRNA.sg,"~/R_project/PDAC/Article/Intermediate/scRNA.sg_231029.RDS")
#rm(scRNA.un,scRNA.imm,scRNA.list.scRNAlist,scRNA)
#
library(Seurat)
#查看PIDA4和GSDMD基因表达
VlnPlot(scRNA.sg,features = "SLC39A4",group.by = "Final_Type2")+
  VlnPlot(scRNA.sg,features = "GSDMD",group.by = "Final_Type2")

### 3. 定量分析细胞的CNV评分来推断细胞是否为恶性  ###
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
infercnv_final_obj = readRDS("~/R_project/PDAC/Article/tmp/InferCNV/PDAC/results/run.final.infercnv_obj")
expr <- infercnv_final_obj@expr.data
expr[1:10,1:4]
infercnv_final_obj@reference_grouped_cell_indices
#直接计算每个细胞的CNVscore
expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
#细胞barcode赋值 
CNV_score$Barcode=rownames(CNV_score)
colnames(CNV_score)[1]<-"CNV_Score"
#把cluster和CNV得分合并
pri_cluster<-as.data.frame(scRNA.cnv$Annotation_Level3)
pri_cluster$Barcode<-colnames(scRNA.cnv)
colnames(pri_cluster)[1]<-"Sce_Cluster"
pri_cluster$Barcode<-gsub("-","_",pri_cluster$Barcode)
#
CNV_Cluster<-left_join(CNV_score,pri_cluster)

#画图
color_v<-my43color
CNV_Cluster%>%ggplot(aes(Sce_Cluster,CNV_Score))+geom_violin(aes(fill=Sce_Cluster),color="NA")+
  scale_fill_manual(values = color_v)+
  theme_bw()+scale_y_continuous(limits = c(0, 0.005), breaks = c(0, 0.0025, 0.005))+#控制Y轴范围
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5
  ))   #让横轴倾斜
rownames(CNV_Cluster)<-CNV_Cluster$Barcode
CNV_Cluster<-CNV_Cluster[gsub("-","_",colnames(scRNA.cnv)),]
scRNA.cnv$CNV_Score<-CNV_Cluster$CNV_Score
table(gsub("-","_",colnames(scRNA.cnv)) == rownames(CNV_Cluster))

### 4. 将潜在的上皮和恶性细胞混合群体分开 ###
scRNA.un<-scRNA.cnv[,scRNA.cnv$Annotation_Level3 %in% c("Acinar_Cells","Epithelial_C12","Epithelial_C16","Epithelial_C17")]
#每次提取新的细胞亚群，都要走标准流程
scRNA.un <- NormalizeData(scRNA.un) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.un <- RunPCA(scRNA.un,npcs = 50,seed.use = 240701292)
scRNA.un <- RunHarmony(scRNA.un,group.by.vars = "orig.ident",max.iter.harmony = 10)
ElbowPlot(scRNA.un,ndims = 50)
pc.num=1:20
scRNA.un <- FindNeighbors(scRNA.un,dims = pc.num,reduction = "harmony")
scRNA.un <- FindClusters(scRNA.un, resolution = seq(5,15,by=1),random.seed = 240701292)
library(clustree)
clustree(scRNA.un)
sel.clust = "RNA_snn_res.10"
scRNA.un <- SetIdent(scRNA.un, value = sel.clust)
scRNA.un = RunUMAP(scRNA.un,reduction = "harmony",dims = pc.num)
#DimPlot(scRNA.un,group.by = "RNA_snn_res.10",label = T)+FeaturePlot(scRNA.un,features = "CNV_Score")
meta <- as.data.frame(scRNA.un@meta.data)
meta%>%ggplot(aes(RNA_snn_res.10,CNV_Score))+geom_boxplot(aes(fill=RNA_snn_res.10),color="NA")+
  #scale_fill_manual(values = my43color)+
  theme_bw()+scale_y_continuous(limits = c(0, 0.0025), breaks = c(0, 0.00125, 0.0025))+#控制Y轴范围
  theme(axis.text.x = element_text(angle = 45#, hjust = 0.5, vjust = 0.5
  )) 
DimPlot(scRNA.un,group.by = "Malignant")
scRNA.un$Malignant<-ifelse(scRNA.un$RNA_snn_res.10 %in% c("1","16","18","19","53","7","8"),"Malignant_Cells",paste0(str_split_fixed(scRNA.un$Annotation_Level3,"_",2)[,1],"_Cells"))
table(scRNA.un$Malignant,scRNA.un$Annotation_Level3)
#可视化比较分类效果
meta <- as.data.frame(scRNA.un@meta.data)
meta%>%ggplot(aes(Malignant,CNV_Score))+geom_violin(aes(fill=Malignant),color="NA")+
  scale_fill_manual(values = my43color)+
  theme_bw()+scale_y_continuous(limits = c(0, 0.005), breaks = c(0, 0.0025, 0.005))+#控制Y轴范围
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5
  ))   #让横轴倾斜

#提出混合细胞的meta信息
meta.mix <- as.data.frame(scRNA.un$Malignant)#这一步千万不能加@meta.data，否则会失去行名
colnames(meta.mix)<-"Annotation_Level4"
meta.mix$Barcode<-rownames(meta.mix)
#提出非混合细胞的meta信息
meta.nonmix<-as.data.frame(scRNA.cnv[,!(scRNA.cnv$Annotation_Level3 %in% c("Acinar_Cells","Epithelial_C12","Epithelial_C16","Epithelial_C17"))]$Annotation_Level3)
colnames(meta.nonmix)<-"Annotation_Level4"
meta.nonmix$Barcode<-rownames(meta.nonmix)
#合并并组合到scRNA.cnv
meta.all<-rbind(meta.nonmix,meta.mix)
meta.all<-meta.all[colnames(scRNA.cnv),]
scRNA.cnv$Annotation_Level4<-meta.all$Annotation_Level4
table(scRNA.cnv$Annotation_Level3,scRNA.cnv$Annotation_Level4)
table(scRNA.cnv$Annotation_Level4)
#检验
meta<-scRNA.cnv@meta.data
meta%>%ggplot(aes(Annotation_Level4,CNV_Score))+geom_violin(aes(fill=Annotation_Level4),color="NA")+
  scale_fill_manual(values = my43color)+
  theme_bw()+scale_y_continuous(limits = c(0, 0.005), breaks = c(0, 0.0025, 0.005))+#控制Y轴范围
  theme(axis.text.x = element_text(angle = 45#, hjust = 0.5, vjust = 0.5
  ))  
#先区分出恶性细胞
scRNA.cnv$Malignant<-ifelse(scRNA.cnv$Annotation_Level4 %in% c("Acinar_Cells","Endothelial_Cells","Epithelial_C18",
                                                               "Epithelial_C9","Epithelial_C8","Epithelial_Cells",
                                                               "Fibroblasts","Periocytes"),"Non_Malignant","Malignant")
#再把非恶性细胞的类型还原
scRNA.cnv$Annotation_Level5<-ifelse(scRNA.cnv$Malignant %in% "Malignant","Malignant_Cells",
                                    ifelse(str_split_fixed(scRNA.cnv$Annotation_Level4,"_",2)[,1] %in% "Epithelial","Epithelial_Cells",scRNA.cnv$Annotation_Level4))
table(scRNA.cnv$Annotation_Level5)
meta<-scRNA.cnv@meta.data
meta%>%ggplot(aes(Annotation_Level5,CNV_Score))+geom_violin(aes(fill=Annotation_Level5),color="NA")+
  scale_fill_manual(values = my43color)+
  theme_bw()+scale_y_continuous(limits = c(0, 0.005), breaks = c(0, 0.0025, 0.005))+#控制Y轴范围
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5
  ))  

### 5. 确认处理得当，开始把所有注释信息归拢到scRNA.sg ###
#提出CNV细胞的meta信息
meta.cnv <- data.frame(row.names = colnames(scRNA.cnv),
                       "Annotation" = scRNA.cnv$Annotation_Level5,
                       "CNV_Score" = scRNA.cnv$CNV_Score)#这一步千万不能加@meta.data，否则会失去行名

#提出免疫细胞的meta信息
table(scRNA.sg$Annotation_Level1,scRNA.sg$Annotation_Level2)
meta.imm<-data.frame(row.names = colnames(scRNA.sg[,scRNA.sg$Annotation_Level1 %in% "Immunocytes"]),
                     "Annotation" =scRNA.sg[,scRNA.sg$Annotation_Level1 %in% "Immunocytes"]$Annotation_Level2,
                     "CNV_Score" =rep(NA,ncol(scRNA.sg[,scRNA.sg$Annotation_Level1 %in% "Immunocytes"])))

#合并并组合到scRNA.sg
meta.all<-rbind(meta.cnv,meta.imm)
meta.all<-meta.all[colnames(scRNA.sg),]
scRNA.sg$Annotation_Level3<-meta.all$Annotation
scRNA.sg$CNV_Score<-meta.all$CNV_Score
table(scRNA.sg$Annotation_Level2,scRNA.sg$Annotation_Level3)
table(scRNA.sg$Annotation_Level4)
#检验
scRNA.epi<-scRNA.sg[,scRNA.sg$Annotation_Level3 %in% c("Acinar_Cells","Epithelial_Cells","Malignant_Cells")]
#走一遍标准流程，把上皮细胞分为多个cluster，以区分其中潜在的恶性细胞
scRNA.epi <- NormalizeData(scRNA.epi) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.epi <- RunPCA(scRNA.epi,npcs = 50,seed.use = 240701292) %>% RunHarmony(group.by.vars = "orig.ident",max.iter.harmony = 10)
ElbowPlot(scRNA.epi,ndims = 50)
pc.num=1:30
scRNA.epi <- FindNeighbors(scRNA.epi,dims = pc.num,reduction = "harmony")
scRNA.epi <- FindClusters(scRNA.epi, resolution = seq(0,1,by=0.1),random.seed = 240701292)
scRNA.epi = RunUMAP(scRNA.epi,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.epi,group.by = "Annotation_Level3",pt.size = 1,cols = my43color[c(2,3,1)])+FeaturePlot(scRNA.epi,"CNV_Score",
                                                                                                     max.cutoff = 0.005,pt.size = 1,cols = c("lightgrey","darkred"))



####将infercnvpy结果保存
metadata <- read.csv("/home/chenzb/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/meta_CNV.csv", row.names = 1)
View(metadata)
scRNA.sg<-qread("~/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/scRNA_Annotationed_AS.qs")
metadata<-metadata[colnames(scRNA.sg),]
scRNA.sg@meta.data<-metadata
View(metadata)
table(scRNA.sg$Annotation_Final,scRNA.sg$cnv_status)
table(scRNA.sg$Annotation_Level2,scRNA.sg$cnv_status)
scRNA.sg <- NormalizeData(scRNA.sg) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.sg <- RunPCA(scRNA.sg,seed.use = 240701292)
ElbowPlot(scRNA.sg, ndims = 50)
###PCA
#scRNA.sg <- RunPCA(scRNA.sg,npcs = 50,seed.use = 240701292,assay = "SCT")
## 整合方法：1. 单个样本间进行整合（推荐，效果更好）
scRNA.sg <- RunHarmony(scRNA.sg,group.by.vars = "Sample",#assay.use = "SCT",
                       max.iter.harmony = 10)
# lambda参数默认是1，唯一一个决定了Harmony的整合力度的参数，值越小，整合力度越大，取0.5-2
pc.num=1:30
###找cluster
scRNA.sg <- FindNeighbors(scRNA.sg,dims = pc.num,reduction = "harmony"#,assay = "SCT"
)#本次不用SCT
scRNA.sg <- FindClusters(scRNA.sg,resolution = seq(0,2,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.sg)
DimPlot(scRNA.sg,group.by = "RNA_snn_res.2",label = T,raster=FALSE)
DimPlot(scRNA.sg,group.by = "Annotation_Final",label = T,raster=FALSE)
DimPlot(scRNA.sg,group.by = "Annotation_Final",raster=FALSE)
DimPlot(scRNA.sg,group.by = "Annotation_level2",raster=FALSE)
View(scRNA.sg)
DimPlot(scRNA.sg,group.by = "Annotation_Level2",raster=FALSE)
sel.clust = "RNA_snn_res.1"
scRNA.sg <- SetIdent(scRNA.sg, value = sel.clust)
table(scRNA.sg@active.ident)
scRNA.sg = RunUMAP(scRNA.sg,reduction = "harmony",dims = pc.num,random.seed = 240701292)
DimPlot(scRNA.sg,group.by = "Annotation_Level2",raster=FALSE)
DimPlot(scRNA.sg,group.by = "RNA_snn_res.1",label = T,raster=FALSE)
DimPlot(scRNA.sg,group.by = "Annotation_Level2",raster=FALSE)
qsave(scRNA.sg,file = "~/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/scRNA.sg_Final_annotation.qs")
scRNA.sg<-qread("~/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/scRNA.sg_Final_annotation.qs")
table(scRNA.sg$Annotation_Final)
rm(list = ls())
gc()



#################----五、数据处理结果的基本可视化展示----###################
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(patchwork)
library(dplyr)
scRNA.sg<-qread("~/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/scRNA.sg_Final_annotation.qs")
table(scRNA.sg$Annotation_Final)
DimPlot(scRNA.sg, group.by = "Annotation_Final",cols = my43color[2:15],
        label = TRUE,shuffle = T,repel = T,seed = 240701292, label.size = 3) +
  theme_minimal() +
  theme(legend.position = "right") +
  #NoLegend() +
  ggtitle("UMAP Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
#一次性获得大量区分度颜色
library(RColorBrewer)
library(viridis)
# Generate colors
set1_colors <- brewer.pal(9, "Set1")
set2_colors <- brewer.pal(8, "Set2")
set3_colors <- brewer.pal(12, "Set3")
viridis_colors <- viridis(14) # 14 colors from the viridis palette
# Combine the color sets
my43color <- c(set1_colors, set2_colors, set3_colors, viridis_colors)
DimPlot(scRNA.sg, group.by = "Annotation_Final",cols = my43color[2:15],
        label = TRUE,shuffle = T,repel = T,seed = 240701292, label.size = 3) +
  theme_minimal() +
  theme(legend.position = "right") +
  #NoLegend() +
  ggtitle("UMAP Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
DimPlot(scRNA.sg, group.by = "Annotation_Final",cols = my43color[2:15],label = TRUE,
        shuffle = T,repel = T,seed = 240701292, label.size = 3,raster = FALSE) +
  theme_minimal() +
  theme(legend.position = "right") +
  #NoLegend() +
  ggtitle("UMAP Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
DimPlot(scRNA.sg, group.by = "Annotation_Level2",cols = my43color[2:15],label = TRUE,
        shuffle = T,repel = T,seed = 240701292, label.size = 3,raster = FALSE) +
  theme_minimal() +
  theme(legend.position = "right") +
  #NoLegend() +
  ggtitle("UMAP Clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
### 2. 画marker图 ###
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "NRP1","CLEC4C","IRF7","IRF8",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "MZB1","IGHG3","JCHAIN",#浆细胞
              "CMC1","GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
)
# Calculate average expression and percent detection
scRNA.sg2 <- AverageExpression(scRNA.sg, assay = "RNA", slot = "data", features = Markergene)
library(ggplot2)
table(scRNA.sg$Annotation_Final)
### 2. 画marker图 ###
Markergene<-c("PTPRC",#免疫细胞
              "LYZ","C1QA","S100A9","CD1C","CD14","TYROBP",#Myeloid_cells
              "NRP1","CLEC4C","IRF7","IRF8",#DC
              "CD3D","CD3E","CD52","CD7",#T细胞
              "BANK1","MS4A1","CD79A","CD79B",#B细胞
              "CMC1","GNLY","GZMB","NKG7",#NK细胞
              "DCN","COL1A1","COL3A1","TAGLN",#成纤维
              "CDH5","VWF","ENG","CLDN5",#内皮
              "CHGA", "SYP","PHOX2A","PHOX2B",#神经内分泌细胞
              "SOX10","S100B","PLP1","EGR2"#施旺细胞
)
# Calculate average expression and percent detection
scRNA.sg2 <- AverageExpression(scRNA.sg, assay = "RNA", slot = "data", features = Markergene)
library(ggplot2)
cluster_order <- c("Myeloid_cells","DCs","T_cells","B_cells","NK_cells","Fibroblasts","Endothelials","NEs","SCs")
scRNA.sg$Annotation_Level2 <- factor(scRNA.sg$Annotation_Level2, levels = cluster_order)
any(duplicated(cluster_order))
all(cluster_order %in% levels(scRNA.sg$Annotation_Level2))
setdiff(levels(scRNA.sg$Annotation_Level2), cluster_order)
DotPlot(scRNA.sg, features = Markergene, group.by = "Annotation_Level2", dot.scale = 6, cols = c("lightgrey", "blue")) +
  theme_minimal() +
  theme(#这一步耗时比较长，建议用Rscript
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  ) +
  scale_color_viridis_c()
# 绘制 DotPlot
DotPlot(scRNA.sg, features = Markergene, group.by = "Annotation_Level2", dot.scale = 6, cols = c("lightgrey", "blue")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  )
# 绘制 DotPlot 使用 viridis 颜色刻度
DotPlot(scRNA.sg, features = Markergene, group.by = "Annotation_Level2", dot.scale = 6) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  ) +
  scale_color_viridis_c()
# Print the DotPlot
print(Dotplot)
#更新一下这个最终注释后的数据
scRNA.sg<-qread("~/scRNA/R_Project/Schwann_project/Process_Data/GSE137804/scRNA.sg_Final_annotation.qs")
rm(list = ls())
gc()

#################----六、开始进行细胞分组和差异分析----###################
### 1. 细胞分组 ###

scRNA.sg<-qread("~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_sg.qs")

### 2.1 走一遍标准流程，把肿瘤细进行展示
#肿瘤细胞
scRNA.mal<-scRNA.sg[,scRNA.sg$Annotation_Level2 %in% "tumor"]
scRNA.mal <- NormalizeData(scRNA.mal) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA.mal <- RunPCA(scRNA.mal,npcs = 50,seed.use = 240701292) %>% RunHarmony(group.by.vars = "orig.ident",max.iter.harmony = 10)
ElbowPlot(scRNA.mal,ndims = 50)
pc.num=1:30
scRNA.mal <- FindNeighbors(scRNA.mal,dims = pc.num,reduction = "harmony")
scRNA.mal <- FindClusters(scRNA.mal, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA.mal)
table(scRNA.mal$RNA_snn_res.0.1)
sel.clust = "RNA_snn_res.0.1"
scRNA.mal <- SetIdent(scRNA.mal, value = sel.clust)
scRNA.mal = RunUMAP(scRNA.mal,reduction = "harmony",dims = pc.num)
DimPlot(scRNA.mal,group.by = "orig.ident")+DimPlot(scRNA.mal,group.by = "RNA_snn_res.0.1")
scRNA.mal$seurat_clusters<-scRNA.mal$RNA_snn_res.0.1
DimPlot(scRNA.mal)

### 2.2 可视化

table(scRNA.mal$cell_type)
names(scRNA.mal@meta.data)
names(scRNA.sg@meta.data)
scRNA.mal$cell_type<-paste0("NB_C",scRNA.mal$seurat_clusters)
meta.NB<-as.data.frame(scRNA.mal$cell_type)
colnames(meta.NB)<-"cell_type"
meta.NB$Barcode<-rownames(meta.NB)
table(scRNA.sg$Annotation_Level2)
meta.non_NB<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Final %in% "tumor")]$Annotation_Final)
colnames(meta.non_NB)<-"cell_type"
table(meta.non_NB$cell_type)
head(meta.non_NB)
meta.non_NB$Barcode<-rownames(meta.non_NB)

#合并并组合到scRNA.sg
meta.all<-rbind(meta.non_NB,meta.NB)
meta.all<-meta.all[colnames(scRNA.sg),]
scRNA.sg$cell_type<-meta.all$cell_type
table(scRNA.sg$cell_type,scRNA.sg$Annotation_Level2)
table(scRNA.sg$cell_type)




#Schwann_cells
scRNA_SCs<-scRNA.sg[,scRNA.sg$Annotation_Final %in% "SCs"]
scRNA_SCs <- NormalizeData(scRNA_SCs) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
scRNA_SCs <- RunPCA(scRNA_SCs,npcs = 50,seed.use = 240701292) %>% RunHarmony(group.by.vars = "Sample",max.iter.harmony = 10)
ElbowPlot(scRNA_SCs,ndims = 50)
pc.num=1:20
scRNA_SCs <- FindNeighbors(scRNA_SCs,dims = pc.num,reduction = "harmony")
scRNA_SCs <- FindClusters(scRNA_SCs, resolution = seq(0,1,by=0.1),random.seed = 240701292)
library(clustree)
clustree(scRNA_SCs)
table(scRNA_SCs$RNA_snn_res.0.1)
sel.clust = "RNA_snn_res.0.1"
scRNA_SCs <- SetIdent(scRNA_SCs, value = sel.clust)
scRNA_SCs = RunUMAP(scRNA_SCs,reduction = "harmony",dims = pc.num)
DimPlot(scRNA_SCs,group.by = "cell_type")+DimPlot(scRNA_SCs,group.by = "RNA_snn_res.0.1")
scRNA_SCs$seurat_clusters<-scRNA_SCs$RNA_snn_res.0.1
DimPlot(scRNA_SCs)

### 2.2 可视化
table(scRNA_SCs$cell_type)
names(scRNA_SCs@meta.data)
names(scRNA.sg@meta.data)
scRNA.mal$SCs_Subtype<-paste0("SCs_C",scRNA.mal$seurat_clusters)
meta.SCs<-as.data.frame(scRNA.mal$SCs_Subtype)
colnames(meta.SCs)<-"SCs_Subtype"
meta.SCs$Barcode<-rownames(meta.SCs)
table(scRNA.sg$Annotation_Level2)
meta.non_SCs<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Final %in% "SCs")]$Annotation_Final)
colnames(meta.non_SCs)<-"SCs_Subtype"
table(meta.non_SCs$SCs_Subtype)
meta.non_SCs$Barcode<-rownames(meta.non_SCs)
head(meta.non_SCs)
#合并并组合到scRNA.sg
meta.all<-rbind(meta.non_SCs,meta.SCs)
meta.all<-meta.all[colnames(scRNA.sg),]
scRNA.sg$SCs_Subtype<-meta.all$SCs_Subtype
table(scRNA.sg$SCs_Subtype,scRNA.sg$Annotation_Final)
table(scRNA.sg$Sample_ID,scRNA.sg$SCs_Subtype)
table(scRNA.sg$SCs_Subtype)

qsave(scRNA.sg,file = "~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA.sg.qs")
qsave(scRNA.mal,"~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_NB.qs")
qsave(scRNA_SCs,"~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_SCs.qs")
rm(list=ls())
gc()


###对施旺细胞进行NMF分群
library(NMF)
library(dplyr)
library(data.table)
library(clusterProfiler)
library(Seurat)

scRNA_SCs<-scRNA.sg[,scRNA.sg$Annotation_Final %in% "SCs"]
#需要首先对单核单细胞亚群矩阵进行归一化
sub_sce=CreateSeuratObject(
  counts = scRNA_SCs@assays$RNA@counts,
  meta.data = scRNA_SCs@meta.data
)

Idents(sub_sce)<-sub_sce$Annotation_Final
table(Idents(sub_sce))
sub_sce$celltype = Idents(sub_sce)
p2=DimPlot(sub_sce,label = T)
library(patchwork)
p1+p2


sub_sce = NormalizeData(sub_sce) %>% FindVariableFeatures() %>% ScaleData(do.center = F)

##非负矩阵分解分析
vm <- sub_sce@assays$RNA@scale.data
result <- nmf(vm,
              rank = 2:10,
              method = "lee",
              nrun = 100,
              seed = 240701292)

plot(result)
res <- nmf(vm,
           rank = 3,
           method = "lee",
           nrun = 100,
           seed = 240701292)


index <- extractFeatures(res,"max") 
sig.order <- unlist(index)
NMF.rank3 <- vm[sig.order,]
NMF.rank3 <- na.omit(NMF.rank3)
group <- data.frame(predict(res))
table(group)
head(group)

group <- predict(res) # 提出亚型
table(group)
sample.order <- names(group[order(group)])
#设置颜色
jco <- c("#2874C5","#EABF00","#C6524A","#868686")

pdf(file = "consensusmap.pdf")
consensusmap(res,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(vm)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3])))
dev.off()

pdf(file = "basismap.pdf" ) 
# 从这张图可以比较清晰地看到各亚型中的驱动signature（深色），与下面的nmf heatmap对应

basismap(res,
         cexCol = 1,
         cexRow = 0.3,
         annColors=list(c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
dev.off()
pdf(file = "NMF_heatmap.pdf" )
aheatmap(as.matrix(vm[,sample.order]), 
         Rowv=NA, 
         Colv=NA, 
         annCol = data.frame("cluster"=group[sample.order]),
         annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3])),
         color=c("#EAF0FA","#6081C3","#3454A7"), # 例文的蓝色渐变
         revC=TRUE, 
         cexCol = 0.3, 
         cexRow = 0.3,
         filename = "NMF_heatmap.pdf")
dev.off()

ac=data.frame(group=paste0('SCs_C',group))
rownames(ac)=colnames(vm)
library(pheatmap)
pdf("NMF_cluster.pdf")
pheatmap(vm,show_colnames = F,annotation_col = ac)
dev.off()
pdf("NMF_cluster.1.pdf")
pheatmap(vm[,sample.order],
         show_colnames = F,cluster_cols = F,annotation_col = ac)
dev.off()

save(result,file = "./NMF/NMF_2-10_result.rda")
save(res,file = "./NMF/NMF_res3.rda")







## 选择后续分析的因子，使用 NMF 运行的结果进行降维和聚类 
sub_sce
sub_sce <- RunPCA(sub_sce,verbose = F)
sub_sce@reductions$nmf <- sub_sce@reductions$pca
sub_sce@reductions$nmf@cell.embeddings <- t(coef(res) )
sub_sce@reductions$nmf@feature.loadings <- basis(res) 
sub_sce <- RunUMAP(sub_sce,reduction = "nmf",dims = 1:2) 
##基于非负矩阵分解后的结果进行分群
sub_sce <- FindNeighbors(sub_sce,reduction = "nmf",dims = 1:3) %>% FindClusters(resolution = 0.1)

sub_sce$cluster <- apply(NMF::coefficients(res)[1:3,],2,which.max)

table( sub_sce$Annotation_Final ,sub_sce$cluster)
table(Idents(sub_sce) ,sub_sce$cluster)
DimPlot(sub_sce,group.by = "cluster")
sub_sce$cluster<-paste0("SCs_C",sub_sce$cluster)
table(sub_sce$cluster)
qsave(sub_sce,file = "~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_SCs_NMF.qs")

sub_sce<-qread("~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_SCs_NMF.qs")

###将NMF分群并入所有细胞中
meta.un<-as.data.frame(sub_sce$cluster)
colnames(meta.un)<-"Annotation"
meta.un$Barcode<-rownames(meta.un)
meta.unun<-as.data.frame(scRNA.sg[,!(scRNA.sg$Annotation_Final %in% "SCs")]$Annotation_Final)
colnames(meta.unun)<-"Annotation"
meta.unun$Barcode<-rownames(meta.unun)
meta<-rbind(meta.unun,meta.un) 
meta<-meta[colnames(scRNA.sg),]
scRNA.sg$SCs_Final_cluster<-meta$Annotation

table(scRNA.sg$SCs_Final_cluster)
table(scRNA.sg$Annotation_Final,scRNA.sg$SCs_Final_cluster)
table(scRNA.sg$NMF_cluster,scRNA.sg$SCs_Final_cluster)
colnames(scRNA.sg@meta.data)[45]<-"NEs_Seurat_clusters"
colnames(scRNA.sg@meta.data)[46]<-"SCs_Seurat_clusters"
colnames(scRNA.sg@meta.data)[47]<-"cell_type"
qsave(scRNA.sg,file = "~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_sg.qs")



### 3. 开始寻找差异基因 ###
######寻找差异基因
table(scRNA.mal$cell_type)
scRNA.mal <- SetIdent(scRNA.mal, value = "cell_type")
mal.markers <- FindAllMarkers(object = scRNA.mal,
                              only.pos = FALSE,
                              min.pct = 0.25, 
                              thresh.use = 0.25) 
#把差异基因拿去cmap预测药物
write.table(mal.markers,file = "~/scRNA/R_Project/Atherosclerosis/scRNA_analysis/middle_data/SMC/SMC_DEG.txt",quote = F,sep = "\t")
#scRNA.neutro <- SetIdent(scRNA.neutro, value = "Article_Group")
#neotro.markers <- FindAllMarkers(object = scRNA.neutro,
#                                 only.pos = TRUE,
#                                 min.pct = 0.25, 
#                                 thresh.use = 0.25)
scRNA.mal@assays$RNA@var.features<-mal.markers$gene
scRNA.mal<-ScaleData(scRNA.mal)
#scRNA.neutro@assays$RNA@var.features<-neotro.markers$gene
#scRNA.neutro<-ScaleData(scRNA.neutro)
#热图
DoHeatmap(scRNA.mal,features = mal.markers$gene)+NoLegend()+scale_fill_gradientn(colors = c("navy","white","firebrick3"))
#DoHeatmap(scRNA.neutro,features = neotro.markers$gene)+NoLegend()+scale_fill_gradientn(colors = c("navy","white","firebrick3"))

#################----七、开始进行功能富集分析----###################
scRNA.mal<-qread("../../Analysis_data/scRNA_SCs.qs")
#saveRDS(scRNA.neutro,"~/R_project/PDAC/Article/Intermediate/scRNA.netro.RDS")
### 1. 跳转irGSEA_GO_Mal.R等等进行功能富集分析 ###
### 2. 跳转7. irGSEA_Heatmap.R进行结果可视化 ###




#################----八、开始进行细胞交互分析----###################
qsave(scRNA.sg,"~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA.sg.qs")
### 1. 跳转8. CellChat_Run.R进行细胞交互分析计算与可视化
### 2. 对运行结果进一步可视化
library(CellChat)
library(ComplexHeatmap)
library(patchwork)
load("~/scRNA/R_Project/Schwann_project/scRNA_analysis/tmp/CellChat/SCs_NMF_cluster_CellChat/SCs_cellchat_before_cycle.Rda")
glimpse(cellchat)
table(cellchat@meta$Article_Group)
netVisual_heatmap(cellchat, measure = "weight",color.heatmap = c("#DDEAF0","#6731CC"),#color.use = my43color,
                  #targets.use = names(table(cellchat@meta$Article_Group)[c(10:16)]), 
                  #sources.use = names(table(cellchat@meta$Article_Group)[c(2:4,7:9,20:22)])#,
                  #                  row.show = c(1:9,17:22),
                  #                  col.show =  c(10:16)
)
netAnalysis_signalingRole_heatmap(cellchat,
                                  pattern = "incoming", #传入
                                  signaling = cellchat@netP$pathways,
                                  #title = "Outgoing",
                                  width = 8,#color.heatmap = "BuPu",
                                  height = 10)
netAnalysis_signalingRole_heatmap(cellchat,
                                  pattern = "outgoing", #传出
                                  signaling = cellchat@netP$pathways,
                                  #title = "Outgoing",
                                  width = 8,#color.heatmap = "BuPu",
                                  height = 10)

#################----九、开始进行细胞拟时序分析----###################
setwd("~/scRNA/R_Project/Schwann_project/scRNA_analysis/tmp/monocle2/")
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle)
library(biomaRt)
library(SeuratObject)
library(qs)
### 1. 输入数据准备 ###
scRNA.mal<-qread("~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_SCs.qs")
expr_matrix<-as(as.matrix(scRNA.mal@assays$RNA@counts),"sparseMatrix")
DefaultAssay(scRNA.mal)<-"RNA"
p_data<-scRNA.mal@meta.data
f_data<-data.frame(gene_short_name = row.names(scRNA.mal),row.names = row.names(scRNA.mal))

#构建CDS对象
pd <-new("AnnotatedDataFrame",data = p_data)
fd<-new("AnnotatedDataFrame",data = f_data)
cds<-newCellDataSet(expr_matrix,
                    phenoData = pd,
                    featureData = fd,
                    lowerDetectionLimit = 0.5,
                    expressionFamily = negbinomial.size())#适合UMI数据

### 2. 估计size factor和离散度 ###
cds <- estimateSizeFactors(cds)#帮助标准化mRNA的差异
cds <- estimateDispersions(cds)#帮助后续差异分析，类似归一化处理

### 3. 轨迹定义基因选择及可视化和构建轨迹 ###
#选择高变基因进行轨迹分析
scRNA.mal<-FindVariableFeatures(scRNA.mal,nfeatures = 2000)
ordergenes<-scRNA.mal@assays$RNA@var.features
cds<-setOrderingFilter(cds,ordergenes)
plot_ordering_genes(cds)
#该图黑点表示构建轨迹的差异基因，灰色背景基因，红色为gene.size和离散度分布趋势，可以看到找到的基因为离散度比较高的基因

### 4. CDS数据降维 ###
#耗时，用Rscript跑
save(cds,file="./SCs_cds.rda")
cds <- reduceDimension(cds,max_components = 2,method="DDRTree",)#可以额外接受sigma参数调控分叉，可以改为0.1减少分叉
save(cds,file = "./SCs_cds_afterRed.Rda")
load("~/scRNA/R_Project/Schwann_project/scRNA_analysis/tmp/monocle2/SCs_cds_afterRed.Rda")
#拟时间轴轨迹构建、排列细胞
cds <- orderCells(cds)
#降低版本后碰到Error in get.diameter(minSpanningTree(cds)) : ，安装igraph即可
#如果无法安装igraph就手动到国内镜像下载glpk：https://mirrors.nju.edu.cn/gnu/glpk/
#随后手动configure, make, make install即可
#随后在R中添加路径：Sys.setenv(LD_LIBRARY_PATH = paste("/home/liangxs/software/lib", Sys.getenv("LD_LIBRARY_PATH"), sep=":"))
#Sys.setenv(LD_LIBRARY_PATH = paste("/home/liangxs/software/lib", Sys.getenv("LD_LIBRARY_PATH"), sep=":"))
#Sys.setenv(PKG_LIBS = "-L/home/liangxs/software/lib")
#日内瓦，发现packages里面是有igraph的，直接载入就行，
library(igraph)
#Error in reducedDimK(cds) <- K_old : 没有"reducedDimK<-"这个函数
#解决：order_cells.R这文件中的1158行（如果您是从该博主这下载的）#reducedDimK(cds) <- K_old注释掉即可，然后重新source('order_cells.R')；我这里还报了reducedDimW(cds) <- old_W这个错误，同样的1162行#reducedDimW(cds) <- old_W注释掉，再重新source('order_cells.R')
#打开函数源码编辑界面(或者选中函数按F2)
source('/home/liangxs/R/x86_64-pc-linux-gnu-library/4.2/monocle_2.21_corrected/R/order_cells.R')
#trace('project2MST',edit = T,where = asNamespace("monocle"))
#打开编辑界面，然后删除报错的行，save就行
#总结一下monocle斗智斗勇的过程：
#1. 下载版本为2.21.1
#2. 去order_cell.R修改：! (projection)..., reducedDimK<-, reducedDimK(cds) <- K_old
#3. 去BEAM.R、utils.R分别修改对应的函数
#4. 去expr_models.R 把153行和155行x@family@vfamily....改为any(x@expressionFamily@vfamil...
#5. 加载igraph
source('/home/liangxs/R/x86_64-pc-linux-gnu-library/4.2/monocle_2.21_corrected/R/expr_models.R')
source('/home/liangxs/R/x86_64-pc-linux-gnu-library/4.2/monocle_2.21_corrected/R/BEAM.R')
### 5. 拟时序轨迹可视化 ###
library(ggsci)
library(RColorBrewer)
# 假设你有9个独特的状态
n_states <- 5
color9 <- brewer.pal(n_states, "Set3") 

plot_cell_trajectory(cds,color_by = "Pseudotime",size=1,show_backbone = T)
plot_cell_trajectory(cds,color_by = "State",size=1,show_backbone = T)+scale_color_manual(values = color9)
### 5.1 拟时序轨迹
color3 <- c("#F5814C", "#87D6C2", "#6A5ACD")
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = T)+scale_color_manual(values = color3)+
  facet_wrap("~cell_type",nrow=1)+
  scale_color_manual(values = color3)

### 5.2 拟时序峰值图
library(ggpubr)
df<-pData(cds)
#View(df)
ggplot(df,aes(Pseudotime,colour = cell_type, fill = cell_type)) +
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2() + scale_color_manual(values =color3) + 
  scale_fill_manual(values = color3)

### 5.3 单个基因表达轨迹展示
cds_subset<-cds[,cds$Article_Group %in% "Mal_SLC39A4_High_GSDMD_High"]
cds_subset<-cds_subset[c("SLC39A4","GSDMD"),]
#可视化：以group/Type展开
plot_genes_in_pseudotime(cds_subset = cds_subset,color_by = "Article_Group")+ scale_color_manual(values = color4)
plot_genes_violin(cds_subset = cds_subset,color_by = "Article_Group")+ scale_color_manual(values = color2)
plot_genes_jitter(cds_subset = cds_subset,color_by = "Article_Group")+ scale_color_manual(values = color4)
plot_cell_trajectory(cds,color_by = "SLC39A4",) + scale_color_gradient("low" = color5_g[1],"high" = color5_g[5])

### 5.4 时序关键基因分析 ###
Time_diff<-differentialGeneTest(cds[ordergenes,],
                                fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                #reducedModelFormulaStr = "~1", relative_expr = TRUE, cores = 4,
                                verbose = FALSE)

Time_genes <-Time_diff[Time_diff$qval<0.05,] %>% pull(gene_short_name) %>% as.character()

save(Time_diff,Time_genes,file = "./Time_diff_time.rda")###后续出现报错
p=plot_pseudotime_heatmap(cds[1:100,],#num_clusters = 4,
                          show_rownames = F,cores = 10,return_heatmap = T#,
                          #hmcols = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(100)
)
source("/home/liangxs/R/x86_64-pc-linux-gnu-library/4.2/monocle_2.21_corrected/R/BEAM.R")
BEAM_res_1 <- BEAM(cds[ordergenes,],
                   branch_states=c(1,2,3),
                   #branch_point = 1,
                   cores = 3,fullModelFormulaStr = "~sm.ns(Pseudotime)" )
myColors <- c("#407BD0", "#B7D0EA", "#E3EAF1", "#F4E7E8", "#F4E7E8", "#DB464E", "#DB464E", "#A32A31")
myColorPalette <- colorRampPalette(myColors)(100)

p1<-plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res_1,qval < 0.01)),],
                                branch_point = 1,
                                branch_labels = c(2,1),
                                return_heatmap = T,
                                num_clusters = 5,
                                cores = 1,hmcols = colorRampPalette(c("#10AEC2","white","#EE3F4D","#DC143C"))(100),
                                use_gene_short_name = T)
#遇到了Monocle无法解决的保存问题，the condition has lenght > 1

p1$ph_res$tree_row
clusters<-cutree(p1$ph_res$tree_row,k=5)
clusters <- data.frame(clusters)
clusters[,1] <- as.character(clusters[,1])
colnames(clusters) <- "Gene_Clusters"
table(clusters)
write.csv(clusters,"~/R_project/TMZ_sc_GBM40/Intermediate/State5_Clusters.csv",row.names = F)
#
branch_deg<-p1$ph_res$tree_row$labels
clusters<-rownames_to_column(clusters,var="ID")
branch_deg_up<-branch_deg[as.numeric(rownames(clusters[clusters$Gene_Clusters %in% c("1"),]))]
library(tidyverse)
library(dbplyr)
library(biomaRt)
library(qs)
library(clusterProfiler)

qsavem(cds,p1,branch_deg_up,file="~/R_project/PDAC/Article/Intermediate/tmp_use.qs")
rm(cds)
qreadm(file="~/R_project/PDAC/Article/Intermediate/tmp_use.qs")
Mart<-useMart("ensembl")
Mart_db<-listDatasets(Mart)
mart_oas=useMart("ensembl", "hsapiens_gene_ensembl")
Mono_branch_symbol=getBM(attributes = c("hgnc_symbol","ensembl_gene_id","entrezgene_id","description","chromosome_name"),
                         filters = "hgnc_symbol",
                         values = branch_deg_up,
                         mart = mart_oas)
Mono_branch_symbol<-Mono_branch_symbol[,c(1,3)]
Mono_branch_symbol[1:4,1:2]

#使用气泡图进行功能富集分析
require(DOSE)
require(clusterProfiler)
require(org.Hs.eg.db)
enrich_gobp <- enrichGO(gene=Mono_branch_symbol$entrezgene_id,'org.Hs.eg.db',
                        ont="BP",
                        pvalueCutoff=0.05,
                        readable=TRUE)
Time_diff
qsavem(Mono_branch_symbol,enrich_gobp,mart_oas,enrich_kegg,file = "~/R_project/PDAC/Article/Intermediate/tmp2.qs")
qreadm(file = "~/R_project/PDAC/Article/Intermediate/tmp2.qs")
#kegg
enrich_kegg <- enrichKEGG(gene = Mono_branch_symbol$entrezgene_id,
                          
                          organism = "hsa",
                          
                          pAdjustMethod = 'BH',
                          
                          pvalueCutoff=0.05,
                          
                          qvalueCutoff = 0.05)
dotplot(enrich_kegg,showCategory=100)
dotplot(enrich_gobp,showCategory=30)
GOBP_Category<-c("cell-substrate adhesion","leukocyte migration","extracellular matrix organization",
                 "ERK1 and ERK2 cascade","response to tumor necrosis factor","cellular response to transforming growth factor beta stimulus",
                 "integrin-mediated signaling pathway","aerobic electron transport chain","granulocyte chemotaxis","ATP synthesis coupled electron transport")

KEGG_Category<-c("Focal adhesion","Regulation of actin cytoskeleton","Oxidative phosphorylation",
                 "Phagosome","ECM-receptor interaction","TNF signaling pathway","Antigen processing and presentation","Relaxin signaling pathway")
dotplot(enrich_kegg,showCategory = KEGG_Category)+ scale_color_continuous(low='#EE3F4D', high='#10AEC2')
dotplot(enrich_gobp,showCategory = terms)+ scale_color_continuous(low='#EE3F4D', high='#10AEC2')

#整一手差异基因
library(Libra)
#install.packages('TMB', type = 'source') #弹出窗口提示选NO，重新安装后解决问题
#Idents(scRNA.mal) <- "NMF_Group" 
scRNA.mal@active.ident
#准备输入文件：矩阵、metadata
expr<-as.matrix(scRNA.mal@assays$RNA@data)
meta<-scRNA.mal@meta.data
meta<-data.frame(cell_type = meta$Article_Group,
                 replicate=meta$orig.ident,
                 label = meta$Annotation_Level1)
#cell_type, replicate, and label
#####开始比较
#由于Libra只能做两两比较，故分别做一组与其他两组的对比，参考https://github.com/neurorestore/Libra/issues/6
meta$cell_type<-ifelse(meta$cell_type %in% "Mal_SLC39A4_High_GSDMD_High","High","Low")
DE_PseudoLimma_TF = run_de(expr,meta,cell_type_col = "cell_type",label_col = "cell_type",
                           replicate="replicate",de_family = "pseudobulk",de_method = "edgeR",de_type = "LRT")
meta$cell_type<-factor(meta$cell_type,ordered = F)
library(Seurat)
table(scRNA.mal@active.ident)
allmarkers<-FindAllMarkers(scRNA.mal,logfc.threshold = 0,min.pct = 0,
                           min.diff.pct = 0,
                           only.pos = T)
sc_df<-allmarkers[allmarkers$cluster %in% "Mal_SLC39A4_High_GSDMD_High",c(7,2)]
colnames(sc_df)<-c("ID","logFC")

enrichdata<-enrich_gobp@result
enrichdata$category<-rep("GO",nrow(enrichdata))
colnames(enrichdata)[c(2,5,8)]<-c("term","adj_pval","genes")
enrichdata<-enrichdata[,c("category","ID","term","genes","adj_pval")]
enrichdata$genes<-gsub("/",", ",enrichdata$genes)
library(GOplot)
circ <- circle_dat(enrichdata, sc_df)
terms<-c("nitric oxide mediated signal transduction", "response to peptide hormone", "granulocyte migration", "alpha-amino acid metabolic process", "cellular metal ion homeostasis", "regulation of lipase activity" )
circ<-circ[circ$adj_pval < 0.05,]
circ<-circ[circ$term %in% terms,]
GOCircle(circ,nsub = unique(circ$ID),title = "Branch1",
         lfc.col = c("#00BFFF", "#FF4500"),zsc.col =  c("#FF4500","white","#00BFFF"))

#################----十、开始进行CytoTRACE评估分化潜能----###################
### 1. 跳转到CytoTRACE_Run
### 2. 读取结果进行可视化
library(CytoTRACE) #加载包
setwd("~/scRNA/R_Project/Schwann_project/scRNA_analysis/tmp/CytoTRACE/")
scRNA.mal<-qread("~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_SCs.qs")
CT_results<-qread("./SCs_CytoTRACE_results.qs")
load("../monocle2/SCs_cds_afterRed.Rda")
phe<-scRNA.mal@meta.data$cell_type
names(phe)<-colnames(scRNA.mal)
#将cytoTRACE用拟时序的坐标展示出来
cds@reducedDimS
cds.emb<-as.data.frame(t(cds@reducedDimS))
colnames(cds.emb)=c("component1","component2")
#定量比较两组的cytoTRACE order分数
cyto_order<-as.data.frame(CT_results$CytoTRACE)
cyto_order$group<-scRNA.mal@meta.data$cell_type
cyto_order$Barcode<-colnames(scRNA.mal)
colnames(cyto_order)[1]<-"CytoTRACE_order"
cyto_order<-pivot_wider(cyto_order,names_from = "group",values_from = "CytoTRACE_order")
cluster1_order<-na.omit(cyto_order$SCs_C1)
cluster2_order<-na.omit(cyto_order$SCs_C2)
cluster3_order<-na.omit(cyto_order$SCs_C3)
wilcox.test(cluster1_order,cluster2_order,paired = F,alternative = "two.sided",exact = F, correct = F)

plotCytoTRACE(CT_results, phenotype = phe, emb = cds.emb, outputDir = "./")
#可视化与 CytoTRACE 相关的基因
plotCytoGenes(CT_results, numOfGenes = 10,outputDir = "./")


positive_gene<-c("TPT1","RPS18","ATP5F1B","PPP1R14B","IGFBP5","RPL13","JUND","RACK1","RPL7A","GNAS")
negative_gene<-c("ITGB8","LGI4","SEMA3B","HLA−DRA","RARRES2","HLA−DPA1","HLA−DRB1","S100B","HLA−DRB5","CD74")
gene_list <- list(positive_gene = positive_gene, negative_gene = negative_gene)
scRNA.mal<-AddModuleScore(
  object = scRNA.mal,
  features = gene_list,
  #ctrl = 100,
  name = 'Features'
)
names(scRNA.mal@meta.data)[32]<-"Negative_Gene"
table(scRNA.mal$cell_type)
plots_list <- list()
for (i in names(scRNA.mal@meta.data[31:32])) {
  p <- VlnPlot(
    scRNA.mal,
    features =i,
    group.by ="cell_type",
    pt.size = 0,
    raster = F
  )+ geom_boxplot(width = .25, fill = "white")+ xlab("") + ylab("hME") + NoLegend()
  plots_list[[i]] <- p
}
# 使用patchwork组合所有图形
combined_plot <- wrap_plots(plots_list, ncol = 2) # 指定每行显示的图形数
combined_plot




rm(list = ls())
gc()

#################----十一、开始进行Tricycle评估细胞周期----###################
#ndevtools::install_github("jokergoo/circlize")
library(tricycle)
library(ggplot2)
library(scattermore)
library(scater)
scRNA.mal<-qread("~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_SCs.qs")
DefaultAssay(scRNA.mal)<-"RNA"
#data(cycle_obj)
### 1. 将Seurat对象转变为singlecellexperiment对象
cycle_obj<-SingleCellExperiment(
  assays = list(
    counts = as.matrix(scRNA.mal@assays$RNA@counts),
    logcounts = as.matrix(scRNA.mal@assays$RNA@data)
  ),
  colData = scRNA.mal@meta.data,
  rowData = rownames(scRNA.mal)
)
#预处理
cycle_obj <- project_cycle_space(cycle_obj,
                                 exprs_values = "logcounts",
                                 species = "human",
                                 gname.type = "SYMBOL"
)
### 2. 画个质控图
scater::plotReducedDim(cycle_obj, dimred = "tricycleEmbedding") +
  labs(x = "PC1", y = "PC2") +
  ggtitle(sprintf("Projected cell cycle space (n=%d)",
                  ncol(cycle_obj))) +
  theme_bw(base_size = 14)

#推断细胞周期位置
cycle_obj <- estimate_cycle_position(cycle_obj)
names(colData(cycle_obj))
#> [1] "TotalUMIs"        "sample"           "tricyclePosition"
top2a.idx <- which(rowData(cycle_obj)$Gene == 'TOP2A')
fit.l <- fit_periodic_loess(cycle_obj$tricyclePosition,
                            assay(cycle_obj, 'logcounts')[top2a.idx,],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                            fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                               ncol(cycle_obj), ")"))
names(fit.l)




top2a.idx <- which(rowData(cycle_obj)$X == 'PDLIM5')
fit.l <- fit_periodic_loess(cycle_obj$tricyclePosition,
                            assay(cycle_obj, 'logcounts')["PDLIM5",],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(PDLIM5)",
                            fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                               ncol(cycle_obj), ")"))
names(fit.l)
#> [1] "fitted"   "residual" "pred.df"  "loess.o"  "rsquared" "fig"
fit.l$fig + theme_bw(base_size = 14)

#Alternative: Infer cell cycle stages
cycle_obj <- estimate_Schwabe_stage(cycle_obj,
                                    gname.type = 'SYMBOL',
                                    species = 'human')
scater::plotReducedDim(cycle_obj, dimred = "tricycleEmbedding",
                       colour_by = "cell_type") +
  labs(x = "Projected PC1", y = "Projected PC2",
       title = paste0("Projected cell cycle space (n=", ncol(cycle_obj), ")")) +
  theme_bw(base_size = 14)


plot_ccposition_den(cycle_obj$tricyclePosition,
                    cycle_obj$cell_type, 'cell_type',
                    bw = 10, fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)

plot_ccposition_den(cycle_obj$tricyclePosition,
                    cycle_obj$cell_type, 'cell_type', type = "circular",
                    bw = 10,  fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)  

plot_ccposition_den(cycle_obj$tricyclePosition,
                    cycle_obj$seurat_clusters,'orig.ident',#color_name  = color,
                    bw = 10, fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)

plot_ccposition_den(cycle_obj$tricyclePosition,
                    cycle_obj$seurat_clusters, 'orig.ident', type = "circular",
                    bw = 10,  fig.title = "Kernel density of \u03b8") +
  theme_bw(base_size = 14)  

#Plot out embedding scater plot colored by cell cycle position
library(cowplot)
p <- plot_emb_circle_scale(cycle_obj, dimred = 1,
                           point.size = 3.5, point.alpha = 0.9) +
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))
p
#save(cycle_obj,file = "~/R_project/NDC80/Intermediate/Tricycle_obj.Rda")
#rm(cycle_obj)
saveRDS(fit.l,file = "~/scRNA/R_Project/Atherosclerosis/scRNA_analysis/tmp/Tricycle/Tricycle_fit.RDS")
qsave(cycle_obj,file =  "~/scRNA/R_Project/Schwann_project/scRNA_analysis/tmp/Tricycle/Tricycle_obj.qs")
rm(list = ls())
gc()

#########准备跑空间
##将数据保存为Scanpy输入
library(SeuratDisk)
scRNA.sg<-qread("~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/scRNA_sg.qs")
table(scRNA.sg$cell_type)
scRNA.scanpy <- CreateSeuratObject(scRNA.sg@assays$RNA@counts,meta.data = scRNA.sg@meta.data)
#1. 去掉scale.data
scRNA.scanpy<-DietSeurat(
  scRNA.sg,
  counts = TRUE,
  #data = TRUE,
  scale.data = FALSE,
  features = rownames(scRNA.sg),
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)
#2. # step 2: factor to character, or else your factor will be number in adata
i <- sapply(scRNA.scanpy@meta.data, is.factor)
scRNA.scanpy@meta.data[i] <- lapply(scRNA.scanpy@meta.data[i], as.character)

# step 3: convert
setwd("~/scRNA/R_Project/Schwann_project/scRNA_analysis/GSE137804/scanpy_data/")
setwd("~/scRNA/R_Project/Schwann_project/scRNA_analysis/Analysis_data/Scanpy_data/")
SaveH5Seurat(scRNA.scanpy, filename = "scRNA_sg_count.h5seurat", overwrite = TRUE)
Convert("scRNA_sg_count.h5seurat", "scRNA_sg_count.h5ad", assay="RNA", overwrite = TRUE)
rm(scRNA.scanpy)
DefaultAssay(scRNA.scanpy)<-"RNA"
