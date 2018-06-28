library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(clusterProfiler)
# load("~/Documents/cohort/ESCC/for_NEJM/nejm.rda")
# Comebine SNV results and CNV results of ESCC single-cell
# ======

# Reading SNV results
# ------
snv <- read.csv("~/Documents/cohort/ESCC/for_NEJM/sf9.csv",row.names=1)
# ------

# Reading CNV results
# ------
cnv_17 <- read.table("~/Documents/cohort/ESCC/for_NEJM/17cm_CNV_matrix_statistic_ann.txt",header=T,sep="\t")
cnv_23 <- read.table("~/Documents/cohort/ESCC/for_NEJM/23cm_CNV_matrix_statistic_ann.txt",header=T,sep="\t")
library(dplyr)
cnv_17 <- select(cnv_17,-therapy_increase)
cnv_23 <- select(cnv_23,-therapy_increase)
cnv <- cbind(select(cnv_17,1:3),select(cnv_17,4:28),select(cnv_23,4:48),
             select(cnv_17,primary_average_17=primary_average,therapy_average_17=therapy_average,wilcox_17=wilcox_p),
             select(cnv_23,primary_average_23=primary_average,therapy_average_23=therapy_average,wilcox_23=wilcox_p),
             select(cnv_17,cytoBand=an17.cytoBand,Genes=an17.Genes))
# ------

# Transforming CNV
# I transform the CNV from regions to cytoBands. The number is the count of 
# regions in one cytoBand.
# ------
temp <- cnv[,4:73]
temp <- temp[,match(colnames(snv),colnames(temp))]
temp <- apply(temp,2,function(x) sapply(x,function(x) {
  if (x > 2) return(1)
  if (x < 2) return(-1)
  if (x == 2) return(0)
  }))
temp <- cbind(cnv[,1:3],temp,cnv[,74:81])
cnv <- matrix(0,nrow=0,ncol=70)
temp <- sapply(levels(temp$cytoBand),function(i) {
  x <- filter(temp,cytoBand == i)
  y <- apply(x[,4:73],2,function(a) {
    if (length(a[a == 1]) > 0)
      j <- length(a[a == 1])
    else
      j <- 0
    if (length(a[a == -1]) > 0)
      k <- -length(a[a == -1])
    else
      k <- 0
    return(c(j,k))
  })
  rownames(y) <- c(paste0("+",i),paste0("-",i))
  colnames(y) <- colnames(temp)[4:73]
  cnv <<- rbind(cnv,y)
})
rm(temp)
# ------

# Heatmap of SNV
# ------
# Final version
temp <- snv[,26:70]
temp <- apply(temp,1,function(x) length(x[x != 0])) %>>% (temp[ . != 0,])
cluster <- t(temp[,1:30]) %>>% dist(method="euclidean") %>>% hclust(method="ward.D2")
temp <- temp[,1:30] %>>% (.[,cluster$order]) %>>% cbind(temp[,31:45])
cluster <- t(temp[,31:45]) %>>% dist(method="euclidean") %>>% hclust(method="ward.D2")
temp <- temp[,31:45] %>>% (.[,cluster$order]) %>>% (cbind(temp[,1:30],.))
# Clustering for rows
# Canberra distance may be better!
cluster <- dist(temp,method="euclidean") %>>% hclust(method="ward.D2")
temp <- temp[rev(cluster$order),]

### SNV table for single-cell 23cm
st3_snv <- read.table("~/Documents/cohort/ESCC/for_NEJM/st3_table.tsv",sep="\t",header=T,as.is=T)
st3_genes <- paste(st3_snv[,4],collapse=",") %>>% (strsplit(.,";",fixed=T)[[1]]) %>>% 
  paste(collapse=",") %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique
sc_variant <- st3_snv[grep(paste(intersect(rownames(temp),st3_genes),collapse="|"),st3_snv[,4]),]
sc_variant <- filter(sc_variant,sc_variant[,4] != "PRAMEF11")
temp <- temp[intersect(rownames(temp),st3_genes),]
t <- rep(NA,nrow(sc_variant))
t[grep(paste(genes_snv_23_single_cell_important,collapse="|"),as.character(sc_variant[,4]))] <- "*"
sc_variant <- cbind(sc_variant,important=t)
write.table(sc_variant,file="~/Documents/cohort/ESCC/for_NEJM/T1. Single-cell_variants.tsv",sep="\t",
            col.names=T,row.names=F,quote=F)
rm(t)
t <- bitr(genes_snv_23_single_cell_important,fromType="SYMBOL",toType="ENTREZID",
          OrgDb="org.Hs.eg.db")$ENTREZID %>>% enrichKEGG %>>% as.data.frame
t$geneID <- sapply(t$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
                     (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>%
                     paste(collapse="/"),simplify=T)
write.csv(t,file="~/Documents/cohort/ESCC/for_NEJM/T1. Single-cell_important_SNV_genes_KEGG.csv",quote=F)
rm(t)
###

# 1:75 are not all good!
# temp <- temp[1:75,]
temp <- temp[1:54,] %>>% (.[apply(temp[1:54,],1,function(x) x[31:45] %>>% 
                                    (length(.[. != 0]) <= 1)),]) %>>% 
  rbind(temp[55:75,])
# Different clusters, different colors
temp[,1:12] <- apply(temp[,1:12],2,function(x) {
  x[x != 0] <- 1
  return(x)
})
temp[,13:30] <- apply(temp[,13:30],2,function(x) {
  x[x != 0] <- 2
  return(x)
})
temp[,31:45] <- apply(temp[,31:45],2,function(x) {
  x[x != 0] <- 1
  return(x)
})
# p <- apply(temp,1,function(x) {
#   a <- x[1:30]
#   b <- x[31:45]
#   a <- length(a[a != 0])
#   b <- length(b[b != 0])
#   fisher.test(matrix(c(a,b,30-a,15-b),nrow=2))$p.val
# })
# rownames(temp)[match(names(p[p >= 0.05]),rownames(temp))] <- paste0(names(p[p >= 0.05])," *")
p <- apply(temp,1,function(x) {
  a <- x[1:30]
  b <- x[31:45]
  a <- length(a[a != 0])
  b <- length(b[b != 0])
  fisher.test(matrix(c(a,b,30-a,15-b),nrow=2),alternative="less")$p.val
})
row_labels <- paste0(names(p[p < 0.05])," **")
row_subsets <- match(names(p[p < 0.05]),rownames(temp))
rownames(temp)[match(names(p[p < 0.05]),rownames(temp))] <- paste0(names(p[p < 0.05])," **")
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",30),rep("Therapy_post",15)),
                        Sub_groups=anno_points(rep(1,45),
                                               pch=c(rep(7,12),rep(8,18),rep(7,15)),
                                               size=unit(12,"pt"),border=F,ylim=c(0.5,1.5),
                                               gp=gpar(col=c(rep("palevioletred",12),rep("cornflowerblue",18),
                                                             rep("palevioletred",15)))),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy_post"="#00DAE0")),
                        # annotation_height=c(0.5,0.5),
                        show_annotation_name=T,
                        # gap=unit(3,"pt"),
                        annotation_name_gp=gpar(fontsize=18,fontface="bold",fontfamily="serif"),
                        annotation_legend_param=list(title_gp=gpar(fontfamily="serif",fontface="bold",fontsize=26),
                                                     labels_gp=gpar(fontfamily="serif",fontsize=22),
                                                     grid_height=unit(10,"mm"),
                                                     grid_width=unit(10,"mm")))
row_ha <- rowAnnotation(link=row_anno_link(row_subsets,labels=row_labels,
                                           labels_gp=gpar(fontfamily="serif",fontsize=22)),
                        width=unit(10,"pt") + max_text_width(row_labels))
ht <- Heatmap(temp,col=colorRamp2(c(0,1,2),c("gray98","palevioletred","cornflowerblue"),space="sRGB"),show_row_names=T,
              show_column_names=F,row_names_gp=gpar(fontsize=16,fontfamily="serif"),cluster_columns=F,cluster_rows=F,
              top_annotation=ha,top_annotation_height=unit(30,"pt"),name="SNV",
              show_heatmap_legend=F,rect_gp=gpar(col="gray98")) #+
              # heatmap_legend_param=list(color_bar="continuous"))
  # row_ha
lgd <- Legend(c("A","B"),title="Sub_groups",type="points",pch=7:8,background="white",size=unit(6,"mm"),
              grid_height=unit(10,"mm"),grid_width=unit(10,"mm"),
              legend_gp=gpar(col=c("palevioletred","cornflowerblue")),
              title_gp=gpar(fontfamily="serif",fontface="bold",fontsize=26),
              labels_gp=gpar(fontfamily="serif",fontsize=22))
draw(ht,annotation_legend_list=list(lgd))
decorate_heatmap_body("SNV",{
  i <- c(12,30)
  x <- i/ncol(temp)
  grid.lines(c(x[1],x[1]),c(0,1),gp=gpar(lwd=2,col="red"))
  grid.lines(c(x[2],x[2]),c(0,1),gp=gpar(lwd=2,col="red"))
  grid.lines(c(0,1),c(21/nrow(temp),21/nrow(temp)),gp=gpar(lwd=2,col="green"))
})
rm(temp)
# ------

# Heatmap of CNV
# ------
cluster <- t(cnv[,1:20]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
cnv <- cnv[,1:20] %>>% (.[,cluster$order]) %>>% cbind(cnv[,21:70])
cluster <- t(cnv[,21:25]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
cnv <- cnv[,21:25] %>>% (.[,cluster$order]) %>>% (cbind(cnv[,1:20],.,cnv[,26:70]))
cluster <- t(cnv[,26:55]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
cnv <- cnv[,26:55] %>>% (.[,cluster$order]) %>>% (cbind(cnv[,1:25],.,cnv[,56:70]))
cluster <- t(cnv[,56:70]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
cnv <- cnv[,56:70] %>>% (.[,cluster$order]) %>>% (cbind(cnv[,1:55],.))
pd_col <- data.frame(Sample=c(rep("17cm-pre",20),rep("17cm-post",5),rep("23cm-pre",30),rep("23cm-post",15)))
rownames(pd_col) <- colnames(cnv)
ann_color <- list(Sample=c("17cm-pre"=blues9[2],
                           "17cm-post"=blues9[4],
                           "23cm-pre"=blues9[6],
                           "23cm-post"=blues9[8]))
# For profile
pheatmap(cnv,cluster_cols=F,clustering_distance_rows="binary",
         clustering_method="ward.D2",show_rownames=F,annotation_col=pd_col,
         annotation_colors=ann_color,
         color=c(colorRampPalette(c("cornflowerblue","gray95"))(127),
                 "white",
                 colorRampPalette(c("gray95","magenta"))(82)))
# For details
pheatmap(cnv,cluster_cols=F,clustering_distance_rows="binary",
         clustering_method="ward.D2",annotation_col=pd_col,
         color=c(colorRampPalette(c("cornflowerblue","gray95"))(127),
                 "white",
                 colorRampPalette(c("gray95","magenta"))(82)),
         annotation_colors=ann_color,cellwidth=5,cellheight=4,fontsize_row=4,fontsize_col=4,fontsize=6)
# Final version
temp <- cnv[,26:70] 
temp <- apply(temp,1,function(x) length(x[x != 0])) %>>% (temp[ . != 0,])

### CNA table for single-cell 23cm
write.csv(temp,file="~/Documents/cohort/ESCC/for_NEJM/T2. Single-cell_CNA.csv",quote=F)
###

temp <- apply(temp,1,function(x) length(x[x != 0])) %>>% (temp[ . <= 40,])
# Manually clustering, 1:632 rows showed amplifications for radio-resistance.
cluster <- dist(temp,method="binary") %>>% hclust(method="ward.D2")
temp <- temp[rev(cluster$order),]
p <- apply(temp,1,function(x) {
  a <- x[1:30]
  b <- x[31:45]
  a <- length(a[a != 0])
  b <- length(b[b != 0])
  fisher.test(matrix(c(a,b,30-a,15-b),nrow=2),alternative="less")$p.val
})
row_labels <- paste0(names(p[p < 0.05])," **")
row_subsets <- match(names(p[p < 0.05]),rownames(temp))
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",30),rep("Therapy_post",15)),
                        Sub_groups=anno_points(rep(1,45),
                                               pch=c(rep(9,14),rep(10,16),rep(9,3),rep(10,12)),
                                               size=unit(6,"pt"),border=F,ylim=c(0.5,1.5),
                                               gp=gpar(col=c(rep("indianred",14),rep("gold",16),
                                                             rep("indianred",3),rep("gold",12)))),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy_post"="#00DAE0")),
                        # annotation_height=c(0.5,0.5),
                        show_annotation_name=T,
                        # gap=unit(3,"pt"),
                        annotation_name_gp=gpar(fontsize=12,fontface="bold",fontfamily="serif"),
                        annotation_legend_param=list(title_gp=gpar(fontface="bold",fontfamily="serif",fontsize=26),
                                                     labels_gp=gpar(fontfamily="serif",fontsize=22),
                                                     grid_height=unit(10,"mm"),
                                                     grid_width=unit(10,"mm")))
row_ha <- rowAnnotation(link=row_anno_link(row_subsets,labels=row_labels,
                                           labels_gp=gpar(fontfamily="serif")),
                        width=unit(10,"pt") + max_text_width(row_labels))
ht <- Heatmap(temp,col=colorRamp2(c(-2,0,2),c("cornflowerblue","gray98","violet"),space="sRGB"),show_row_names=F,
              show_column_names=F,row_names_gp=gpar(fontsize=3),cluster_columns=F,cluster_rows=F,
              top_annotation=ha,top_annotation_height=unit(20,"pt"),name="CNA",
              heatmap_legend_param=list(color_bar="discrete",labels=c("-2","-1","0","1",">= 2"),
                                        title_gp=gpar(fontface="bold",fontfamily="serif",fontsize=26),
                                        labels_gp=gpar(fontfamily="serif",fontsize=22),
                                        grid_height=unit(10,"mm"),grid_width=unit(10,"mm"))) +
  row_ha
lgd <- Legend(c("A","B"),title="Sub_groups",type="points",pch=9:10,background="white",size=unit(6,"mm"),
              grid_height=unit(10,"mm"),grid_width=unit(10,"mm"),
              legend_gp=gpar(col=c("indianred","gold")),labels_gp=gpar(fontfamily="serif",fontsize=22),
              title_gp=gpar(fontface="bold",fontfamily="serif",fontsize=26))
draw(ht,annotation_legend_list=list(lgd))
decorate_heatmap_body("CNA",{
  i <- c(14,30,33)
  x <- i/ncol(temp)
  grid.lines(c(x[1],x[1]),c(0,1),gp=gpar(lwd=2,col="red"))
  grid.lines(c(x[2],x[2]),c(0,1),gp=gpar(lwd=2,col="red"))
  grid.lines(c(x[3],x[3]),c(0,1),gp=gpar(lwd=2,col="red"))
  # Top 632 (1910 - 632 = 1278)
  grid.lines(c(0,1),c(1278/nrow(temp),1278/nrow(temp)),gp=gpar(lwd=2,col="blue"))
  # Bottom 692
  grid.lines(c(0,1),c(692/nrow(temp),692/nrow(temp)),gp=gpar(lwd=2,col="green"))
})
# Genes of important CNV
cyto <- rownames(temp)[c(1:632,1219:1910)] %>>% (sub("^[+|-](.*)","\\1",.))
genes_cnv_23_single_cell_important <- lapply(cyto,function(x) which(cnv_23$an23.cytoBand == x)) %>>% unlist %>>% 
  unique %>>% (cnv_23[.,]) %>>% (paste(.$an23.Genes,collapse=",")) %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique

### Genes covered by important CNA
temp <- cnv[,26:70] 
temp <- apply(temp,1,function(x) length(x[x != 0])) %>>% (temp[ . != 0,])
temp <- sub("^[+|-](.*)","\\1",rownames(temp)) %>>% 
  lapply(function(x) which(cnv_23$an23.cytoBand == x)) %>>% unlist %>>% unique %>>% (cnv_23[.,]) %>>% 
  (paste(.$an23.Genes,collapse=",")) %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique %>>% 
  (bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID) %>>% enrichKEGG %>>% 
  as.data.frame
temp$geneID <- sapply(temp$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
                        (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>%
                        paste(collapse="/"),simplify=T)
write.csv(temp,file="~/Documents/cohort/ESCC/for_NEJM/T2. Single-cell_CNA_genes_KEGG.csv",quote=F)

lapply(cyto[1:632],function(x) which(cnv_23$an23.cytoBand == x)) %>>% unlist %>>% unique %>>% 
  (cnv_23[.,]) %>>% (paste(.$an23.Genes,collapse=",")) %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique %>>%
  writeLines("~/Documents/cohort/ESCC/for_NEJM/T3. Single-cell_important_CNA_genes_1.txt")
lapply(cyto[633:1324],function(x) which(cnv_23$an23.cytoBand == x)) %>>% unlist %>>% unique %>>% 
  (cnv_23[.,]) %>>% (paste(.$an23.Genes,collapse=",")) %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique %>>%
  writeLines("~/Documents/cohort/ESCC/for_NEJM/T3. Single-cell_important_CNA_genes_2.txt")
temp <- lapply(cyto,function(x) which(cnv_23$an23.cytoBand == x)) %>>% unlist %>>% unique %>>% 
  (cnv_23[.,]) %>>% (paste(.$an23.Genes,collapse=",")) %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique %>>%
  (bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID) %>>% enrichKEGG %>>% 
  as.data.frame
temp$geneID <- sapply(temp$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
                        (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>%
                        paste(collapse="/"),simplify=T)
write.csv(temp,file="~/Documents/cohort/ESCC/for_NEJM/T3. Single-cell_important_CNA_genes_KEGG.csv",quote=F)
temp <- lapply(cyto[1:632],function(x) which(cnv_23$an23.cytoBand == x)) %>>% unlist %>>% unique %>>% 
  (cnv_23[.,]) %>>% (paste(.$an23.Genes,collapse=",")) %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique %>>%
  (bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID) %>>% enrichKEGG %>>% 
  as.data.frame
temp$geneID <- sapply(temp$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
                        (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>%
                        paste(collapse="/"),simplify=T)
write.csv(temp,file="~/Documents/cohort/ESCC/for_NEJM/T3. Single-cell_important_CNA_genes_1_KEGG.csv",quote=F)
temp <- lapply(cyto[633:1324],function(x) which(cnv_23$an23.cytoBand == x)) %>>% unlist %>>% unique %>>% 
  (cnv_23[.,]) %>>% (paste(.$an23.Genes,collapse=",")) %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique %>>%
  (bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID) %>>% enrichKEGG %>>% 
  as.data.frame
temp$geneID <- sapply(temp$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
                        (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>%
                        paste(collapse="/"),simplify=T)
write.csv(temp,file="~/Documents/cohort/ESCC/for_NEJM/T3. Single-cell_important_CNA_genes_2_KEGG.csv",quote=F)
###

genes_snv_23_single_cell_important <- c("CFAP46","OBSCN","MUC6","OR8U1","OR8U8","HLA-DRB1","TRIM64C","CCDC8","ZNF717",
                                        "NPIPB6","MUC3A","PDE4DIP","MADCAM1","CDC27","MYO7A","AHNAK2","DNAH9","ODF4",
                                        "ZNF79","AP1G1","EXD3")
intersect(genes_cnv_23_single_cell_important,genes_snv_23_single_cell_important)
# Genes of higher frequency
genes_cnv_23_single_cell_higher <- sapply(row_labels,str_split," ") %>>% sapply(function(x) x[1]) %>>% str_remove("^[+|-]") %>>% 
  (cnv_23[(as.character(cnv_23$an23.cytoBand) %in% .),]) %>>% 
  (.$an23.Genes) %>>% paste(collapse=",") %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unlist %>>% unique
rm(temp)
# ------

# Heatmap of MUC6, PTEN, (881), (884) for single-cell
# ------
snv_muc6_pten <- snv[c("MUC6","PTEN"),]
cnv_881_884 <- cbind(cnv_17[c(10006,10049:10051),1:28],cnv_23[c(10006,10049:10051),c(4:48,52:53)])
temp <- cnv_881_884[,4:73]
temp <- temp[,match(colnames(snv_muc6_pten),colnames(temp))]
temp <- apply(temp,2,function(x) sapply(x,function(x) {
  if (x > 2) return(1)
  if (x < 2) return(-1)
  if (x == 2) return(0)
}))
temp <- cbind(cnv_881_884[,1:3],temp,cnv_881_884[,74:75])
cnv_881_884 <- matrix(0,nrow=0,ncol=70)
temp <- sapply(unique(as.character(temp$an23.cytoBand)),function(i) {
  x <- filter(temp,an23.cytoBand == i)
  y <- apply(x[,4:73],2,function(a) {
    if (length(a[a == 1]) > 0)
      j <- length(a[a == 1])
    else
      j <- 0
    if (length(a[a == -1]) > 0)
      k <- -length(a[a == -1])
    else
      k <- 0
    return(c(j,k))
  })
  rownames(y) <- c(paste0("+",i),paste0("-",i))
  colnames(y) <- colnames(temp)[4:73]
  cnv_881_884 <<- rbind(cnv_881_884,y)
})
rm(temp)
mutation <- rbind(snv_muc6_pten,cnv_881_884)
cluster <- t(mutation[,1:20]) %>>% dist(method="canberra") %>>% hclust(method="ward.D")
mutation <- mutation[,1:20] %>>% (.[,cluster$order]) %>>% cbind(mutation[,21:70])
cluster <- t(mutation[,21:25]) %>>% dist(method="canberra") %>>% hclust(method="ward.D")
mutation <- mutation[,21:25] %>>% (.[,cluster$order]) %>>% (cbind(mutation[,1:20],.,mutation[,26:70]))
cluster <- t(mutation[,26:55]) %>>% dist(method="canberra") %>>% hclust(method="ward.D")
mutation <- mutation[,26:55] %>>% (.[,cluster$order]) %>>% (cbind(mutation[,1:25],.,mutation[,56:70]))
cluster <- t(mutation[,56:70]) %>>% dist(method="canberra") %>>% hclust(method="ward.D")
mutation <- mutation[,56:70] %>>% (.[,cluster$order]) %>>% (cbind(mutation[,1:55],.))
pheatmap(mutation,cluster_cols=F,clustering_distance_rows="canberra",
         clustering_method="ward.D2",
         annotation_col=pd_col,
         color=c(colorRampPalette(c("blue","cornflowerblue"))(2),
                 "white",
                 colorRampPalette(c("pink","red"))(4)),
         annotation_colors=ann_color,cellwidth=5,cellheight=4,fontsize_row=4,fontsize_col=4,fontsize=6)
# Final version
# Removing 17cm and 17cm_post
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",30),rep("Therapy_post",15)),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy_post"="#00DAE0")),
                        annotation_height=0.5,
                        show_annotation_name=T,
                        annotation_name_gp=gpar(fontsize=11))
Heatmap(mutation[,26:70],col=colorRamp2(c(-3,0,4),c("cornflowerblue","gray98","violet"),space="sRGB"),
        show_row_names=T,show_column_names=F,row_names_gp=gpar(fontsize=11),cluster_columns=F,
        clustering_distance_rows="canberra",clustering_method_rows="ward.D2",row_dend_reorder=F,
        show_row_dend=F,top_annotation=ha,top_annotation_height=unit(10,"pt"),name="Mutation",
        rect_gp=gpar(col="gray98"),heatmap_legend_param=list(color_bar="continuous"))
# ------

# Heatmap of MUC6 and PTEN, (881), (884) for single-cell with modification
# ------
snv_23cm_muc6_and_pten <- apply(snv_muc6_pten[,26:70],2,function(x) 
  if (length(x[x != 0]) == 2) return(1) else return(0)) %>>% t %>>% as.data.frame
rownames(snv_23cm_muc6_and_pten) <- "MUC6 and PTEN"
cnv_23cm_881_884 <- apply(cnv_881_884[,26:70],2,function(x) {
  if (x[1] == 0 & x[2] == 0)
    a <- 0
  if (x[1] == 0 & x[2] != 0)
    a <- x[2]
  if (x[1] != 0 & x[2] == 0)
    a <- x[1]
  if (x[3] == 0 & x[4] == 0)
    b <- 0
  if (x[3] == 0 & x[4] != 0)
    b <- x[4]
  if (x[3] != 0 & x[4] == 0)
    b <- x[3]
  return(c(a,b))
})
rownames(cnv_23cm_881_884) <- c("19q13.13","19q13.31")
mutation_23cm <- rbind(snv_23cm_muc6_and_pten,cnv_23cm_881_884)
cluster <- t(mutation_23cm[,1:30]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_23cm <- mutation_23cm[,1:30] %>>% (.[,cluster$order]) %>>% cbind(mutation_23cm[,31:45])
cluster <- t(mutation_23cm[,31:45]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_23cm <- mutation_23cm[,31:45] %>>% (.[,cluster$order]) %>>% (cbind(mutation_23cm[,1:30],.))
pd_col_23cm <- as.data.frame(pd_col[26:70,])
pd_col_23cm[,1] <- as.character(pd_col_23cm[,1])
colnames(pd_col_23cm) <- "Sample"
rownames(pd_col_23cm) <- colnames(mutation_23cm)
pheatmap(mutation_23cm,cluster_cols=F,clustering_distance_rows="binary",
         clustering_method="ward.D2",
         annotation_col=pd_col_23cm,
         color=c(colorRampPalette(c("blue","cornflowerblue"))(2),
                 "white",
                 colorRampPalette(c("pink","red"))(2)),
         annotation_colors=list(Sample=c("23cm-pre"=blues9[6],"23cm-post"=blues9[8])),
         cellwidth=5,cellheight=4,fontsize_row=4,fontsize_col=4,fontsize=6)
# New...
mutation_23cm_combine <- rbind(snv_23cm_muc6_and_pten,cnv_881_884[,26:70])
cluster <- t(mutation_23cm_combine[,1:30]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_23cm_combine <- mutation_23cm_combine[,1:30] %>>% (.[,cluster$order]) %>>% cbind(mutation_23cm_combine[,31:45])
cluster <- t(mutation_23cm_combine[,31:45]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_23cm_combine <- mutation_23cm_combine[,31:45] %>>% (.[,cluster$order]) %>>% (cbind(mutation_23cm_combine[,1:30],.))
mutation_23cm_combine[c(3,5),] <- adply(mutation_23cm_combine[c(3,5),],1,function(x) {
  x[x != 0] <- 1
  return(x)
})
mutation_23cm_combine[1,] <- sapply(mutation_23cm_combine[1,],function(x) {
  x[x != 0] <- 2
  return(x)
})
mutation_23cm_combine[c(2,4),] <- adply(mutation_23cm_combine[c(2,4),],1,function(x) {
  x[x != 0] <- 3
  return(x)
})
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",30),rep("Therapy_post",15)),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy_post"="#00DAE0")),
                        annotation_height=0.5,
                        show_annotation_name=T,
                        annotation_name_gp=gpar(fontsize=11,fontface="bold",fontfamily="serif"),
                        annotation_legend_param=list(title_gp=gpar(fontface="bold",fontfamily="serif"),
                                                     labels_gp=gpar(fontfamily="serif")))
Heatmap(mutation_23cm_combine,
        col=colorRamp2(c(0,1,2,3),c("gray98","cornflowerblue","yellowgreen","violetred"),space="sRGB"),
        show_row_names=T,show_column_names=F,row_names_gp=gpar(fontsize=11,fontfamily="serif"),cluster_columns=F,
        clustering_distance_rows="binary",clustering_method_rows="ward.D2",row_dend_reorder=F,
        show_row_dend=F,top_annotation=ha,top_annotation_height=unit(10,"pt"),name="Mutation",
        rect_gp=gpar(col="gray98"),heatmap_legend_param=list(color_bar="discrete",
                                                             labels=c("","Deletion","SNV","Amplification"),
                                                             labels_gp=gpar(fontfamily="serif"),
                                                             title_gp=gpar(fontface="bold",fontfamily="serif")))
# ------

# SNV features of single-cell by CNV profile
# ------
groups_23cm_pre <- data.frame(Sample=colnames(cnv)[26:55],Group=c(rep("Group1",14),rep("Group2",16)))
groups_23cm_post <- data.frame(Sample=colnames(cnv)[56:70],Group=c(rep("Group1",3),rep("Group2",12)))
p_less_23cm_pre <- cbind(snv[,filter(groups_23cm_pre,Group == "Group1")[,1]],
                         snv[,filter(groups_23cm_pre,Group == "Group2")[,1]]) %>>% 
  apply(1,function(x) {
    a <- x[1:14]
    b <- x[15:30]
    c <- matrix(c(length(a[a != 0]),length(b[b != 0]),length(a[a == 0]),length(b[b == 0])),nrow=2) %>>% 
      fisher.test(alternative="less")
    return(c$p.val)
  })
p_greater_23cm_pre <- cbind(snv[,filter(groups_23cm_pre,Group == "Group1")[,1]],
                            snv[,filter(groups_23cm_pre,Group == "Group2")[,1]]) %>>% 
  apply(1,function(x) {
    a <- x[1:14]
    b <- x[15:30]
    c <- matrix(c(length(a[a != 0]),length(b[b != 0]),length(a[a == 0]),length(b[b == 0])),nrow=2) %>>% 
      fisher.test(alternative="greater")
    return(c$p.val)
  })
p_greater_23cm_post <- cbind(snv[,filter(groups_23cm_post,Group == "Group1")[,1]],
                             snv[,filter(groups_23cm_post,Group == "Group2")[,1]]) %>>% 
  apply(1,function(x) {
    a <- x[1:3]
    b <- x[4:15]
    c <- matrix(c(length(a[a != 0]),length(b[b != 0]),length(a[a == 0]),length(b[b == 0])),nrow=2) %>>% 
      fisher.test(alternative="greater")
    return(c$p.val)
  })
# ------

# CNV features and list of single-cell by CNV profile
# ------
p_23cm_pre_cnv <- cbind(cnv[,filter(groups_23cm_pre,Group == "Group1")[,1]],
                        cnv[,filter(groups_23cm_pre,Group == "Group2")[,1]]) %>>% 
  apply(1,function(x) {
    a <- x[1:14]
    b <- x[15:30]
    c <- wilcox.test(a,b)
    return(c$p.val)
  })
diff_23cm_pre_cnv <- p_23cm_pre_cnv[!is.na(p_23cm_pre_cnv) & p_23cm_pre_cnv < 0.05] %>>% names %>>% 
  (sub("[+|-](.*)","\\1",.)) %>>% (cnv_23[match(.,as.character(cnv_23$an23.cytoBand)),])
write.csv(diff_23cm_pre_cnv,file="~/Documents/cohort/ESCC/for_NEJM/diff_23cm_pre_cnv.csv")
library(clusterProfiler)
kegg_23cm_pre_cnv <- diff_23cm_pre_cnv$an23.Genes %>>% as.character %>>% paste(collapse=",") %>>% 
  (strsplit(.,",")[[1]]) %>>% unique %>>% 
  (bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID) %>>% 
  enrichKEGG(pvalueCutoff=1,minGSSize=1,qvalueCutoff=1) %>>% as.data.frame
gene <- diff_23cm_pre_cnv$an23.Genes %>>% as.character %>>% paste(collapse=",") %>>% (strsplit(.,",")[[1]]) %>>% 
  unique %>>% bitr(fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
kegg_23cm_pre_cnv$geneID <- apply(kegg_23cm_pre_cnv,1,function(x) 
  paste(gene[match(strsplit(x[8],"/")[[1]],gene$ENTREZID),"SYMBOL"],collapse="/"))
write.csv(kegg_23cm_pre_cnv,file="~/Documents/cohort/ESCC/for_NEJM/diff_23cm_pre_cnv_gene_kegg.csv")
# ------

# ============================================================================

# Comebine SNV results and CNV results of ESCC cohort
# ======

# Reading patients information
# ------
pi <- read.csv("~/Documents/cohort/ESCC/final_need_fixed_20160808.csv")
# ------

# Reading SNV results (Somatic)
# ------
snv_cohort <- read.table("~/Documents/cohort/ESCC/for_NEJM/muts_matrix_20170118_cui_from_g.vcf_remove_5_neighbours.tsv",
                         sep="\t",header=T,row.names=1)
temp <- snv_cohort[,c("R1509285859","R1509285860")] %>>% apply(1,function(x) {
  if (length(x[x == "."]) == 1)
    if (x[1] == ".")
      a <- x[2]
    else
      a <- x[1]
  else
    if (length(x[x == "."]) == 2)
      a <- "."
    else
      if (x[1] == "0/0")
        a <- x[2]
      else
        a <- x[1]
  return(a)
})
which(colnames(snv_cohort) == "R1509285859")
snv_cohort[,78] <- temp
snv_cohort <- snv_cohort[,-79]
snv_cohort <- snv_cohort[,-155]
rm(temp)
# Transform SNV from sites to genes
snv_cohort_gene <- sapply(levels(snv_cohort$symbol),function(x) {
  a <- filter(snv_cohort,symbol == x)
  b <- apply(a[,1:153],2,function(y) length(y[y != "." & y != "0/0"]))
  return(b)
}) %>>% t
temp_group <- c(as.character(filter(pi,level == "primary")$病人姓名),
                as.character(filter(pi,level == "therapy")$病人姓名),
                as.character(filter(pi,level == "therapy_after")$病人姓名),
                as.character(filter(pi,level == "relapse")$病人姓名))
temp_group <- temp_group[temp_group != "R1509285860"]
snv_cohort <- snv_cohort[,match(temp_group,colnames(snv_cohort))]
snv_cohort_gene <- as.data.frame(snv_cohort_gene[,match(temp_group,colnames(snv_cohort_gene))])

### Table of cohort snv (somatic)
write.csv(snv_cohort,file="~/Documents/cohort/ESCC/for_NEJM/T4. Cohort_snv_variants.csv",quote=F)
write.csv(snv_cohort_gene,file="~/Documents/cohort/ESCC/for_NEJM/T4. Cohort_snv_genes.csv",quote=F)
###

rm(temp_group)
# ------

# Reading CNV results
# ------
cnv_cohort <- read.table("~/Documents/cohort/ESCC/cnv/cnvr_table.tsv",sep="\t",header=T,row.names=1)
cnv_cohort <- apply(cnv_cohort[,1:3],1,function(x) paste(x,collapse="_")) %>>% 
  sapply(function(x) gsub(" +","",x),simplify=T) %>>% (cbind(cnv_cohort,.))
cnv_cohort_region <- cnv_cohort[,7:160]
rownames(cnv_cohort_region) <- cnv_cohort[,161]
cytoBand <- read.table("~/Documents/cohort/ESCC/cnv/cnvr_annovar.hg19_cytoBand",sep="\t",header=F)
cytoBand <- apply(cytoBand[,3:5],1,function(x) paste0("chr",paste(x,collapse="_"))) %>>% 
  sapply(function(x) gsub(" +","",x),simplify=T) %>>% (cbind(cytoBand,.))
cnv_cohort <- cbind(cnv_cohort,cytoBand[match(cnv_cohort[,161],cytoBand[,9]),c(2,9)])
cnv_cohort <- cnv_cohort[!is.na(cnv_cohort[,162]),c(1:160,162)]
# Transforming CNV data frame
cnv_cohort <- apply(cnv_cohort[,7:160],2,function(x) sapply(x,function(y) as.integer(strsplit(y,"N")[[1]][2]))) %>>% 
  (cbind(cnv_cohort[,1:6],.,cnv_cohort[,161]))
cnv_cohort_region <- apply(cnv_cohort_region,2,function(x) sapply(x,function(y) as.integer(strsplit(y,"N")[[1]][2])))
# Table of region to table of cytoband
temp <- cnv_cohort[,7:160]
temp_group <- c(as.character(filter(pi,level == "primary")$病人姓名),
                as.character(filter(pi,level == "therapy")$病人姓名),
                as.character(filter(pi,level == "therapy_after")$病人姓名),
                as.character(filter(pi,level == "relapse")$病人姓名))
temp <- temp[,match(temp_group,sub("(.*).KY0(.*)","\\1",colnames(temp)))]
temp <- apply(temp,2,function(x) sapply(x,function(x) {
  if (x > 2) return(1)
  if (x < 2) return(-1)
  if (x == 2) return(0)
}))
temp <- cbind(cnv_cohort[,1:6],temp,cnv_cohort[,161])
cnv_cohort <- matrix(0,nrow=0,ncol=154)
temp <- sapply(levels(temp[,161]),function(i) {
  x <- filter(temp,temp[,161] == i)
  y <- apply(x[,7:160],2,function(a) {
    if (length(a[a == 1]) > 0)
      j <- length(a[a == 1])
    else
      j <- 0
    if (length(a[a == -1]) > 0)
      k <- -length(a[a == -1])
    else
      k <- 0
    return(c(j,k))
  })
  rownames(y) <- c(paste0("+",i),paste0("-",i))
  colnames(y) <- colnames(temp)[7:160]
  cnv_cohort <<- rbind(cnv_cohort,y)
})
rm(temp)
rm(temp_group)
colnames(cnv_cohort) <- sub("(.*).KY0(.*)","\\1",colnames(cnv_cohort))
colnames(cnv_cohort_region) <- sub("(.*).KY0(.*)","\\1",colnames(cnv_cohort_region))
cnv_cohort_region <- cnv_cohort_region[,match(colnames(cnv_cohort),colnames(cnv_cohort_region))]
temp <- cnv_cohort[,65:66] %>>% apply(1,function(x) {
  if (x[1] == x[2])
    a <- x[1]
  else
  {
    a <- mean(x)
    set.seed(1984)
    if (rbinom(1,size=1,prob=0.5) == 0)
      a <- floor(a)
    else
      a <- ceiling(a)
  }
})
cnv_cohort <- cbind(cnv_cohort[,1:64],R1509285859=temp,cnv_cohort[,67:154])
temp <- cnv_cohort_region[,65:66] %>>% apply(1,function(x) {
  if (x[1] == x[2])
    a <- x[1]
  else
  {
    a <- mean(x)
    set.seed(1984)
    if (rbinom(1,size=1,prob=0.5) == 0)
      a <- floor(a)
    else
      a <- ceiling(a)
  }
})
cnv_cohort_region <- cbind(cnv_cohort_region[,1:64],R1509285859=temp,cnv_cohort_region[,67:154])
rm(temp)
cnv_cohort <- as.data.frame(cnv_cohort)
cnv_cohort_region <- as.data.frame(cnv_cohort_region)

### Table of cohort CNA
write.csv(cnv_cohort,file="~/Documents/cohort/ESCC/for_NEJM/T5. Cohort_cna_cytobands.csv",quote=F)
write.csv(cnv_cohort_region,file="~/Documents/cohort/ESCC/for_NEJM/T5. Cohort_cna_regions.csv",quote=F)
###

# ------

# Heatmap of SNV
# ------
library(pheatmap)
cluster <- t(snv_cohort_gene[,1:95]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
snv_cohort_gene <- snv_cohort_gene[,1:95] %>>% (.[,cluster$order]) %>>% cbind(snv_cohort_gene[,96:153])
cluster <- t(snv_cohort_gene[,96:126]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
snv_cohort_gene <- snv_cohort_gene[,96:126] %>>% (.[,cluster$order]) %>>%
  (cbind(snv_cohort_gene[,1:95],.,snv_cohort_gene[,127:153]))
cluster <- t(snv_cohort_gene[,127:149]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
snv_cohort_gene <- snv_cohort_gene[,127:149] %>>% (.[,cluster$order]) %>>% 
  (cbind(snv_cohort_gene[,1:126],.,snv_cohort_gene[,150:153]))
cluster <- t(snv_cohort_gene[,150:153]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
snv_cohort_gene <- snv_cohort_gene[,150:153] %>>% (.[,cluster$order]) %>>% (cbind(snv_cohort_gene[,1:149],.))
pd <- data.frame(Samples=c(rep("Primary",95),rep("Therapy",31),rep("Therapy_post",23),rep("Relapse",4)))
rownames(pd) <- colnames(snv_cohort_gene)
ann_color <- list(Samples=c("Primary"=blues9[2],
                           "Therapy"=blues9[4],
                           "Therapy_post"=blues9[6],
                           "Relapse"=blues9[8]))
pheatmap(snv_cohort_gene,cluster_cols=F,clustering_distance_rows="maximum",treeheight_row=0,
         clustering_method="ward.D2",show_rownames=F,annotation_col=pd,legend_breaks=c(0,1,2,3),
         color=c("white",colorRampPalette(c("cornflowerblue","magenta"))(3)),annotation_colors=ann_color)
# Final version
# Different clusters, different colors
temp <- snv_cohort_gene
temp[,1:27] <- apply(temp[,1:27],2,function(x) {
  x[x != 0] <- 1
  return(x)
})
temp[,28:105] <- apply(temp[,28:105],2,function(x) {
  x[x != 0] <- 2
  return(x)
})
temp[,106:107] <- apply(temp[,106:107],2,function(x) {
  x[x != 0] <- 1
  return(x)
})
temp[,108:129] <- apply(temp[,108:129],2,function(x) {
  x[x != 0] <- 2
  return(x)
})
temp[,130:132] <- apply(temp[,130:132],2,function(x) {
  x[x != 0] <- 1
  return(x)
})
temp[,133:149] <- apply(temp[,133:149],2,function(x) {
  x[x != 0] <- 2
  return(x)
})
a <- temp[,150]
a[a != 0] <- 1
temp[,150] <- a
rm(a)
temp[,151:153] <- apply(temp[,151:153],2,function(x) {
  x[x != 0] <- 2
  return(x)
})
cluster <- dist(temp,method="maximum") %>>% hclust(method="ward.D2")
temp <- temp[rev(cluster$order),]
p <- apply(temp,1,function(x) {
  a <- x[1:95]
  b <- x[96:149]
  a <- length(a[a != 0])
  b <- length(b[b != 0])
  fisher.test(matrix(c(a,b,95-a,54-b),nrow=2),alternative="less")$p.val
})
row_labels <- paste0(names(p[p < 0.05])," **")
row_subsets <- match(names(p[p < 0.05]),rownames(temp))
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",95),rep("Therapy",31),rep("Therapy_post",23),rep("Relapse",4)),
                        Sub_groups=anno_points(rep(1,153),
                                               pch=c(rep(7,27),rep(8,68),rep(8,10),rep(7,2),rep(8,19),
                                                     rep(8,3),rep(7,3),rep(8,17),rep(7,1),rep(8,3)),
                                               size=unit(6,"pt"),border=F,ylim=c(0.5,1.5),
                                               gp=gpar(col=c(rep("palevioletred",27),rep("cornflowerblue",68),
                                                             rep("cornflowerblue",10),rep("palevioletred",2),
                                                             rep("cornflowerblue",19),rep("cornflowerblue",3),
                                                             rep("palevioletred",3),rep("cornflowerblue",17),
                                                             rep("palevioletred",1),rep("cornflowerblue",3)))),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy"="#E199FF",
                                               "Therapy_post"="#00DAE0",
                                               "Relapse"="#FF9289")),
                        # annotation_height=c(0.5,0.5),
                        show_annotation_name=T,
                        # gap=unit(3,"pt"),
                        annotation_name_gp=gpar(fontsize=12,fontface="bold",fontfamily="serif"),
                        annotation_legend_param=list(title_gp=gpar(fontface="bold",fontfamily="serif",fontsize=26),
                                                     labels_gp=gpar(fontfamily="serif",fontsize=22),
                                                     grid_height=unit(10,"mm"),
                                                     grid_width=unit(10,"mm")))
row_ha <- rowAnnotation(link=row_anno_link(row_subsets,labels=row_labels,
                                           labels_gp=gpar(fontfamily="serif",fontsize=22)),
                        width=unit(10,"pt") + max_text_width(row_labels))
ht <- Heatmap(temp,col=colorRamp2(c(0,1,2),c("gray98","palevioletred","cornflowerblue"),space="sRGB"),
              show_row_names=F,show_column_names=F,row_names_gp=gpar(fontsize=1),cluster_columns=F,
              show_row_dend=F,top_annotation=ha,top_annotation_height=unit(20,"pt"),name="SNV",
              show_heatmap_legend=F) + row_ha
lgd <- Legend(c("A","B"),title="Sub_groups",type="points",pch=7:8,background="white",size=unit(6,"mm"),
              grid_height=unit(10,"mm"),grid_width=unit(10,"mm"),
              legend_gp=gpar(col=c("palevioletred","cornflowerblue")),
              title_gp=gpar(fontface="bold",fontfamily="serif",fontsize=26),
              labels_gp=gpar(fontfamily="serif",fontsize=22))
# 8 * 16 pdf
draw(ht,annotation_legend_list=list(lgd))
decorate_heatmap_body("SNV",{
  i <- c(27,95,105,107,126,129,132,149,150)
  x <- i/ncol(temp)
  grid.lines(c(x[1],x[1]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[2],x[2]),c(0,1),gp=gpar(lwd=2,col="green"))
  grid.lines(c(x[3],x[3]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[4],x[4]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[5],x[5]),c(0,1),gp=gpar(lwd=2,col="green"))
  grid.lines(c(x[6],x[6]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[7],x[7]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[8],x[8]),c(0,1),gp=gpar(lwd=2,col="green"))
  grid.lines(c(x[9],x[9]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
})
rm(temp)
# ------

# Heatmap of CNV
# ------
cluster <- t(cnv_cohort[,1:95]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
cnv_cohort <- cnv_cohort[,1:95] %>>% (.[,cluster$order]) %>>% cbind(cnv_cohort[,96:153])
cluster <- t(cnv_cohort[,96:126]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
cnv_cohort <- cnv_cohort[,96:126] %>>% (.[,cluster$order]) %>>%
  (cbind(cnv_cohort[,1:95],.,cnv_cohort[,127:153]))
cluster <- t(cnv_cohort[,127:149]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
cnv_cohort <- cnv_cohort[,127:149] %>>% (.[,rev(cluster$order)]) %>>% 
  (cbind(cnv_cohort[,1:126],.,cnv_cohort[,150:153]))
cluster <- t(cnv_cohort[,150:153]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
cnv_cohort <- cnv_cohort[,150:153] %>>% (.[,cluster$order]) %>>% (cbind(cnv_cohort[,1:149],.))
rownames(pd) <- colnames(cnv_cohort)
# Removing empty rows
temp <- cnv_cohort[apply(cnv_cohort,1,function(x) length(x[x != 0])) %>>% (. != 0),]
pheatmap(temp,cluster_cols=F,clustering_distance_rows="euclidean",treeheight_row=0,
         clustering_method="ward.D",show_rownames=F,annotation_col=pd,
         legend_breaks=c(-10,0,11),
         color=c(colorRampPalette(c("cornflowerblue","#AFC7F3"))(9),
                 "gray98",colorRampPalette(c("#FCA98A","orangered"))(10)),
         annotation_colors=ann_color)
# Final version
# Manually clustering, 1:632 rows showed amplifications for radio-resistance.
# Different weights for members parameter can reorder the clustering.
cluster <- dist(temp,method="euclidean") %>>% hclust(method="ward.D2") %>>% 
  reorder(wts=apply(temp,1,function(x) length(x[x != 0])),agglo.FUN="uwmean")
temp <- temp[rev(cluster$order),]
# Show the top 135 CNAs
temp <- temp[1:135,]
p <- apply(temp,1,function(x) {
  a <- x[1:95]
  b <- x[96:149]
  a <- length(a[a != 0])
  b <- length(b[b != 0])
  fisher.test(matrix(c(a,b,95-a,54-b),nrow=2),alternative="less")$p.val
})
row_labels <- paste0(names(p[p < 0.05])," **")
row_subsets <- match(names(p[p < 0.05]),rownames(temp))
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",95),rep("Therapy",31),rep("Therapy_post",23),rep("Relapse",4)),
                        Sub_groups=anno_points(rep(1,153),
                                               pch=c(rep(7,18),rep(8,17),rep(9,12),rep(10,46),rep(11,2),
                                                     rep(10,9),rep(7,22),rep(10,10),rep(7,13),rep(7,4)),
                                               size=unit(5,"pt"),border=F,ylim=c(0.5,1.5),
                                               gp=gpar(col=c(rep("cornflowerblue",18),rep("lightpink",17),
                                                             rep("indianred",12),rep("gold",46),rep("forestgreen",2),
                                                             rep("gold",9),rep("cornflowerblue",22),
                                                             rep("gold",10),rep("cornflowerblue",13),
                                                             rep("cornflowerblue",4)))),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy"="#E199FF",
                                               "Therapy_post"="#00DAE0",
                                               "Relapse"="#FF9289")),
                        # annotation_height=c(0.5,0.5),
                        show_annotation_name=T,
                        # gap=unit(3,"pt"),
                        annotation_name_gp=gpar(fontsize=18,fontface="bold",fontfamily="serif"),
                        annotation_legend_param=list(title_gp=gpar(fontface="bold",fontfamily="serif",fontsize=26),
                                                     labels_gp=gpar(fontfamily="serif",fontsize=22),
                                                     grid_height=unit(10,"mm"),
                                                     grid_width=unit(10,"mm")))
row_ha <- rowAnnotation(link=row_anno_link(row_subsets,labels=row_labels,
                                           labels_gp=gpar(fontfamily="serif",fontsize=18)),
                        width=unit(2,"cm") + max_text_width(row_labels))
ht <- Heatmap(temp,col=colorRamp2(c(-2,0,2),c("cornflowerblue","gray98","violet"),space="sRGB"),show_row_names=F,
              show_column_names=F,row_names_gp=gpar(fontsize=3),cluster_columns=F,cluster_rows=F,
              show_row_dend=F,top_annotation=ha,top_annotation_height=unit(30,"pt"),name="CNA",
              # rect_gp=gpar(col="gray98"),
              heatmap_legend_param=list(color_bar="discrete",
                                        title_gp=gpar(fontface="bold",fontfamily="serif",fontsize=26),
                                        labels=c("-2","-1","0","1",">= 2"),
                                        grid_height=unit(10,"mm"),grid_width=unit(10,"mm"),
                                        labels_gp=gpar(fontfamily="serif",fontsize=22))) + row_ha
lgd <- Legend(c("A","B","C","D","E"),title="Sub_groups",type="points",pch=7:11,background="white",size=unit(6,"mm"),
              grid_height=unit(10,"mm"),grid_width=unit(10,"mm"),
              legend_gp=gpar(col=c("cornflowerblue","lightpink","indianred","gold","forestgreen")),
              title_gp=gpar(fontface="bold",fontfamily="serif",fontsize=26),
              labels_gp=gpar(fontfamily="serif",fontsize=22))
draw(ht,annotation_legend_list=list(lgd))
decorate_heatmap_body("CNA",{
  i <- c(18,35,47,93,95,104,126,136,149)
  x <- i/ncol(temp)
  grid.lines(c(x[1],x[1]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[2],x[2]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[3],x[3]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[4],x[4]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[5],x[5]),c(0,1),gp=gpar(lwd=2,col="green"))
  grid.lines(c(x[6],x[6]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[7],x[7]),c(0,1),gp=gpar(lwd=2,col="green"))
  grid.lines(c(x[8],x[8]),c(0,1),gp=gpar(lwd=1,lty=2,col="red"))
  grid.lines(c(x[9],x[9]),c(0,1),gp=gpar(lwd=2,col="green"))
  # Top 28 CNV
  grid.lines(c(0,1),c(107/nrow(temp),107/nrow(temp)),gp=gpar(lwd=1,lty=2,col="blue"))
})
# Changes of 12p12.1/16q21 in each patient
pi$level <- relevel(pi$level,"therapy_after")
pi$level <- relevel(pi$level,"therapy")
pi$level <- relevel(pi$level,"primary")
patient_12p12.1 <- sapply(unique(pi$姓名),function(x) {
  pipi <- pi[order(pi$level),]
  a <- as.character(pipi[which(pipi$姓名 == x),3])
  if (length(a) == 1) {
    b <- as.data.frame(temp[26,a]) # 12p12.1 in row 26 of temp
    colnames(b) <- pipi[which(pipi$姓名 == x),5]
    rownames(b) <- rownames(temp)[26]
    b <- pipi[which(pipi$姓名 == x)[1],2] %>>% 
      str_remove("[A|B|C|D]") %>>% (cbind(b,label=.))
    return(b)
  } else {
    b <- temp[26,a]
    colnames(b) <- pipi[which(pipi$姓名 == x),5]
    b <- pipi[which(pipi$姓名 == x)[1],2] %>>% 
      str_remove("[A|B|C|D]") %>>% (cbind(b,label=.))
    return(b)
  }
}) %>>% rbind.fill %>>% {. <- .[,c(4,1,5,2,3)]} # melt %>>% dcast(label~variable)
patient_16q21 <- sapply(unique(pi$姓名),function(x) {
  pipi <- pi[order(pi$level),]
  a <- as.character(pipi[which(pipi$姓名 == x),3])
  if (length(a) == 1) {
    b <- as.data.frame(temp[12,a]) # 16q21 in row 12 of temp
    colnames(b) <- pipi[which(pipi$姓名 == x),5]
    rownames(b) <- rownames(temp)[12]
    b <- pipi[which(pipi$姓名 == x)[1],2] %>>% 
      str_remove("[A|B|C|D]") %>>% (cbind(b,label=.))
    return(b)
  } else {
    b <- temp[12,a]
    colnames(b) <- pipi[which(pipi$姓名 == x),5]
    b <- pipi[which(pipi$姓名 == x)[1],2] %>>% 
      str_remove("[A|B|C|D]") %>>% (cbind(b,label=.))
    return(b)
  }
}) %>>% rbind.fill %>>% {. <- .[,c(4,1,5,2,3)]}
# Genes of p < 0.05 cohort CNVs
cyto_cohort_sig <- names(p[p < 0.05]) %>>% str_remove("\\+") %>>% 
  sapply(function(x) as.character(cytoBand[which(as.character(cytoBand[,2]) == x),9])) %>>% 
  unlist %>>% unique
paste(sub("chr(.*)_(.*)_(.*)","\\1",cyto_cohort_sig),
      sub("chr(.*)_(.*)_(.*)","\\2",cyto_cohort_sig),
      sub("chr(.*)_(.*)_(.*)","\\3",cyto_cohort_sig),
      "0","-",sep="\t") %>>% 
  writeLines("~/Documents/cohort/ESCC/for_NEJM/cyto_cnv_cohort_significant.avinput")
genes_cnv_cohort_sig <- 
  read.table("~/Documents/cohort/ESCC/for_NEJM/cyto_cnv_cohort_significant.avoutput.hg19_bed",
             header=F,sep="\t",as.is=T) %>>% (.[,2]) %>>% str_remove("Name=") %>>% 
  str_split(",") %>>% unlist %>>% unique
bitr(genes_cnv_cohort_sig,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>% 
  enrichKEGG %>>% as.data.frame # n2s
# Genes of 28 important cohort CNVs
cyto_cohort <- rownames(temp)[1:28] %>>% (sub("^[+|-](.*)","\\1",.))
cyto_cohort_regions <- sapply(cyto_cohort,function(x) 
  as.character(cytoBand[which(as.character(cytoBand[,2]) == x),9])) %>>% 
  unlist %>>% unique
paste(sub("chr(.*)_(.*)_(.*)","\\1",cyto_cohort_regions),
      sub("chr(.*)_(.*)_(.*)","\\2",cyto_cohort_regions),
      sub("chr(.*)_(.*)_(.*)","\\3",cyto_cohort_regions),
      "0","-",sep="\t") %>>% 
  writeLines("~/Documents/cohort/ESCC/for_NEJM/cyto_cnv_cohort_important.avinput")
genes_cnv_cohort_important <- 
  read.table("~/Documents/cohort/ESCC/for_NEJM/cyto_cnv_cohort_important.avoutput.variant_function",
             header=F,sep="\t",as.is=T) %>>% (.[grep("exonic",.[,1]),2]) %>>% paste(collapse=",") %>>% 
  (strsplit(.,",",fixed=T)[[1]]) %>>% unique

### Table of important genes in cohort cna
writeLines(genes_cnv_cohort_important,"~/Documents/cohort/ESCC/for_NEJM/T6. Cohort_important_cna_genes.txt")
temp <- bitr(genes_cnv_cohort_important,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>%
  enrichKEGG %>>% as.data.frame
temp$geneID <- sapply(temp$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
                        (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
                        paste(collapse="/"),simplify=T)
write.csv(temp,file="~/Documents/cohort/ESCC/for_NEJM/T6. Cohort_important_cna_genes_KEGG.txt",quote=F)
###

intersect(genes_cnv_23_single_cell_important,genes_cnv_cohort_important) %>>%
#,genes_snv_23_single_cell_important)
  writeLines("~/Documents/cohort/ESCC/for_NEJM/T7. Cohort_and_Single-cell_important_cna_genes.txt")
rm(temp)
# ------

# Heatmap of MUC6, PTEN, (881), (884) for cohort
# ------
snv_muc6_pten_cohort <- snv_cohort_gene[c("MUC6","PTEN"),]
temp <- read.table("~/Documents/cohort/ESCC/cnv/cnvr_table.tsv",sep="\t",header=T,row.names=1)
temp <- apply(temp[,1:3],1,function(x) paste(x,collapse="_")) %>>% 
  sapply(function(x) gsub(" +","",x),simplify=T) %>>% (cbind(temp,.))
temp_cyto <- read.table("~/Documents/cohort/ESCC/cnv/cnvr_annovar.hg19_cytoBand",sep="\t",header=F)
temp_cyto <- apply(temp_cyto[,3:5],1,function(x) paste0("chr",paste(x,collapse="_"))) %>>% 
  sapply(function(x) gsub(" +","",x),simplify=T) %>>% (cbind(temp_cyto,.))
temp <- cbind(temp,temp_cyto[match(temp[,161],temp_cyto[,9]),c(2,9)])
temp <- temp[!is.na(temp[,162]),c(1:160,162)]
rownames(temp) <- apply(temp[,1:3],1,function(x) paste(x,collapse="_")) %>>% (gsub(" +","",.))
a <- temp[which(rownames(temp) == "chr19_37443901_38976746"),161]
b <- temp[which(rownames(temp) == "chr19_43589739_45922981"),161]
cnv_881_884_cohort <- rbind(cnv_cohort[grep(a,rownames(cnv_cohort)),],cnv_cohort[grep(b,rownames(cnv_cohort)),])
temp_group <- c(as.character(filter(pi,level == "primary")$病人姓名),
                as.character(filter(pi,level == "therapy")$病人姓名),
                as.character(filter(pi,level == "therapy_after")$病人姓名),
                as.character(filter(pi,level == "relapse")$病人姓名))
temp_group <- temp_group[temp_group != "R1509285860"]
snv_muc6_pten_cohort <- snv_muc6_pten_cohort[,match(temp_group,colnames(snv_muc6_pten_cohort))]
cnv_881_884_cohort <- cnv_881_884_cohort[,match(temp_group,colnames(cnv_881_884_cohort))]
mutation_cohort <- rbind(snv_muc6_pten_cohort,cnv_881_884_cohort)
cluster <- t(mutation_cohort[,1:95]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_cohort <- mutation_cohort[,1:95] %>>% (.[,cluster$order]) %>>% cbind(mutation_cohort[,96:153])
cluster <- t(mutation_cohort[,96:126]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_cohort <- mutation_cohort[,96:126] %>>% (.[,cluster$order]) %>>% 
  (cbind(mutation_cohort[,1:95],.,mutation_cohort[,127:153]))
cluster <- t(mutation_cohort[,127:149]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_cohort <- mutation_cohort[,127:149] %>>% (.[,cluster$order]) %>>% 
  (cbind(mutation_cohort[,1:126],.,mutation_cohort[,150:153]))
cluster <- t(mutation_cohort[,150:153]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_cohort <- mutation_cohort[,150:153] %>>% (.[,cluster$order]) %>>% (cbind(mutation_cohort[,1:149],.))
# Final version
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",95),rep("Therapy",31),rep("Therapy_post",23),rep("Relapse",4)),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy"="#E199FF",
                                               "Therapy_post"="#00DAE0",
                                               "Relapse"="#FF9289")),
                        annotation_height=0.5,
                        show_annotation_name=T,
                        annotation_name_gp=gpar(fontsize=11))
Heatmap(mutation_cohort,col=colorRamp2(c(-1,0,1,2),c("cornflowerblue","gray98","violet","violetred"),space="sRGB"),
        show_row_names=T,show_column_names=F,row_names_gp=gpar(fontsize=11),cluster_columns=F,
        clustering_distance_rows="canberra",clustering_method_rows="ward.D2",row_dend_reorder=F,
        show_row_dend=F,top_annotation=ha,top_annotation_height=unit(10,"pt"),name="mutation",
        rect_gp=gpar(col="gray98"),heatmap_legend_param=list(color_bar="continuous"))
rm(temp)
rm(temp_cyto)
rm(temp_group)
rm(a)
rm(b)
# ------

# Reading SNV/INDEL results (Germline + Somatic)
# ------
snv_indel_cohort <- read.csv("~/Documents/cohort/ESCC/snv_gt.csv",header=T,row.names=1)
temp <- read.csv("~/Documents/cohort/ESCC/indel_gt.csv",header=T,row.names=1)
snv_indel_cohort <- rbind(snv_indel_cohort,temp)
rm(temp)
# Transform SNV from sites to genes
snv_indel_cohort_gene <- sapply(unique(sub("(.*)_(.*)_(.*)","\\3",rownames(snv_indel_cohort))),function(x) {
  a <- snv_indel_cohort[grep(x,rownames(snv_indel_cohort)),]
  b <- apply(a,2,function(y) length(y[y != "." & y != "0/0"]))
  return(b)
}) %>>% t
temp_group <- c(as.character(filter(pi,level == "primary")$病人姓名),
                as.character(filter(pi,level == "therapy")$病人姓名),
                as.character(filter(pi,level == "therapy_after")$病人姓名),
                as.character(filter(pi,level == "relapse")$病人姓名))
temp_group <- temp_group[temp_group != "R1509285860"]
snv_indel_cohort <- snv_indel_cohort[,match(temp_group,colnames(snv_indel_cohort))]
snv_indel_cohort_gene <- as.data.frame(snv_indel_cohort_gene[,match(temp_group,
                                                                    colnames(snv_indel_cohort_gene))])
rm(temp_group)
# ------

# Heatmap of SNV/INDEL
# ------
cluster <- t(snv_indel_cohort_gene[,1:95]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
snv_indel_cohort_gene <- snv_indel_cohort_gene[,1:95] %>>% (.[,cluster$order]) %>>% cbind(snv_indel_cohort_gene[,96:153])
cluster <- t(snv_indel_cohort_gene[,96:126]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
snv_indel_cohort_gene <- snv_indel_cohort_gene[,96:126] %>>% (.[,cluster$order]) %>>%
  (cbind(snv_indel_cohort_gene[,1:95],.,snv_indel_cohort_gene[,127:153]))
cluster <- t(snv_indel_cohort_gene[,127:149]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
snv_indel_cohort_gene <- snv_indel_cohort_gene[,127:149] %>>% (.[,cluster$order]) %>>% 
  (cbind(snv_indel_cohort_gene[,1:126],.,snv_indel_cohort_gene[,150:153]))
cluster <- t(snv_indel_cohort_gene[,150:153]) %>>% dist(method="canberra") %>>% hclust(method="ward.D2")
snv_indel_cohort_gene <- snv_indel_cohort_gene[,150:153] %>>% (.[,cluster$order]) %>>% (cbind(snv_indel_cohort_gene[,1:149],.))
# Final version
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",95),rep("Therapy",31),rep("Therapy_post",23),rep("Relapse",4)),
                        # Sub_groups=anno_points(rep(1,153),
                        #                        pch=c(rep(7,27),rep(8,68),rep(8,10),rep(7,2),rep(8,19),
                        #                              rep(8,3),rep(7,3),rep(8,17),rep(7,1),rep(8,3)),
                        #                        size=unit(5,"pt"),border=F,ylim=c(0.5,1.5),
                        #                        gp=gpar(col=c(rep("lightpink",27),rep("forestgreen",68),
                        #                                      rep("forestgreen",10),rep("lightpink",2),
                        #                                      rep("forestgreen",19),rep("forestgreen",3),
                        #                                      rep("lightpink",3),rep("forestgreen",17),
                        #                                      rep("lightpink",1),rep("forestgreen",3)))),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy"="#E199FF",
                                               "Therapy_post"="#00DAE0",
                                               "Relapse"="#FF9289")),
                        annotation_height=0.5,
                        show_annotation_name=T,
                        # gap=unit(3,"pt"),
                        annotation_name_gp=gpar(fontsize=8))
ht <- Heatmap(snv_indel_cohort_gene,col=colorRamp2(c(0,3),c("gray98","red"),space="sRGB"),
              show_row_names=F,show_column_names=F,row_names_gp=gpar(fontsize=1),cluster_columns=F,
              clustering_distance_rows="maximum",clustering_method_rows="ward.D2",row_dend_reorder=F,
              show_row_dend=F,top_annotation=ha,top_annotation_height=unit(20,"pt"),name="SNV",
              heatmap_legend_param=list(color_bar="continuous"))
# lgd <- Legend(c("A","B"),title="Sub_groups",type="points",pch=7:8,background="white",
              # legend_gp=gpar(col=c("lightpink","forestgreen")))
# draw(ht,annotation_legend_list=list(lgd))
draw(ht)
# ------

# Heatmap of (SNV/INDEL) MUC6, PTEN, (881), (884) for cohort
# ------
snv_indel_muc6_pten_cohort <- snv_indel_cohort_gene[c("MUC6","PTEN"),]
temp_group <- c(as.character(filter(pi,level == "primary")$病人姓名),
                as.character(filter(pi,level == "therapy")$病人姓名),
                as.character(filter(pi,level == "therapy_after")$病人姓名),
                as.character(filter(pi,level == "relapse")$病人姓名))
temp_group <- temp_group[temp_group != "R1509285860"]
snv_indel_muc6_pten_cohort <- snv_indel_muc6_pten_cohort[,match(temp_group,colnames(snv_indel_muc6_pten_cohort))]
cnv_881_884_cohort <- cnv_881_884_cohort[,match(temp_group,colnames(cnv_881_884_cohort))]
mutation_indel_cohort <- rbind(snv_indel_muc6_pten_cohort,cnv_881_884_cohort)
cluster <- t(mutation_indel_cohort[,1:95]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_indel_cohort <- mutation_indel_cohort[,1:95] %>>% (.[,cluster$order]) %>>% cbind(mutation_indel_cohort[,96:153])
cluster <- t(mutation_indel_cohort[,96:126]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_indel_cohort <- mutation_indel_cohort[,96:126] %>>% (.[,cluster$order]) %>>% 
  (cbind(mutation_indel_cohort[,1:95],.,mutation_indel_cohort[,127:153]))
cluster <- t(mutation_indel_cohort[,127:149]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_indel_cohort <- mutation_indel_cohort[,127:149] %>>% (.[,cluster$order]) %>>% 
  (cbind(mutation_indel_cohort[,1:126],.,mutation_indel_cohort[,150:153]))
cluster <- t(mutation_indel_cohort[,150:153]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_indel_cohort <- mutation_indel_cohort[,150:153] %>>% (.[,cluster$order]) %>>% 
  (cbind(mutation_indel_cohort[,1:149],.))
# Final version
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",95),rep("Therapy",31),rep("Therapy_post",23),rep("Relapse",4)),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy"="#E199FF",
                                               "Therapy_post"="#00DAE0",
                                               "Relapse"="#FF9289")),
                        annotation_height=0.5,
                        show_annotation_name=T,
                        annotation_name_gp=gpar(fontsize=11))
Heatmap(mutation_indel_cohort,col=colorRamp2(c(-1,0,1,2),c("cornflowerblue","gray98","violet","violetred"),space="sRGB"),
        show_row_names=T,show_column_names=F,row_names_gp=gpar(fontsize=11),cluster_columns=F,
        clustering_distance_rows="binary",clustering_method_rows="ward.D2",row_dend_reorder=F,
        show_row_dend=F,top_annotation=ha,top_annotation_height=unit(10,"pt"),name="mutation",
        rect_gp=gpar(col="gray98"),heatmap_legend_param=list(color_bar="continuous"))
rm(temp)
rm(temp_cyto)
rm(temp_group)
rm(a)
rm(b)
# Modification
mutation_indel_cohort_combine <- rbind(MUC6_PTEN=apply(mutation_indel_cohort[1:2,],
                                                       2,
                                                       function(x) 
                                                         if (length(x[x != 0]) == 2) 
                                                           return(1)
                                                         else
                                                           return(0)),
                                       mutation_indel_cohort[3:6,])
cluster <- t(mutation_indel_cohort_combine[,1:95]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_indel_cohort_combine <- mutation_indel_cohort_combine[,1:95] %>>% (.[,cluster$order]) %>>% 
  cbind(mutation_indel_cohort_combine[,96:153])
cluster <- t(mutation_indel_cohort_combine[,96:126]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_indel_cohort_combine <- mutation_indel_cohort_combine[,96:126] %>>% (.[,cluster$order]) %>>% 
  (cbind(mutation_indel_cohort_combine[,1:95],.,mutation_indel_cohort_combine[,127:153]))
cluster <- t(mutation_indel_cohort_combine[,127:149]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_indel_cohort_combine <- mutation_indel_cohort_combine[,127:149] %>>% (.[,cluster$order]) %>>% 
  (cbind(mutation_indel_cohort_combine[,1:126],.,mutation_indel_cohort_combine[,150:153]))
cluster <- t(mutation_indel_cohort_combine[,150:153]) %>>% dist(method="binary") %>>% hclust(method="ward.D2")
mutation_indel_cohort_combine <- mutation_indel_cohort_combine[,150:153] %>>% (.[,cluster$order]) %>>% 
  (cbind(mutation_indel_cohort_combine[,1:149],.))
temp <- mutation_indel_cohort_combine
cluster <- dist(temp,method="euclidean") %>>% hclust(method="ward.D") 
temp <- temp[cluster$order,]
temp[c(3,5),] <- adply(temp[c(3,5),],1,function(x) {
  x[x != 0] <- 1
  return(x)
})
temp[1,] <- sapply(temp[1,],function(x) {
  x[x != 0] <- 2
  return(x)
})
temp[c(2,4),] <- adply(temp[c(2,4),],1,function(x) {
  x[x != 0] <- 3
  return(x)
})
ha <- HeatmapAnnotation(Irradiation=c(rep("Primary",95),rep("Therapy",31),rep("Therapy_post",23),rep("Relapse",4)),
                        col=list(Irradiation=c("Primary"="#96CA00",
                                               "Therapy"="#E199FF",
                                               "Therapy_post"="#00DAE0",
                                               "Relapse"="#FF9289")),
                        annotation_height=0.5,
                        show_annotation_name=T,
                        annotation_name_gp=gpar(fontsize=11,fontface="bold",fontfamily="serif"),
                        annotation_legend_param=list(title_gp=gpar(fontface="bold",fontfamily="serif"),
                                                     labels_gp=gpar(fontfamily="serif")))
Heatmap(temp,
        col=colorRamp2(c(0,1,2,3),c("gray98","cornflowerblue","yellowgreen","violetred"),space="sRGB"),
        show_row_names=T,show_column_names=F,row_names_gp=gpar(fontsize=11,fontfamily="serif"),cluster_columns=F,
        cluster_rows=F,show_row_dend=F,top_annotation=ha,top_annotation_height=unit(10,"pt"),name="Mutation",
        rect_gp=gpar(col="gray98"),heatmap_legend_param=list(color_bar="discrete",
                                                             labels=c("","Deletion","SNV","Amplification"),
                                                             labels_gp=gpar(fontfamily="serif"),
                                                             title_gp=gpar(fontface="bold",fontfamily="serif")))
# ------

# Venn diagrams
# 18*18 pdf
# ------
library(VennDiagram)
draw.triple.venn(21,15113,3279,16,2407,8,7,
                 c("Radiotherapy-related genes in single-cell SNV",
                   "Radiotherapy-related genes in single-cell CNA",
                   "Radiotherapy-related genes in cohort CNA"),
                 fill=c("#00468BFF","#ED0000FF","#42B540FF"),
                 cat.col=c("#00468BFF","#ED0000FF","#42B540FF"),
                 cex=6,cat.cex=2.5,cat.pos=c(-15,18,180),lty=1:3)
# 6 * 6 pdf
draw.quad.venn(21,15113,3279,14694,16,8,18,2407,10562,2457,7,14,7,1853,6,
               c("Radiotherapy-related genes in single-cell SNV",
                 "Radiotherapy-related genes in single-cell CNA",
                 "Radiotherapy-related genes in cohort CNA",
                 "All mutant genes in cohort SNV"),
               fill=c("cornflowerblue","violetred","yellowgreen","orangered1"),
               cat.col=c("cornflowerblue","violetred","yellowgreen","orangered1"),
               cat.pos=c(-212,202,0,45),cat.dist=c(0.50,0.50,0.11,0.11),cat.cex=rep(1.2,4),
               lty=1:4,fontfamily=rep("sans",15),cat.fontfamily=rep("sans",4))
# ------

# ============================================================================

# Unifying data/method of analysis
# Which cnv in 23 cm was covered by cnv_cohort_region
# ------
rownames(cnv_23) <- apply(cnv_23[,1:3],1,function(x) paste(x[1:3],collapse="_") %>>% (gsub(" +","",.)))
c23 <- strsplit(rownames(cnv_23),"_",fixed=T) %>>% as.data.frame %>>% t
colnames(c23) <- c("Chr","Start","End")
c23 <- as.data.frame(c23)
c23$Start <- as.character(c23$Start) %>>% as.integer
c23$End <- as.character(c23$End) %>>% as.integer
# Parallelizing calculation
library(doSNOW)
cl <- makeCluster(4,type="SOCK")
clusterExport(cl,list("cnv_23","c23"))
registerDoSNOW(cl)
system("date")
result <- parLapply(cl,rownames(cnv_cohort_region),function(x) {
  chr <- strsplit(x,"_",fixed=T)[[1]][1]
  start <- strsplit(x,"_",fixed=T)[[1]][2] %>>% as.integer
  end <- strsplit(x,"_",fixed=T)[[1]][3] %>>% as.integer
  r <- filter(c23,Chr == chr & Start >= start & End <= end) %>>% list
  names(r) <- paste(chr,start,end,sep="_")
  return(r)
})
system("date")
stopCluster(cl)
# ------

# Which cnv in cohort was covered by cnv_23
# ------
c_cohort <- strsplit(rownames(cnv_cohort_region),"_",fixed=T) %>>% as.data.frame %>>% t
colnames(c_cohort) <- c("Chr","Start","End")
c_cohort <- as.data.frame(c_cohort)
c_cohort$Start <- as.character(c_cohort$Start) %>>% as.integer
c_cohort$End <- as.character(c_cohort$End) %>>% as.integer
# Parallelizing calculation
library(doSNOW)
cl <- makeCluster(4,type="SOCK")
clusterExport(cl,list("cnv_cohort_region","c_cohort"))
registerDoSNOW(cl)
system("date")
result2 <- parLapply(cl,rownames(cnv_23),function(x) {
  chr <- strsplit(x,"_",fixed=T)[[1]][1]
  start <- strsplit(x,"_",fixed=T)[[1]][2] %>>% as.integer
  end <- strsplit(x,"_",fixed=T)[[1]][3] %>>% as.integer
  r <- filter(c_cohort,Chr == chr & Start >= start & End <= end) %>>% list
  names(r) <- paste(chr,start,end,sep="_")
  return(r)
})
system("date")
stopCluster(cl)
# ------

# ------
result <- sapply(result,function(x) nrow(x[1][[1]])) %>>% (result[which(. != 0)])
result2 <- sapply(result2,function(x) nrow(x[1][[1]])) %>>% (result2[which(. != 0)])
cohort_cover_sc23 <- 
  ldply(result,function(x) return(c(names(x[1]),
                                    paste(apply(x[[1]],1,function(y) paste(y,collapse="_")),collapse=","))))
sc23_cover_cohort <- 
  ldply(result2,function(x) return(c(names(x[1]),
                                     paste(apply(x[[1]],1,function(y) paste(y,collapse="_")),collapse=","))))
colnames(cohort_cover_sc23) <- c("CNA_cohort","CNA_covered_in_sc23")
colnames(sc23_cover_cohort) <- c("CNA_sc23","CNA_covered_in_cohort")
write.csv(cohort_cover_sc23,file="~/Documents/cohort/ESCC/for_NEJM/CNA_cohort_cover_sc23.csv")
write.csv(sc23_cover_cohort,file="~/Documents/cohort/ESCC/for_NEJM/CNA_sc23_cover_cohort.csv")
# ------

# Schematic drawing
# Single-cell, male
# ------
library(showtext)
library(cowplot)
font.add("xf","~/Public/fonts/Hexafont.ttf")
sc23 <- data.frame(Sample=factor(c(rep("Control",10),rep("Dysplasia",15),rep("Cancer",30))),
                   Cell=c(1:10,1:15,1:30),Char=c(rep("O",10),rep("Q",15),rep("R",30)))
sc23$Sample <- relevel(sc23$Sample,"Dysplasia")
sc23$Sample <- relevel(sc23$Sample,"Cancer")
pdf("test.pdf",width=3,height=7)
showtext.begin()
ggplot(sc23,aes(Cell,Sample)) + 
  geom_text(aes(label=Char,color=Sample,family="xf",size=56),position="jitter") +
  theme(legend.position="none",axis.text.x=element_text(angle=45,vjust=0.5)) + 
  scale_x_continuous(name="Cell numbers",trans=scales::reverse_trans()) + 
  scale_y_discrete(name="",labels=c("Cancer cells","Dysplasia cells","Control cells")) +
  coord_flip()
showtext.end()
dev.off()
# ------ 
# halted

# Common SNV genes (Sensitive, Resistant, De novo) in cohort, pathways...
# ------
pi <- pi[pi[,3] != "R1509285860",]
common_genes <- 
  read.csv("~/Documents/cohort/ESCC/for_NEJM/SNV_sc23 and cohort_germline_and_somatic.csv",row.names=1)[,1] %>>% 
  as.character
common_genes_primary <- snv_indel_cohort_gene[common_genes,] %>>% (.[,match(pi[pi$level == "primary",3],colnames(.))])
common_genes_therapy_aftertherapy <- 
  snv_indel_cohort_gene[common_genes,] %>>% 
  (.[,match(pi[pi$level == "therapy" | pi$level == "therapy_after",3],colnames(.))])
common_genes_up <- sapply(common_genes,function(x) {
  a <- common_genes_primary[x,] %>>% as.integer %>>% (length(.[. != 0]))
  b <- common_genes_primary[x,] %>>% as.integer %>>% (length(.[. == 0]))
  c <- common_genes_therapy_aftertherapy[x,] %>>% as.integer %>>% (length(.[. != 0]))
  d <- common_genes_therapy_aftertherapy[x,] %>>% as.integer %>>% (length(.[. == 0]))
  fisher.test(matrix(c(a,c,b,d),nrow=2),alternative="less") %>>% (return(.$p.val))
}) %>>% (which(. < 0.05))
write.csv(names(common_genes_up),file="~/Documents/cohort/ESCC/for_NEJM/common_genes_up.csv")
common_genes_down <- sapply(common_genes,function(x) {
  a <- common_genes_primary[x,] %>>% as.integer %>>% (length(.[. != 0]))
  b <- common_genes_primary[x,] %>>% as.integer %>>% (length(.[. == 0]))
  c <- common_genes_therapy_aftertherapy[x,] %>>% as.integer %>>% (length(.[. != 0]))
  d <- common_genes_therapy_aftertherapy[x,] %>>% as.integer %>>% (length(.[. == 0]))
  fisher.test(matrix(c(a,c,b,d),nrow=2),alternative="greater") %>>% (return(.$p.val))
}) %>>% (which(. < 0.05))
write.csv(names(common_genes_down),file="~/Documents/cohort/ESCC/for_NEJM/common_genes_down.csv")
common_genes_insignificant <- common_genes[-c(common_genes_up,common_genes_down)]
write.csv(common_genes_insignificant,file="~/Documents/cohort/ESCC/for_NEJM/common_genes_insignificant.csv")
k <- bitr(common_genes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>% enrichKEGG %>>% 
  as.data.frame
k$geneID <- sapply(k$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
               (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
               paste(collapse="/"))
write.csv(k,file="~/Documents/cohort/ESCC/for_NEJM/common_genes_kegg.csv")
k <- bitr(names(common_genes_up),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>% 
  enrichKEGG %>>% as.data.frame
k$geneID <- sapply(k$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
               (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
               paste(collapse="/"))
write.csv(k,file="~/Documents/cohort/ESCC/for_NEJM/common_genes_up_kegg.csv")
k <- bitr(names(common_genes_down),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>% 
  enrichKEGG %>>% as.data.frame
k$geneID <- sapply(k$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
               (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
               paste(collapse="/"))
write.csv(k,file="~/Documents/cohort/ESCC/for_NEJM/common_genes_down_kegg.csv")
k <- bitr(common_genes_insignificant,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>% 
  enrichKEGG %>>% as.data.frame
k$geneID <- sapply(k$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
               (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
               paste(collapse="/"))
write.csv(k,file="~/Documents/cohort/ESCC/for_NEJM/common_genes_insignificant_kegg.csv")
# ------

# Common CNV analysis
# ------
cna_cohort_cover_23 <- read.csv("~/Documents/cohort/ESCC/for_NEJM/CNA_cohort_cover_sc23.csv",as.is=T,row.names=1) %>>% 
  rev
colnames(cna_cohort_cover_23) <- c("sc23","cohort")
cna_23_cover_cohort <- read.csv("~/Documents/cohort/ESCC/for_NEJM/CNA_sc23_cover_cohort.csv",as.is=T,row.names=1)
colnames(cna_23_cover_cohort) <- c("sc23","cohort")
cna_common <- rbind(cna_cohort_cover_23,cna_23_cover_cohort)
write.csv(c(cna_common[1:346,1],cna_common[347:nrow(cna_common),2]),
          file="~/Documents/cohort/ESCC/for_NEJM/cna_common_for_annovar.avinput",row.names=F,quote=F)
# After annotation by ANNOVAR
cna_common_genes <- 
  read.table("~/Documents/cohort/ESCC/for_NEJM/cna_common_for_annovar.avoutput.hg19_bed",header=F,sep="\t") %>>% 
  (sub("Name=(.*)","\\1",.[,2])) %>>% sapply(function(x) strsplit(x,",",fixed=T)[[1]]) %>>% unlist %>>% unique
k <- bitr(cna_common_genes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>% enrichKEGG %>>% 
  as.data.frame
k$geneID <- sapply(k$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
               (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
               paste(collapse="/"))
write.csv(k,file="~/Documents/cohort/ESCC/for_NEJM/cna_common_genes_kegg.csv")
# Fisher.test
common_cna_in_cohort <- lapply(cna_common[,2],function(x) strsplit(x,",",fixed=T)[[1]]) %>>% unlist
common_cna_primary <- cnv_cohort_region[common_cna_in_cohort,] %>>% 
  (.[,match(pi[pi$level == "primary",3],colnames(.)) %>>% (.[!is.na(.)])])
common_cna_therapy_aftertherapy <- cnv_cohort_region[common_cna_in_cohort,] %>>% 
  (.[,match(pi[pi$level == "therapy" | pi$level == "therapy_after",3],colnames(.))])
common_cna_up <- sapply(common_cna_in_cohort,function(x) {
  a <- common_cna_primary[x,] %>>% as.integer %>>% (length(.[. != 0]))
  b <- common_cna_primary[x,] %>>% as.integer %>>% (length(.[. == 0]))
  c <- common_cna_therapy_aftertherapy[x,] %>>% as.integer %>>% (length(.[. != 0]))
  d <- common_cna_therapy_aftertherapy[x,] %>>% as.integer %>>% (length(.[. == 0]))
  fisher.test(matrix(c(a,c,b,d),nrow=2),alternative="less") %>>% (return(.$p.val))
}) %>>% (which(. < 0.05))
# Looking for the overlapping region in both cohort and single-cell
common_cna_up2 <- cna_common[match(names(common_cna_up),cna_common[,2]),1] %>>% (gsub(" +","",.)) %>>% 
  (strsplit(.,",",fixed=T)[[1]])
common_cna_down <- sapply(common_cna_in_cohort,function(x) {
  a <- common_cna_primary[x,] %>>% as.integer %>>% (length(.[. != 0]))
  b <- common_cna_primary[x,] %>>% as.integer %>>% (length(.[. == 0]))
  c <- common_cna_therapy_aftertherapy[x,] %>>% as.integer %>>% (length(.[. != 0]))
  d <- common_cna_therapy_aftertherapy[x,] %>>% as.integer %>>% (length(.[. == 0]))
  fisher.test(matrix(c(a,c,b,d),nrow=2),alternative="greater") %>>% (return(.$p.val))
}) %>>% (which(. < 0.05))
common_cna_insignificant <- common_cna_in_cohort[-c(common_cna_up,common_cna_down)]
# Looking for the overlapping region in both cohort and single-cell
common_cna_insignificant <- lapply(common_cna_insignificant,function(x) {
  a <- grep(x,cna_common[,2])
  if (a <= 346) {
    return(cna_common[a,1])
  } else {
    return(x)
  }
}) %>>% unlist %>>% paste(collapse=",") %>>% (strsplit(.,",",fixed=T)[[1]]) %>>% unique
write.csv(common_cna_up2,file="~/Documents/cohort/ESCC/for_NEJM/common_cna_up.avinput",
          quote=F,row.names=F)
# Nothing for enrichment
# write.csv(names(common_cna_down),file="~/Documents/cohort/ESCC/for_NEJM/common_cna_down.avinput",
#           quote=F,row.names=F)
write.csv(common_cna_insignificant,file="~/Documents/cohort/ESCC/for_NEJM/common_cna_insignificant.avinput",
          quote=F,row.names=F)
common_cna_up_genes <- 
  read.table("~/Documents/cohort/ESCC/for_NEJM/common_cna_up.avoutput.hg19_bed",header=F,sep="\t") %>>% 
  (sub("Name=(.*)","\\1",.[,2])) %>>% lapply(function(x) strsplit(x,",",fixed=T)[[1]]) %>>% unlist %>>% unique
common_cna_insignificant_genes <- 
  read.table("~/Documents/cohort/ESCC/for_NEJM/common_cna_insignificant.avoutput.hg19_bed",header=F,sep="\t") %>>% 
  (sub("Name=(.*)","\\1",.[,2])) %>>% lapply(function(x) strsplit(x,",",fixed=T)[[1]]) %>>% unlist %>>% unique
k <- bitr(common_cna_up_genes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>% enrichKEGG %>>% 
  as.data.frame
k$geneID <- sapply(k$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
               (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
               paste(collapse="/"))
write.csv(k,file="~/Documents/cohort/ESCC/for_NEJM/common_cna_up_genes_kegg.csv")
k <- bitr(common_cna_insignificant_genes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID %>>% 
  enrichKEGG %>>% as.data.frame
k$geneID <- sapply(k$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
               (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
               paste(collapse="/"))
write.csv(k,file="~/Documents/cohort/ESCC/for_NEJM/common_cna_insignificant_genes_kegg.csv")
# ------

# Oncogene or TSG analysis
# ------
cosmic_census <- read.csv("~/Documents/cohort/ESCC/for_NEJM/Census_allThu Jul 20 04-34-56 2017.csv")
cosmic_census[match(common_genes,cosmic_census$Gene.Symbol),c("Gene.Symbol","Role.in.Cancer")] %>>% 
  (.[!is.na(.[,2]) & .[,2] != "",]) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_genes_oncogene_TGS.csv",row.names=F,quote=F)
cosmic_census[match(names(common_genes_up),cosmic_census$Gene.Symbol),c("Gene.Symbol","Role.in.Cancer")] %>>% 
  (.[!is.na(.[,2]) & .[,2] != "",]) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_genes_up_oncogene_TGS.csv",row.names=F,quote=F)
cosmic_census[match(names(common_genes_down),cosmic_census$Gene.Symbol),c("Gene.Symbol","Role.in.Cancer")] %>>% 
  (.[!is.na(.[,2]) & .[,2] != "",]) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_genes_down_oncogene_TGS.csv",row.names=F,quote=F)
cosmic_census[match(common_genes_insignificant,cosmic_census$Gene.Symbol),
              c("Gene.Symbol","Role.in.Cancer")] %>>% (.[!is.na(.[,2]) & .[,2] != "",]) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_genes_insignificant_oncogene_TGS.csv",row.names=F,quote=F)
cosmic_census[match(cna_common_genes,cosmic_census$Gene.Symbol),c("Gene.Symbol","Role.in.Cancer")] %>>% 
  (.[!is.na(.[,2]) & .[,2] != "",]) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_cna_oncogene_TGS.csv",row.names=F,quote=F)
cosmic_census[match(common_cna_up_genes,cosmic_census$Gene.Symbol),c("Gene.Symbol","Role.in.Cancer")] %>>% 
  (.[!is.na(.[,2]) & .[,2] != "",]) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_cna_up_oncogene_TGS.csv",row.names=F,quote=F)
cosmic_census[match(common_cna_insignificant_genes,cosmic_census$Gene.Symbol),c("Gene.Symbol","Role.in.Cancer")] %>>% 
  (.[!is.na(.[,2]) & .[,2] != "",]) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_cna_insignificant_oncogene_TGS.csv",row.names=F,quote=F)
# ------

# Oncogenes/TGS in CNA
# ------
cna_common_cover_genes <- 
  read.table("~/Documents/cohort/ESCC/for_NEJM/cna_common_for_annovar.avoutput.hg19_bed",header=F,sep="\t")
cna_common_cover_genes <- cna_common_cover_genes[,2:5]
cna_common_cover_genes[,5] <- apply(cna_common_cover_genes[,2:4],1,
                                    function(x) paste0("chr",paste(x,collapse="_")))
cna_common_cover_genes <- cna_common_cover_genes[!duplicated(cna_common_cover_genes[,5]),c(5,1)]
library(stringr)
cna_common_cover_genes[,2] <- sapply(cna_common_cover_genes[,2],function(x) str_sub(x,6))
# tail(cna_common_cover_genes)
# cnv_cohort_region["chr9_67269918_67271285",] %>>% as.integer %>>% table

common_cna_in_cohort_genes <- read.table("~/Documents/cohort/ESCC/for_NEJM/common_cna_in_cohort.avoutput.hg19_bed",
                                         sep="\t",header=F,as.is=T)
common_cna_in_cohort_genes <- common_cna_in_cohort_genes[,c(3:5,2)]
rownames(common_cna_in_cohort_genes) <- 
  apply(common_cna_in_cohort_genes[,1:3],1,function(x) paste0("chr",paste(x,collapse="_")))
common_cna_in_cohort_genes[,4] <- sapply(common_cna_in_cohort_genes[,4],function(x) str_sub(x,6))
temp <- lapply(common_cna_in_cohort_genes[,4],function(x) {
  a <- strsplit(x,",",fixed=T)[[1]]
  try(b <- bitr(a,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID,TRUE)
  if (exists("b")) b <- enrichKEGG(b) %>>% as.data.frame else return(NA)
  if (nrow(b) == 0) return(NA)
  b$geneID <- sapply(b$geneID,function(x) strsplit(x,"/",fixed=T)[[1]] %>>% 
                       (bitr(.,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) %>>% 
                       paste(collapse="/") %>>% return)
  return(b)
})
common_cna_in_cohort_genes[,5] <- lapply(temp,function(x) if 
                                         (is.data.frame(x)) apply(x,1,function(x) paste(x[2],x[8],sep="_")) %>>% 
                                           unlist %>>% paste(collapse=",") else return(NA)) %>>% unlist
common_cna_in_cohort_genes[,6] <- lapply(common_cna_in_cohort_genes[,4],function(x) {
  a <- strsplit(x,",",fixed=T)[[1]]
  b <- cosmic_census[match(a,cosmic_census$Gene.Symbol),c("Gene.Symbol","Role.in.Cancer")] %>>% 
    (.[!is.na(.[,2]) & .[,2] != "",]) %>>% apply(1,function(x) paste(x,collapse="_")) %>>% unlist %>>% 
    paste(collapse=",")
  if (length(b) == 0) return(NA) else return(b)
}) %>>% unlist
common_cna_in_cohort <- data.frame(CNA=common_cna_in_cohort)
common_cna_in_cohort <- cbind(common_cna_in_cohort$CNA,
                              common_cna_in_cohort_genes[match(common_cna_in_cohort$CNA,
                                                               rownames(common_cna_in_cohort_genes)),])
common_cna_in_cohort <- common_cna_in_cohort[,c(1,5:7)]
common_cna_in_cohort[,4] <- sapply(common_cna_in_cohort[,4],
                                   function(x) if (!is.na(x) & x == "") return(NA) else return(x))
colnames(common_cna_in_cohort) <- c("CNA","Gene","KEGG","Cosmic")
write.table(common_cna_in_cohort,file="~/Documents/cohort/ESCC/for_NEJM/common_cna_in_cohort.tsv",
            quote=F,sep="\t",row.names=F)
# ------
# halted

# CNA plot for 12p12.1/16q21 
# ======
cnv_cohort_region_chr12 <- grep("^chr12",rownames(cnv_cohort_region)) %>>% 
  (cnv_cohort_region[.,]) 
# Remove out-regions
# ---
cnv_cohort_region_chr12 <- 
  cbind(cnv_cohort_region_chr12,
        start=(str_split(rownames(cnv_cohort_region_chr12),"_") %>>% 
          sapply(function(x) x[2]) %>>% as.integer)) %>>%
  cbind(end=(str_split(rownames(cnv_cohort_region_chr12),"_") %>>% 
               sapply(function(x) x[3]) %>>% as.integer))
cnv_cohort_region_chr12 <- 
  cnv_cohort_region_chr12[-which((cnv_cohort_region_chr12$start < 60000 &
                                    cnv_cohort_region_chr12$end < 60000) | 
                                   (cnv_cohort_region_chr12$start > 133841895 &
                                      cnv_cohort_region_chr12$end > 133841895)),
                          ] # 133851895 - 10000 = 133841895
cnv_cohort_region_chr12$start <- cnv_cohort_region_chr12$start %>>%
  {.[. < 60000] <- 60000;.}
cnv_cohort_region_chr12$end <- cnv_cohort_region_chr12$end %>>%
{.[. > 133841895] <- 133841895;.}
rownames(cnv_cohort_region_chr12) <- paste("chr12",
                                           cnv_cohort_region_chr12$start,
                                           cnv_cohort_region_chr12$end,
                                           sep="_")
cnv_cohort_region_chr12 <- cnv_cohort_region_chr12[,1:153]
# ---
sp <- split(1:133851895,ceiling(seq_along(1:133851895)/1338519))
sp <- sapply(sp,function(x) return(c(x[1],x[length(x)]))) %>>% t %>>% 
  as.data.frame
colnames(sp) <- c("start","end")
average_copyNumber <- apply(sp,1,function(x) {
  a <- sapply(rownames(cnv_cohort_region_chr12),
              function(x) str_split(x,"_")[[1]][2:3] %>>% as.integer) %>>% t
  a <- cnv_cohort_region_chr12[which((a[,1] >= x[1] & a[,2] <= x[2]) |
                                       (a[,1] < x[1] & a[,2] >= x[1]) |
                                       (a[,1] <= x[2] & a[,2] > x[2])),]
  b <- as.matrix(a[,as.character(subset(pi,level == "primary")$病人姓名)]) %>>% 
    mean
  c <- as.matrix(a[,as.character(subset(pi,level == "therapy")$病人姓名)]) %>>% 
    mean
  d <- a[,as.character(subset(pi,level == "therapy_after")$病人姓名)] %>>% 
    as.matrix %>>% mean
  e <- as.matrix(a[,as.character(subset(pi,level == "relapse")$病人姓名)]) %>>% 
    mean
  c(b,c,d,e)
}) %>>% t
average_copyNumber <- as.data.frame(average_copyNumber)
colnames(average_copyNumber) <- c("Primary","Therapy","PostTherapy","Relapse")
average_copyNumber <- apply(average_copyNumber,1,function(x) {
  x[is.nan(x)] <- 2
  return(x)
}) %>>% t %>>% as.data.frame
sp <- cbind(sp,average_copyNumber)
group <- function(x) {
  a <- rep(3,length(which(sp[,x] > 2)))
  names(a) <- which(sp[,x] > 2)
  b <- rep(2,length(which(sp[,x] == 2)))
  names(b) <- which(sp[,x] == 2)
  c <- rep(1,length(which(sp[,x] < 2)))
  names(c) <- which(sp[,x] < 2)
  y <- c(a,b,c)
  y <- y[order(as.integer(names(y)))] %>>% as.character
  y[y == "1"] <- "Deletion"
  y[y == "2"] <- "Normal"
  y[y == "3"] <- "Amplification"
  y <- as.factor(y)
  y <- relevel(y,"Normal")
  return(y)
}
library(ggthemes)
theme_set(theme_bw(base_family="serif"))
g1 <- ggplot(cbind(sp,group=group("Primary"))) + 
  geom_point(aes(start,Primary,color=group)) + 
  geom_vline(xintercept=24093343,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_lancet() + labs(x=NULL) + 
  scale_x_continuous(breaks=c(1,132513382),labels=NULL) + 
  annotate("text",x=32124457,y=3.25,label="12p12.1",family="serif",size=6,
           color="cornflowerblue") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g2 <- ggplot(cbind(sp,group=group("Therapy"))) + 
  geom_point(aes(start,Therapy,color=group)) + 
  geom_vline(xintercept=24093343,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_lancet() + labs(x=NULL) + 
  scale_x_continuous(breaks=c(1,132513382),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g3 <- ggplot(cbind(sp,group=group("PostTherapy"))) + 
  geom_point(aes(start,PostTherapy,color=group)) + 
  geom_vline(xintercept=24093343,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + 
  scale_color_lancet() + labs(x=NULL) + 
  scale_x_continuous(breaks=c(1,132513382),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g4 <- ggplot(cbind(sp,group=group("Relapse"))) + 
  geom_point(aes(start,Relapse,color=group)) + 
  geom_vline(xintercept=24093343,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_lancet() + 
  labs(x="Chromosome 12") + 
  scale_x_continuous(breaks=c(1,132513382),labels=c(1,133851895)) + 
  scale_y_continuous(breaks=2:4,labels=2:4) + 
  guides(color=guide_legend(title.theme=element_text(face="bold",size=20,
                                                     angle=0,family="serif"),
                            label.theme=element_text(size=18,angle=0,
                                                     family="serif"),
                            title="Copy Number")) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
legend <- get_legend(g4)
prow <- plot_grid(g1 + theme(legend.position="none"),
                  g2 + theme(legend.position="none"),
                  g3 + theme(legend.position="none"),
                  g4 + theme(legend.position="none"),align="v",ncol=1)
plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7))
ggsave("~/Documents/cohort/ESCC/for_NEJM/cna_12p12.1.pdf",
       plot=plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7)),width=10)

# ---

cnv_cohort_region_chr16 <- grep("^chr16",rownames(cnv_cohort_region)) %>>% 
  (cnv_cohort_region[.,])
# Remove out-regions
# ---
cnv_cohort_region_chr16 <- 
  cbind(cnv_cohort_region_chr16,
        start=(str_split(rownames(cnv_cohort_region_chr16),"_") %>>% 
                 sapply(function(x) x[2]) %>>% as.integer)) %>>%
  cbind(end=(str_split(rownames(cnv_cohort_region_chr16),"_") %>>% 
               sapply(function(x) x[3]) %>>% as.integer))
cnv_cohort_region_chr16 <- 
  cnv_cohort_region_chr16[-which((cnv_cohort_region_chr16$start < 60000 &
                                    cnv_cohort_region_chr16$end < 60000) | 
                                   (cnv_cohort_region_chr16$start > 90294753 &
                                      cnv_cohort_region_chr16$end > 90294753)),
                          ] # 90354753 - 60000 = 90294753
cnv_cohort_region_chr16$start <- cnv_cohort_region_chr16$start %>>%
{.[. < 60000] <- 60000;.}
cnv_cohort_region_chr16$end <- cnv_cohort_region_chr16$end %>>%
{.[. > 90294753] <- 90294753;.}
rownames(cnv_cohort_region_chr16) <- paste("chr16",
                                           cnv_cohort_region_chr16$start,
                                           cnv_cohort_region_chr16$end,
                                           sep="_")
cnv_cohort_region_chr16 <- cnv_cohort_region_chr16[,1:153]
# ---
sp <- split(1:90354753,ceiling(seq_along(1:90354753)/903548))
sp <- sapply(sp,function(x) return(c(x[1],x[length(x)]))) %>>% t %>>% 
  as.data.frame
colnames(sp) <- c("start","end")
average_copyNumber <- apply(sp,1,function(x) {
  a <- sapply(rownames(cnv_cohort_region_chr16),
              function(x) str_split(x,"_")[[1]][2:3] %>>% as.integer) %>>% t
  a <- cnv_cohort_region_chr16[which((a[,1] >= x[1] & a[,2] <= x[2]) |
                                       (a[,1] < x[1] & a[,2] >= x[1]) |
                                       (a[,1] <= x[2] & a[,2] > x[2])),]
  b <- as.matrix(a[,as.character(subset(pi,level == "primary")$病人姓名)]) %>>% 
    mean
  c <- as.matrix(a[,as.character(subset(pi,level == "therapy")$病人姓名)]) %>>% 
    mean
  d <- a[,as.character(subset(pi,level == "therapy_after")$病人姓名)] %>>% 
    as.matrix %>>% mean
  e <- as.matrix(a[,as.character(subset(pi,level == "relapse")$病人姓名)]) %>>% 
    mean
  c(b,c,d,e)
}) %>>% t
average_copyNumber <- as.data.frame(average_copyNumber)
colnames(average_copyNumber) <- c("Primary","Therapy","PostTherapy","Relapse")
average_copyNumber <- apply(average_copyNumber,1,function(x) {
  x[is.nan(x)] <- 2
  return(x)
}) %>>% t %>>% as.data.frame
sp <- cbind(sp,average_copyNumber)
group <- function(x) {
  a <- rep(3,length(which(sp[,x] > 2)))
  names(a) <- which(sp[,x] > 2)
  b <- rep(2,length(which(sp[,x] == 2)))
  names(b) <- which(sp[,x] == 2)
  c <- rep(1,length(which(sp[,x] < 2)))
  names(c) <- which(sp[,x] < 2)
  y <- c(a,b,c)
  y <- y[order(as.integer(names(y)))] %>>% as.character
  y[y == "1"] <- "Deletion"
  y[y == "2"] <- "Normal"
  y[y == "3"] <- "Amplification"
  y <- as.factor(y)
  y <- relevel(y,"Normal")
  return(y)
}
theme_set(theme_bw(base_family="serif"))
g1 <- ggplot(cbind(sp,group=group("Primary"))) + 
  geom_point(aes(start,Primary,color=group)) + 
  geom_vline(xintercept=57827073,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_lancet() + labs(x=NULL) + 
  scale_x_continuous(breaks=c(1,89451253),labels=NULL) + 
  annotate("text",x=62344813,y=2.5,label="16q21",family="serif",size=6,
           color="cornflowerblue") + 
  guides(color=guide_legend(title.theme=element_text(face="bold",size=20,
                                                     angle=0,family="serif"),
                            label.theme=element_text(size=18,angle=0,
                                                     family="serif"),
                            title="Copy Number")) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g2 <- ggplot(cbind(sp,group=group("Therapy"))) + 
  geom_point(aes(start,Therapy,color=group)) + 
  geom_vline(xintercept=57827073,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_lancet() + labs(x=NULL) + 
  scale_x_continuous(breaks=c(1,89451253),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g3 <- ggplot(cbind(sp,group=group("PostTherapy"))) + 
  geom_point(aes(start,PostTherapy,color=group)) + 
  geom_vline(xintercept=57827073,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + 
  scale_color_lancet() + labs(x=NULL) + 
  scale_x_continuous(breaks=c(1,89451253),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g4 <- ggplot(cbind(sp,group=group("Relapse"))) + 
  geom_point(aes(start,Relapse,color=group)) + 
  geom_vline(xintercept=57827073,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_lancet() + 
  labs(x="Chromosome 16") + 
  scale_x_continuous(breaks=c(1,89451253),labels=c(1,90354753)) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
legend <- get_legend(g1)
prow <- plot_grid(g1 + theme(legend.position="none"),
                  g2 + theme(legend.position="none"),
                  g3 + theme(legend.position="none"),
                  g4 + theme(legend.position="none"),align="v",ncol=1)
plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7))
ggsave("~/Documents/cohort/ESCC/for_NEJM/cna_16q21.pdf",
       plot=plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7)),width=10)
# ======

# Each patient
# ======
ns <- cnv_cohort_region_chr12["chr12_24980841_25799873",] %>>% 
  (.[which(pi$level == "primary")]) %>>% (.[which(. > 2)]) %>>% names
sapply(ns,function(x) {
  y <- as.character(pi[which(pi$姓名 == pi[which(pi$病人姓名 == x),1]),3])
  cnv_cohort_region_chr12["chr12_24980841_25799873",y]
}) %>>% (.[which(sapply(.,length) > 1)])
# ======

# New Fig. 3A
# ======
# Make a data.frame named "temp" by the codes of single-cell CNV heatmap
intersect(cna_common_genes,rownames(temp)) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_genes_NoBulkSNV_INTERSECT_cna_common_genes.csv")
intersect(cna_common_genes,rownames(temp)) %>>% 
  (bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")$ENTREZID) %>>%
  enrichGO(OrgDb="org.Hs.eg.db",ont="CC",readable=T,pvalueCutoff=1,qvalueCutoff=1,minGSSize=1) %>>% 
  write.csv(file="~/Documents/cohort/ESCC/for_NEJM/common_genes_NoBulkSNV_INTERSECT_cna_common_genes_GO_CC.csv")
c(cna_common[1:346,1],cna_common[347:nrow(cna_common),2]) %>>% 
  paste(collapse=",") %>>% (str_split(.,",")[[1]]) %>>% unique %>>% length
cosmic_census[match(intersect(cna_common_genes,rownames(temp)),
                    cosmic_census$Gene.Symbol),
              c("Gene.Symbol","Role.in.Cancer")] %>>% 
  (.[!is.na(.[,2]) & .[,2] != "",]) %>>% 
  write.csv(
    file="~/Documents/cohort/ESCC/for_NEJM/common_genes_NoBulkSNV_cna_oncogene_TGS_intersection.csv",
    row.names=F,quote=F)
# =======

# CNA plot for new regions 
# ======
cnv_cohort_region_chr1 <- grep("^chr1",rownames(cnv_cohort_region)) %>>% 
  (cnv_cohort_region[.,]) 
# Remove out-regions
# ---
cnv_cohort_region_chr1 <- 
  cbind(cnv_cohort_region_chr1,
        start=(str_split(rownames(cnv_cohort_region_chr1),"_") %>>% 
          sapply(function(x) x[2]) %>>% as.integer)) %>>%
  cbind(end=(str_split(rownames(cnv_cohort_region_chr1),"_") %>>% 
               sapply(function(x) x[3]) %>>% as.integer))
# 249250621 - 10000 = 249240621
out <- which((cnv_cohort_region_chr1$start < 10000 & 
                cnv_cohort_region_chr1$end < 10000) | 
               (cnv_cohort_region_chr1$start > 249240621 &
                  cnv_cohort_region_chr1$end > 249240621))
if (length(out) > 1)
  cnv_cohort_region_chr1 <- 
  cnv_cohort_region_chr1[-out,] 
cnv_cohort_region_chr1$start <- cnv_cohort_region_chr1$start %>>%
  {.[. < 10000] <- 10000;.}
cnv_cohort_region_chr1$end <- cnv_cohort_region_chr1$end %>>%
{.[. > 249240621] <- 249240621;.}
rownames(cnv_cohort_region_chr1) <- paste("chr1",
                                          cnv_cohort_region_chr1$start,
                                          cnv_cohort_region_chr1$end,
                                          sep="_")
cnv_cohort_region_chr1 <- cnv_cohort_region_chr1[,1:153]
# ---
sp <- split(1:249250621,ceiling(seq_along(1:249250621)/2492507))
sp <- sapply(sp,function(x) return(c(x[1],x[length(x)]))) %>>% t %>>% 
  as.data.frame
colnames(sp) <- c("start","end")
average_copyNumber <- apply(sp,1,function(x) {
  i <- cnv_23[which(cnv_23$Chr == "chr1"),] %>>%
    (.[which((.$Start >= x[1] & .$End <= x[2]) |
               (.$Start < x[1] & .$End >= x[1]) |
               (.$Start <= x[2] & .$End > x[2])),])
  j <- as.matrix(i[,4:33]) %>>% mean
  k <- as.matrix(i[,34:48]) %>>% mean
  a <- sapply(rownames(cnv_cohort_region_chr1),
              function(x) str_split(x,"_")[[1]][2:3] %>>% as.integer) %>>% t
  a <- cnv_cohort_region_chr1[which((a[,1] >= x[1] & a[,2] <= x[2]) |
                                      (a[,1] < x[1] & a[,2] >= x[1]) |
                                      (a[,1] <= x[2] & a[,2] > x[2])),]
  b <- as.matrix(a[,as.character(subset(pi,level == "primary")$病人姓名)]) %>>% 
    mean
  c <- as.matrix(a[,as.character(subset(pi,level == "therapy")$病人姓名)]) %>>% 
    mean
  d <- a[,as.character(subset(pi,level == "therapy_after")$病人姓名)] %>>% 
    as.matrix %>>% mean
  e <- as.matrix(a[,as.character(subset(pi,level == "relapse")$病人姓名)]) %>>% 
    mean
  c(j,k,b,c,d,e)
}) %>>% t
average_copyNumber <- as.data.frame(average_copyNumber)
colnames(average_copyNumber) <- c("Primary_sc","PostTherapy_sc","Primary",
                                  "Therapy","PostTherapy","Relapse")
average_copyNumber <- apply(average_copyNumber,1,function(x) {
  x[is.nan(x)] <- 2
  return(x)
}) %>>% t %>>% as.data.frame
sp <- cbind(sp,average_copyNumber)
group <- function(x) {
  a <- rep(3,length(which(sp[,x] > 2)))
  names(a) <- which(sp[,x] > 2)
  b <- rep(2,length(which(sp[,x] == 2)))
  names(b) <- which(sp[,x] == 2)
  c <- rep(1,length(which(sp[,x] < 2)))
  names(c) <- which(sp[,x] < 2)
  y <- c(a,b,c)
  y <- y[order(as.integer(names(y)))] %>>% as.character
  y[y == "1"] <- "Deletion"
  y[y == "2"] <- "Normal"
  y[y == "3"] <- "Amplification"
  y <- as.factor(y)
  if (!"Normal" %in% levels(y))
    levels(y) <- c(levels(y),"Normal")
  if (!"Amplification" %in% levels(y))
    levels(y) <- c(levels(y),"Amplification")
  if (!"Deletion" %in% levels(y))
    levels(y) <- c(levels(y),"Deletion")
  y <- relevel(y,"Amplification")
  y <- relevel(y,"Normal")
  return(y)
}
library(ggthemes)
color_pal <- pal_lancet()(3)
names(color_pal) <- c("Normal","Amplification","Deletion")
theme_set(theme_bw(base_family="serif"))
g1 <- ggplot(cbind(sp,group=group("Primary_sc"))) + 
  geom_point(aes(start,Primary_sc,color=group)) + 
  geom_vline(xintercept=147415551,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,249250621),labels=NULL) + 
  annotate("text",x=100000000,y=3.25,label="chr1:147415551-15107969",
           family="serif",size=5,color="cornflowerblue",fontface="italic") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g2 <- ggplot(cbind(sp,group=group("PostTherapy_sc"))) + 
  geom_point(aes(start,PostTherapy_sc,color=group)) + 
  geom_vline(xintercept=147415551,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,249250621),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g3 <- ggplot(cbind(sp,group=group("Primary"))) + 
  geom_point(aes(start,Primary,color=group)) + 
  geom_vline(xintercept=147415551,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,249250621),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g4 <- ggplot(cbind(sp,group=group("Therapy"))) + 
  geom_point(aes(start,Therapy,color=group)) + 
  geom_vline(xintercept=147415551,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,249250621),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g5 <- ggplot(cbind(sp,group=group("PostTherapy"))) + 
  geom_point(aes(start,PostTherapy,color=group)) + 
  geom_vline(xintercept=147415551,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,249250621),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g6 <- ggplot(cbind(sp,group=group("Relapse"))) + 
  geom_point(aes(start,Relapse,color=group)) + 
  geom_vline(xintercept=147415551,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x="Chromosome 1") + 
  scale_x_continuous(breaks=c(1,249250621),labels=c(1,249250621)) + 
  guides(color=guide_legend(title.theme=element_text(face="bold",size=20,
                                                     angle=0,family="serif"),
                            label.theme=element_text(size=18,angle=0,
                                                     family="serif"),
                            title="Copy Number")) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
legend <- get_legend(g6)
prow <- plot_grid(g1 + theme(legend.position="none"),
                  g2 + theme(legend.position="none"),
                  g3 + theme(legend.position="none"),
                  g4 + theme(legend.position="none"),
                  g5 + theme(legend.position="none"),
                  g6 + theme(legend.position="none"),align="v",ncol=1)
plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7))
ggsave("~/Documents/cohort/ESCC/for_NEJM/cna_chr1_147415551_15107969.pdf",
       plot=plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7)),width=10,
       height=8)

# ---
# Almost deletion in single-cell T_T
cnv_cohort_region_chr3 <- grep("^chr3",rownames(cnv_cohort_region)) %>>% 
  (cnv_cohort_region[.,]) 
# Remove out-regions
# ---
cnv_cohort_region_chr3 <- 
  cbind(cnv_cohort_region_chr3,
        start=(str_split(rownames(cnv_cohort_region_chr3),"_") %>>% 
                 sapply(function(x) x[2]) %>>% as.integer)) %>>%
  cbind(end=(str_split(rownames(cnv_cohort_region_chr3),"_") %>>% 
               sapply(function(x) x[3]) %>>% as.integer))
# 198022430 - 60000 = 197962430
out <- which((cnv_cohort_region_chr3$start < 60000 & 
                cnv_cohort_region_chr3$end < 60000) | 
               (cnv_cohort_region_chr3$start > 197962430 &
                  cnv_cohort_region_chr3$end > 197962430))
if (length(out) > 1)
  cnv_cohort_region_chr3 <- 
  cnv_cohort_region_chr3[-out,] 
cnv_cohort_region_chr3$start <- cnv_cohort_region_chr3$start %>>%
{.[. < 60000] <- 60000;.}
cnv_cohort_region_chr3$end <- cnv_cohort_region_chr3$end %>>%
{.[. > 197962430] <- 197962430;.}
rownames(cnv_cohort_region_chr3) <- paste("chr3",
                                          cnv_cohort_region_chr3$start,
                                          cnv_cohort_region_chr3$end,
                                          sep="_")
cnv_cohort_region_chr3 <- cnv_cohort_region_chr3[,1:153]
# ---
sp <- split(1:198022430,ceiling(seq_along(1:198022430)/1980225))
sp <- sapply(sp,function(x) return(c(x[1],x[length(x)]))) %>>% t %>>% 
  as.data.frame
colnames(sp) <- c("start","end")
average_copyNumber <- apply(sp,1,function(x) {
  i <- cnv_23[which(cnv_23$Chr == "chr3"),] %>>%
    (.[which((.$Start >= x[1] & .$End <= x[2]) |
               (.$Start < x[1] & .$End >= x[1]) |
               (.$Start <= x[2] & .$End > x[2])),])
  j <- as.matrix(i[,4:33]) %>>% mean
  k <- as.matrix(i[,34:48]) %>>% mean
  a <- sapply(rownames(cnv_cohort_region_chr3),
              function(x) str_split(x,"_")[[1]][2:3] %>>% as.integer) %>>% t
  a <- cnv_cohort_region_chr3[which((a[,1] >= x[1] & a[,2] <= x[2]) |
                                      (a[,1] < x[1] & a[,2] >= x[1]) |
                                      (a[,1] <= x[2] & a[,2] > x[2])),]
  b <- as.matrix(a[,as.character(subset(pi,level == "primary")$病人姓名)]) %>>% 
    mean
  c <- as.matrix(a[,as.character(subset(pi,level == "therapy")$病人姓名)]) %>>% 
    mean
  d <- a[,as.character(subset(pi,level == "therapy_after")$病人姓名)] %>>% 
    as.matrix %>>% mean
  e <- as.matrix(a[,as.character(subset(pi,level == "relapse")$病人姓名)]) %>>% 
    mean
  c(j,k,b,c,d,e)
}) %>>% t
average_copyNumber <- as.data.frame(average_copyNumber)
colnames(average_copyNumber) <- c("Primary_sc","PostTherapy_sc","Primary",
                                  "Therapy","PostTherapy","Relapse")
average_copyNumber <- apply(average_copyNumber,1,function(x) {
  x[is.nan(x)] <- 2
  return(x)
}) %>>% t %>>% as.data.frame
sp <- cbind(sp,average_copyNumber)
group <- function(x) {
  a <- rep(3,length(which(sp[,x] > 2)))
  names(a) <- which(sp[,x] > 2)
  b <- rep(2,length(which(sp[,x] == 2)))
  names(b) <- which(sp[,x] == 2)
  c <- rep(1,length(which(sp[,x] < 2)))
  names(c) <- which(sp[,x] < 2)
  y <- c(a,b,c)
  y <- y[order(as.integer(names(y)))] %>>% as.character
  y[y == "1"] <- "Deletion"
  y[y == "2"] <- "Normal"
  y[y == "3"] <- "Amplification"
  y <- as.factor(y)
  if (!"Normal" %in% levels(y))
    levels(y) <- c(levels(y),"Normal")
  if (!"Amplification" %in% levels(y))
    levels(y) <- c(levels(y),"Amplification")
  if (!"Deletion" %in% levels(y))
    levels(y) <- c(levels(y),"Deletion")
  y <- relevel(y,"Amplification")
  y <- relevel(y,"Normal")
  return(y)
}
library(ggthemes)
color_pal <- pal_lancet()(3)
names(color_pal) <- c("Normal","Amplification","Deletion")
theme_set(theme_bw(base_family="serif"))
# For showing "Good" picture
sp_11 <- sp[11,]
sp[11,3:4] <- sp[11,3:4] + 1
g1 <- ggplot(cbind(sp,group=group("Primary_sc"))) + 
  geom_point(aes(start,Primary_sc,color=group)) + 
  geom_vline(xintercept=20044101,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,198022430),labels=NULL) + 
  annotate("text",x=52000000,y=3.25,label="chr3:20044101-31703623",
           family="serif",size=5,color="cornflowerblue",fontface="italic") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g2 <- ggplot(cbind(sp,group=group("PostTherapy_sc"))) + 
  geom_point(aes(start,PostTherapy_sc,color=group)) + 
  geom_vline(xintercept=20044101,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,198022430),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g3 <- ggplot(cbind(sp,group=group("Primary"))) + 
  geom_point(aes(start,Primary,color=group)) + 
  geom_vline(xintercept=20044101,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,198022430),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g4 <- ggplot(cbind(sp,group=group("Therapy"))) + 
  geom_point(aes(start,Therapy,color=group)) + 
  geom_vline(xintercept=20044101,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,198022430),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g5 <- ggplot(cbind(sp,group=group("PostTherapy"))) + 
  geom_point(aes(start,PostTherapy,color=group)) + 
  geom_vline(xintercept=20044101,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,198022430),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g6 <- ggplot(cbind(sp,group=group("Relapse"))) + 
  geom_point(aes(start,Relapse,color=group)) + 
  geom_vline(xintercept=20044101,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x="Chromosome 3") + 
  scale_x_continuous(breaks=c(1,198022430),labels=c(1,198022430)) + 
  guides(color=guide_legend(title.theme=element_text(face="bold",size=20,
                                                     angle=0,family="serif"),
                            label.theme=element_text(size=18,angle=0,
                                                     family="serif"),
                            title="Copy Number")) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
legend <- get_legend(g6)
prow <- plot_grid(g1 + theme(legend.position="none"),
                  g2 + theme(legend.position="none"),
                  g3 + theme(legend.position="none"),
                  g4 + theme(legend.position="none"),
                  g5 + theme(legend.position="none"),
                  g6 + theme(legend.position="none"),align="v",ncol=1)
plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7))
ggsave("~/Documents/cohort/ESCC/for_NEJM/cna_chr3_20044101_31703623.pdf",
       plot=plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7)),width=10,
       height=8)

# ---
cnv_cohort_region_chr7 <- grep("^chr7",rownames(cnv_cohort_region)) %>>% 
  (cnv_cohort_region[.,]) 
# Remove out-regions
# ---
cnv_cohort_region_chr7 <- 
  cbind(cnv_cohort_region_chr7,
        start=(str_split(rownames(cnv_cohort_region_chr7),"_") %>>% 
                 sapply(function(x) x[2]) %>>% as.integer)) %>>%
  cbind(end=(str_split(rownames(cnv_cohort_region_chr7),"_") %>>% 
               sapply(function(x) x[3]) %>>% as.integer))
# 159138663 - 10000 = 159128663
out <- which((cnv_cohort_region_chr7$start < 10000 & 
                cnv_cohort_region_chr7$end < 10000) | 
               (cnv_cohort_region_chr7$start > 159128663 &
                  cnv_cohort_region_chr7$end > 159128663))
if (length(out) > 1)
  cnv_cohort_region_chr7 <- 
  cnv_cohort_region_chr7[-out,] 
cnv_cohort_region_chr7$start <- cnv_cohort_region_chr7$start %>>%
{.[. < 10000] <- 10000;.}
cnv_cohort_region_chr7$end <- cnv_cohort_region_chr7$end %>>%
{.[. > 159128663] <- 159128663;.}
rownames(cnv_cohort_region_chr7) <- paste("chr7",
                                          cnv_cohort_region_chr7$start,
                                          cnv_cohort_region_chr7$end,
                                          sep="_")
cnv_cohort_region_chr7 <- cnv_cohort_region_chr7[,1:153]
# ---
sp <- split(1:159138663,ceiling(seq_along(1:159138663)/1591387))
sp <- sapply(sp,function(x) return(c(x[1],x[length(x)]))) %>>% t %>>% 
  as.data.frame
colnames(sp) <- c("start","end")
average_copyNumber <- apply(sp,1,function(x) {
  i <- cnv_23[which(cnv_23$Chr == "chr7"),] %>>%
    (.[which((.$Start >= x[1] & .$End <= x[2]) |
               (.$Start < x[1] & .$End >= x[1]) |
               (.$Start <= x[2] & .$End > x[2])),])
  j <- as.matrix(i[,4:33]) %>>% mean
  k <- as.matrix(i[,34:48]) %>>% mean
  a <- sapply(rownames(cnv_cohort_region_chr7),
              function(x) str_split(x,"_")[[1]][2:3] %>>% as.integer) %>>% t
  a <- cnv_cohort_region_chr7[which((a[,1] >= x[1] & a[,2] <= x[2]) |
                                      (a[,1] < x[1] & a[,2] >= x[1]) |
                                      (a[,1] <= x[2] & a[,2] > x[2])),]
  b <- as.matrix(a[,as.character(subset(pi,level == "primary")$病人姓名)]) %>>% 
    mean
  c <- as.matrix(a[,as.character(subset(pi,level == "therapy")$病人姓名)]) %>>% 
    mean
  d <- a[,as.character(subset(pi,level == "therapy_after")$病人姓名)] %>>% 
    as.matrix %>>% mean
  e <- as.matrix(a[,as.character(subset(pi,level == "relapse")$病人姓名)]) %>>% 
    mean
  c(j,k,b,c,d,e)
}) %>>% t
average_copyNumber <- as.data.frame(average_copyNumber)
colnames(average_copyNumber) <- c("Primary_sc","PostTherapy_sc","Primary",
                                  "Therapy","PostTherapy","Relapse")
average_copyNumber <- apply(average_copyNumber,1,function(x) {
  x[is.nan(x)] <- 2
  return(x)
}) %>>% t %>>% as.data.frame
sp <- cbind(sp,average_copyNumber)
group <- function(x) {
  a <- rep(3,length(which(sp[,x] > 2)))
  names(a) <- which(sp[,x] > 2)
  b <- rep(2,length(which(sp[,x] == 2)))
  names(b) <- which(sp[,x] == 2)
  c <- rep(1,length(which(sp[,x] < 2)))
  names(c) <- which(sp[,x] < 2)
  y <- c(a,b,c)
  y <- y[order(as.integer(names(y)))] %>>% as.character
  y[y == "1"] <- "Deletion"
  y[y == "2"] <- "Normal"
  y[y == "3"] <- "Amplification"
  y <- as.factor(y)
  if (!"Normal" %in% levels(y))
    levels(y) <- c(levels(y),"Normal")
  if (!"Amplification" %in% levels(y))
    levels(y) <- c(levels(y),"Amplification")
  if (!"Deletion" %in% levels(y))
    levels(y) <- c(levels(y),"Deletion")
  y <- relevel(y,"Amplification")
  y <- relevel(y,"Normal")
  return(y)
}
library(ggthemes)
color_pal <- pal_lancet()(3)
names(color_pal) <- c("Normal","Amplification","Deletion")
theme_set(theme_bw(base_family="serif"))
g1 <- ggplot(cbind(sp,group=group("Primary_sc"))) + 
  geom_point(aes(start,Primary_sc,color=group)) + 
  geom_vline(xintercept=85934899,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,159138663),labels=NULL) + 
  annotate("text",x=58000000,y=5,label="chr7:84324513-102463262",
           family="serif",size=5,color="cornflowerblue",fontface="italic") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g2 <- ggplot(cbind(sp,group=group("PostTherapy_sc"))) + 
  geom_point(aes(start,PostTherapy_sc,color=group)) + 
  geom_vline(xintercept=85934899,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,159138663),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g3 <- ggplot(cbind(sp,group=group("Primary"))) + 
  geom_point(aes(start,Primary,color=group)) + 
  geom_vline(xintercept=85934899,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,159138663),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g4 <- ggplot(cbind(sp,group=group("Therapy"))) + 
  geom_point(aes(start,Therapy,color=group)) + 
  geom_vline(xintercept=85934899,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,159138663),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g5 <- ggplot(cbind(sp,group=group("PostTherapy"))) + 
  geom_point(aes(start,PostTherapy,color=group)) + 
  geom_vline(xintercept=85934899,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,159138663),labels=NULL) + 
  guides(color=guide_legend(title.theme=element_text(face="bold",size=20,
                                                     angle=0,family="serif"),
                            label.theme=element_text(size=18,angle=0,
                                                     family="serif"),
                            title="Copy Number")) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g6 <- ggplot(cbind(sp,group=group("Relapse"))) + 
  geom_point(aes(start,Relapse,color=group)) + 
  geom_vline(xintercept=85934899,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x="Chromosome 7") + 
  scale_x_continuous(breaks=c(1,159138663),labels=c(1,159138663)) + 
    theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
legend <- get_legend(g5)
prow <- plot_grid(g1 + theme(legend.position="none"),
                  g2 + theme(legend.position="none"),
                  g3 + theme(legend.position="none"),
                  g4 + theme(legend.position="none"),
                  g5 + theme(legend.position="none"),
                  g6 + theme(legend.position="none"),align="v",ncol=1)
plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7))
ggsave("~/Documents/cohort/ESCC/for_NEJM/cna_chr7_84324513_102463262.pdf",
       plot=plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7)),width=10,
       height=8)

# ---
cnv_cohort_region_chr12 <- grep("^chr12",rownames(cnv_cohort_region)) %>>% 
  (cnv_cohort_region[.,]) 
# Remove out-regions
# ---
cnv_cohort_region_chr12 <- 
  cbind(cnv_cohort_region_chr12,
        start=(str_split(rownames(cnv_cohort_region_chr12),"_") %>>% 
                 sapply(function(x) x[2]) %>>% as.integer)) %>>%
  cbind(end=(str_split(rownames(cnv_cohort_region_chr12),"_") %>>% 
               sapply(function(x) x[3]) %>>% as.integer))
# 133851895 - 60000 = 133791895
out <- which((cnv_cohort_region_chr12$start < 60000 & 
                cnv_cohort_region_chr12$end < 60000) | 
               (cnv_cohort_region_chr12$start > 133791895 &
                  cnv_cohort_region_chr12$end > 133791895))
if (length(out) > 1)
  cnv_cohort_region_chr12 <- 
  cnv_cohort_region_chr12[-out,] 
cnv_cohort_region_chr12$start <- cnv_cohort_region_chr12$start %>>%
{.[. < 60000] <- 60000;.}
cnv_cohort_region_chr12$end <- cnv_cohort_region_chr12$end %>>%
{.[. > 133791895] <- 133791895;.}
rownames(cnv_cohort_region_chr12) <- paste("chr12",
                                          cnv_cohort_region_chr12$start,
                                          cnv_cohort_region_chr12$end,
                                          sep="_")
cnv_cohort_region_chr12 <- cnv_cohort_region_chr12[,1:153]
# ---
sp <- split(1:133851895,ceiling(seq_along(1:133851895)/1338519))
sp <- sapply(sp,function(x) return(c(x[1],x[length(x)]))) %>>% t %>>% 
  as.data.frame
colnames(sp) <- c("start","end")
average_copyNumber <- apply(sp,1,function(x) {
  i <- cnv_23[which(cnv_23$Chr == "chr12"),] %>>%
    (.[which((.$Start >= x[1] & .$End <= x[2]) |
               (.$Start < x[1] & .$End >= x[1]) |
               (.$Start <= x[2] & .$End > x[2])),])
  j <- as.matrix(i[,4:33]) %>>% mean
  k <- as.matrix(i[,34:48]) %>>% mean
  a <- sapply(rownames(cnv_cohort_region_chr12),
              function(x) str_split(x,"_")[[1]][2:3] %>>% as.integer) %>>% t
  a <- cnv_cohort_region_chr12[which((a[,1] >= x[1] & a[,2] <= x[2]) |
                                      (a[,1] < x[1] & a[,2] >= x[1]) |
                                      (a[,1] <= x[2] & a[,2] > x[2])),]
  b <- as.matrix(a[,as.character(subset(pi,level == "primary")$病人姓名)]) %>>% 
    mean
  c <- as.matrix(a[,as.character(subset(pi,level == "therapy")$病人姓名)]) %>>% 
    mean
  d <- a[,as.character(subset(pi,level == "therapy_after")$病人姓名)] %>>% 
    as.matrix %>>% mean
  e <- as.matrix(a[,as.character(subset(pi,level == "relapse")$病人姓名)]) %>>% 
    mean
  c(j,k,b,c,d,e)
}) %>>% t
average_copyNumber <- as.data.frame(average_copyNumber)
colnames(average_copyNumber) <- c("Primary_sc","PostTherapy_sc","Primary",
                                  "Therapy","PostTherapy","Relapse")
average_copyNumber <- apply(average_copyNumber,1,function(x) {
  x[is.nan(x)] <- 2
  return(x)
}) %>>% t %>>% as.data.frame
sp <- cbind(sp,average_copyNumber)
group <- function(x) {
  a <- rep(3,length(which(sp[,x] > 2)))
  names(a) <- which(sp[,x] > 2)
  b <- rep(2,length(which(sp[,x] == 2)))
  names(b) <- which(sp[,x] == 2)
  c <- rep(1,length(which(sp[,x] < 2)))
  names(c) <- which(sp[,x] < 2)
  y <- c(a,b,c)
  y <- y[order(as.integer(names(y)))] %>>% as.character
  y[y == "1"] <- "Deletion"
  y[y == "2"] <- "Normal"
  y[y == "3"] <- "Amplification"
  y <- as.factor(y)
  if (!"Normal" %in% levels(y))
    levels(y) <- c(levels(y),"Normal")
  if (!"Amplification" %in% levels(y))
    levels(y) <- c(levels(y),"Amplification")
  if (!"Deletion" %in% levels(y))
    levels(y) <- c(levels(y),"Deletion")
  y <- relevel(y,"Amplification")
  y <- relevel(y,"Normal")
  return(y)
}
library(ggthemes)
color_pal <- pal_lancet()(3)
names(color_pal) <- c("Normal","Amplification","Deletion")
theme_set(theme_bw(base_family="serif"))
g1 <- ggplot(cbind(sp,group=group("Primary_sc"))) + 
  geom_point(aes(start,Primary_sc,color=group)) + 
  geom_vline(xintercept=20077786,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,133851895),labels=NULL) + 
  annotate("text",x=45000000,y=6,label="chr12:19552035-23974957",
           family="serif",size=5,color="cornflowerblue",fontface="italic") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g2 <- ggplot(cbind(sp,group=group("PostTherapy_sc"))) + 
  geom_point(aes(start,PostTherapy_sc,color=group)) + 
  geom_vline(xintercept=20077786,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,133851895),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g3 <- ggplot(cbind(sp,group=group("Primary"))) + 
  geom_point(aes(start,Primary,color=group)) + 
  geom_vline(xintercept=20077786,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,133851895),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g4 <- ggplot(cbind(sp,group=group("Therapy"))) + 
  geom_point(aes(start,Therapy,color=group)) + 
  geom_vline(xintercept=20077786,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,133851895),labels=NULL) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g5 <- ggplot(cbind(sp,group=group("PostTherapy"))) + 
  geom_point(aes(start,PostTherapy,color=group)) + 
  geom_vline(xintercept=20077786,linetype=3) + 
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x=NULL) + scale_x_continuous(breaks=c(1,133851895),labels=NULL) + 
    theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
g6 <- ggplot(cbind(sp,group=group("Relapse"))) + 
  geom_point(aes(start,Relapse,color=group)) + 
  geom_vline(xintercept=20077786,linetype=3) +
  geom_hline(yintercept=2,linetype=2) + scale_color_manual(values=color_pal) + 
  labs(x="Chromosome 12") + 
  scale_x_continuous(breaks=c(1,133851895),labels=c(1,133851895)) + 
  guides(color=guide_legend(title.theme=element_text(face="bold",size=20,
                                                     angle=0,family="serif"),
                            label.theme=element_text(size=18,angle=0,
                                                     family="serif"),
                            title="Copy Number")) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
legend <- get_legend(g6)
prow <- plot_grid(g1 + theme(legend.position="none"),
                  g2 + theme(legend.position="none"),
                  g3 + theme(legend.position="none"),
                  g4 + theme(legend.position="none"),
                  g5 + theme(legend.position="none"),
                  g6 + theme(legend.position="none"),align="v",ncol=1)
plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7))
ggsave("~/Documents/cohort/ESCC/for_NEJM/cna_chr12_19552035_23974957.pdf",
       plot=plot_grid(prow,legend,ncol=2,rel_widths=c(3,0.7)),width=10,
       height=8)
# ======

# CNA table for new regions
# ======
i <- cnv_23[which(cnv_23$Chr == "chr1"),] %>>% 
  (.[which((.$Start >= 147415551 & .$End <= 151079698) |
             (.$Start < 147415551 & .$End >= 147415551) |
             (.$Start <= 151079698 & .$End > 151079698)),]) %>>% 
  apply(1,function(x) {a <- x[4:33] %>>% as.integer %>>% mean;b <- x[34:48] %>>%
    as.integer %>>% mean;c(a,b)}) %>>% t %>>% (.[.[,1] > 2 & .[,2] > 2,]) %>>% 
  rownames %>>% (cnv_23[.,52:53]) %>>% {
    a <- paste(rownames(.),collapse=",")
    b <- paste(str_split(.$an23.cytoBand,",") %>>% unlist %>>% unique,
               collapse=",")
    c <- paste(str_split(.$an23.Genes,",") %>>% unlist %>>% unique,collapse=",")
    d <- str_split(.$an23.Genes,",") %>>% unlist %>>% unique %>>% 
      (cosmic_census[match(.,cosmic_census$Gene.Symbol),
                     c("Gene.Symbol","Role.in.Cancer")]) %>>% 
      (.[!is.na(.[,2]) & .[,2] != "",]) %>>% apply(1,paste,collapse=":") %>>% 
      paste(collapse=",")
    c(a,b,c,d)
  }
j <- cnv_23[which(cnv_23$Chr == "chr3"),] %>>% 
  (.[which((.$Start >= 20044101 & .$End <= 31703623) |
             (.$Start < 20044101 & .$End >= 20044101) |
             (.$Start <= 31703623 & .$End > 31703623)),]) %>>% 
  apply(1,function(x) {a <- x[4:33] %>>% as.integer %>>% mean;b <- x[34:48] %>>%
    as.integer %>>% mean;c(a,b)}) %>>% t %>>% (.[.[,1] > 2 & .[,2] > 2,]) %>>% 
  rownames %>>% (cnv_23[.,52:53]) %>>% {
    a <- paste(rownames(.),collapse=",")
    b <- paste(str_split(.$an23.cytoBand,",") %>>% unlist %>>% unique,
               collapse=",")
    c <- paste(str_split(.$an23.Genes,",") %>>% unlist %>>% unique,collapse=",")
    d <- str_split(.$an23.Genes,",") %>>% unlist %>>% unique %>>% 
      (cosmic_census[match(.,cosmic_census$Gene.Symbol),
                     c("Gene.Symbol","Role.in.Cancer")]) %>>% 
      (.[!is.na(.[,2]) & .[,2] != "",]) %>>% apply(1,paste,collapse=":") %>>% 
      paste(collapse=",")
    c(a,b,c,d)
  }
k <- cnv_23[which(cnv_23$Chr == "chr7"),] %>>% 
  (.[which((.$Start >= 84324513 & .$End <= 102463262) |
             (.$Start < 84324513 & .$End >= 84324513) |
             (.$Start <= 102463262 & .$End > 102463262)),]) %>>% 
  apply(1,function(x) {a <- x[4:33] %>>% as.integer %>>% mean;b <- x[34:48] %>>%
    as.integer %>>% mean;c(a,b)}) %>>% t %>>% (.[.[,1] > 2 & .[,2] > 2,]) %>>% 
  rownames %>>% (cnv_23[.,52:53]) %>>% {
    a <- paste(rownames(.),collapse=",")
    b <- paste(str_split(.$an23.cytoBand,",") %>>% unlist %>>% unique,
               collapse=",")
    c <- paste(str_split(.$an23.Genes,",") %>>% unlist %>>% unique,collapse=",")
    d <- str_split(.$an23.Genes,",") %>>% unlist %>>% unique %>>% 
      (cosmic_census[match(.,cosmic_census$Gene.Symbol),
                     c("Gene.Symbol","Role.in.Cancer")]) %>>% 
      (.[!is.na(.[,2]) & .[,2] != "",]) %>>% apply(1,paste,collapse=":") %>>% 
      paste(collapse=",")
    c(a,b,c,d)
  }
l <- cnv_23[which(cnv_23$Chr == "chr12"),] %>>% 
  (.[which((.$Start >= 19552035 & .$End <= 23974957) |
             (.$Start < 19552035 & .$End >= 19552035) |
             (.$Start <= 23974957 & .$End > 23974957)),]) %>>% 
  apply(1,function(x) {a <- x[4:33] %>>% as.integer %>>% mean;b <- x[34:48] %>>%
    as.integer %>>% mean;c(a,b)}) %>>% t %>>% (.[.[,1] > 2 & .[,2] > 2,]) %>>% 
  rownames %>>% (cnv_23[.,52:53]) %>>% {
    a <- paste(rownames(.),collapse=",")
    b <- paste(str_split(.$an23.cytoBand,",") %>>% unlist %>>% unique,
               collapse=",")
    c <- paste(str_split(.$an23.Genes,",") %>>% unlist %>>% unique,collapse=",")
    d <- str_split(.$an23.Genes,",") %>>% unlist %>>% unique %>>% 
      (cosmic_census[match(.,cosmic_census$Gene.Symbol),
                     c("Gene.Symbol","Role.in.Cancer")]) %>>% 
      (.[!is.na(.[,2]) & .[,2] != "",]) %>>% apply(1,paste,collapse=":") %>>% 
      paste(collapse=",")
    c(a,b,c,d)
  }
a <- rbind(i,j,k,l)
colnames(a) <- c("Position","CytoBand","Genes","Oncogene/TSG")
rownames(a) <- c("chr1_147415551_151079698","chr3_20044101_31703623",
                 "chr7_84324513_102463262","chr12_19552035_23974957")
data.table::fwrite(as.data.frame(a),
                   file="~/Documents/cohort/ESCC/for_NEJM/new_region.tsv",
                   row.names=T,sep="\t")
