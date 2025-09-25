# load libraries 
suppressMessages(library(edgeR))
suppressMessages(library(data.table))
suppressMessages(library(DESeq2))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(VennDiagram))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(viridis))
suppressMessages(library(pals))
suppressMessages(library(philentropy))
suppressMessages(library(IntClust))
suppressMessages(library(pheatmap))
suppressMessages(library(assertr))
suppressMessages(library(remotes))
suppressMessages(library(GeneOverlap))
suppressMessages(library(stringr))
suppressMessages(library(ggrepel))
suppressMessages(library(msigdbr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Mm.eg.db)) 
suppressMessages(library(dplyr))
suppressMessages(library(proxy))



my_filecount = 'all.counts.gz'

metadata = read.xlsx("Metadata.xlsx")
rownames(metadata) <- metadata$Sample.Name

annotation <- c('GeneID','Chr','Start','End','Strand','Length')
fCounts <- read.delim(file=my_filecount, header=TRUE, check.names = F)

fCountsData <- fCounts[
  , 
  -which(
    tolower(names(fCounts))
    %in% 
      tolower(annotation))]

fCountsAnnotation <- fCounts[
  , 
  which(
    tolower(names(fCounts))
    %in% 
      tolower(annotation))]

geneidColname <- 'Geneid'
geneidIdx <- which(tolower(annotation) %in% tolower(geneidColname))
rownames(fCountsData) <- fCounts[[geneidIdx]]

fCountsData <- fCountsData[,match(metadata[,1], colnames(fCountsData))] 

fCountsData_d = fCountsData[,row.names(metadata)]

# all_counts
y <- DGEList(counts=fCountsData_d, genes = fCountsAnnotation)
SAVE_variable[[variable2save_names[1]]] <- as.data.frame(y$counts)

Nreplica = 3
# expGENES_counts
keep <- rowSums(cpm(y)>1)>=Nreplica

y <- DGEList(counts=fCountsData_d, genes = fCountsAnnotation)
# rpkm
fCountsRPKM = rpkm(y, log=T, gene.length =y$genes$Length)
# filter for expression
keep <- rowSums(cpm(y)>1)>=Nreplica
yf <- y[keep,]



SAVE_variable <- list()

variable2save_names <- c('all_counts', 'expGenes_counts','expGenes_LogCPM', 'expGenes_LogRPKM', 'expGenes_RPKM', 'allGenes_LogRPKM')

SAVE_variable[[variable2save_names[1]]] <- as.data.frame(y$counts)

SAVE_variable[[variable2save_names[2]]] <- as.data.frame(yf$counts)

SAVE_variable[[variable2save_names[3]]] <- as.data.frame(cpm(yf, log=T))


SAVE_variable[[variable2save_names[4]]] <- as.data.frame(rpkm(yf, log=T, gene.length =yf$genes$Length))


SAVE_variable[[variable2save_names[5]]] <- as.data.frame(rpkm(yf, log=F, gene.length =yf$genes$Length))


SAVE_variable[[variable2save_names[6]]] <- as.data.frame(rpkm(y, log=T, gene.length =y$genes$Length))

filename_xls <- "COUNTS_CPM_RPKM_PM_PH.xlsx"

write.xlsx(SAVE_variable,
           file = filename_xls, 
           rowNames = T,
           asTable = T, 
           sheetName =variable2save_names)



N=500
vary <- apply(fCountsRPKM[keep,],1,var)
vary_s <- sort(vary, decreasing = T)
TOP_N <- names(vary_s[1:N])
yTOP <-  y[TOP_N,]
fCountsRPKMTOP <- fCountsRPKM[TOP_N,]

#### Perform PCA of top 500 most variable genes #### 

pcx = 1
pcy = 2
centering = TRUE
scaling = TRUE

pca = prcomp(t(fCountsRPKMTOP), center=centering, scale=scaling)
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
score = as.data.frame(pca$x)

xlab = paste("PC", pcx, " (",var[pcx],"%)", sep="")
ylab = paste("PC", pcy, " (",var[pcy],"%)", sep="")
cum = var[pcx]+var[pcy]
names = rownames(pca$x)

score$Sample.Name = metadata$Sample.Name
score$Condition = metadata$Condition
score$Sample = as.character(metadata$Sample)

pca <- ggplot(score, aes(x=score[,pcx], y=score[,pcy],
                         color=Condition, shape=Sample)) +
  geom_label_repel(data=score, aes(x=score[,pcx], y=score[,pcy],
                                   color=Sample, label=Sample.Name),
                   size=6, box.padding=unit(0.55, "lines"), point.padding=unit(0.55, "lines"),
                   segment.color='grey50') +
  geom_point(size=7) +
  labs(x=xlab, y=ylab, title=paste("PC", pcx, " vs PC", pcy, " scoreplot", sep="")) +
  geom_hline(yintercept=0, linetype="dashed", color="darkgrey") +
  geom_vline(xintercept=0, linetype="dashed", color="darkgrey") +
  scale_color_manual(values=c("LPS"="#FFFF66", "PH"="orange", "PM" = "pink", "SE" = "#66B2FF", "Vehicle" = "darkgreen")) +  # Manually setting colors for Condition
  theme(plot.title=element_text(color="black", size=26, face="bold.italic"),
        axis.text.x=element_text(angle=0, face="bold", color="black", size=22, hjust=.5),
        axis.title.x=element_text(face="bold", color="black", size=24),
        axis.text.y=element_text(angle=0, face="bold", color="black", size=22),
        axis.title.y=element_text(face="bold", color="black", size=24),
        legend.text=element_text(face="bold", color="black", size=18),
        legend.position="right",
        panel.background=element_rect(fill="white", colour="black", size=1, linetype="solid"))

pdf('PCA_top500rpkm_PM_PH.pdf', width = 8, height = 6)
pca
dev.off()

#### Perform heatmap on top 500 most variable genes #### 

annotation_column <- metadata[,2:(dim(metadata)[2])]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 2
cols = gg_color_hue(n)

dev.new(width = 4, height = 4)

mycolors_b <- c('red','blue'); names(mycolors_b) = levels(annotation_column$Condition)

ann_colors = list()
ann_colors = list(
  Condition = mycolors_b
)


options(repr.plot.width=12, repr.plot.height=10)

annotation_column <- metadata
annotation_column$Sample.Name <- NULL
annotation_column$Sample <- as.character(annotation_column$Sample)

ann_colors = list(
  Sample = c("1" = "forestgreen",
             "3" = "brown",
             "4" = "#999900"),
  Condition = c("PH"="orange", "PM" = "pink"))


crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         annotation_colors = ann_colors,
                         cluster_rows = T,
                         cluster_cols = T,
                         show_rownames = F,
                         fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                         display_numbers = F,
                         col=colors,
                         filename = 'Heatmap_500rpkm.pdf',
                         width = 14, height = 11)


#### Perform differential gene expression analysis #### 

f = "Condition"
seqc_pvalue = 0.01
dgeResults = list()

metadata_s <- metadata[metadata$Condition %in% c("PM", "PH"),]
metadata_s$Sample <- as.character(metadata_s$Sample)
comparison = list(Condition = c("PM","PH"))

fCountsData_s = fCountsData_d[,row.names(metadata_s)]
metadata_s$Sample <- as.character(metadata_s$Sample)

for (i in names(comparison)) {
  print(i)
  fCountsData_s = fCountsData_d[,row.names(metadata_s)]
  dds <- DESeqDataSetFromMatrix(
    countData = fCountsData_s,
    colData  = metadata_s,
    design   = as.formula('~Sample + Condition'))
  filter <- rowSums(cpm(counts(dds)) >= 1) >= Nreplica
  table(filter)
  ddsFiltered <- dds[filter,]
  dga <- DESeq(
    object = ddsFiltered,
    test = "Wald",
    fitType = "parametric",
    betaPrior = FALSE,
    minReplicatesForReplace = Inf)

  
  alpha = 0.05
  print(paste(comparison[[i]][1],"_vs_",comparison[[i]][2],sep=''))
  dgeResults.tmp <- results(dga,
                            contrast             = c(f,comparison[[i]][1],comparison[[i]][2]),
                            cooksCutoff          = Inf,
                            independentFiltering = TRUE,
                            alpha                = alpha,
                            pAdjustMethod        = "BH")
  summary(dgeResults.tmp)
  dgeResults[[i]] <- dgeResults.tmp[order(dgeResults.tmp$pvalue, decreasing = F),]
  
  #PCA
  vsd <- vst(dga, blind=FALSE)
  main.factor = "Condition"
  pcaData <- plotPCA(vsd, intgroup=c(main.factor),returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  PCA = ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
    geom_label_repel(data= pcaData, aes(PC1, PC2, color=Condition, label = name),
                     size = 6,  box.padding = unit(0.55, "lines"), point.padding = unit(0.55, "lines"),
                     segment.color = 'grey50') +
    geom_point(size=6) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(paste("PCA",i)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    scale_color_manual(values = c('salmon','grey'))
  pdf(paste('pca_DEGs.pdf'),width=8, height=6)
  print(PCA)
  dev.off()
  
  seqcUP = row.names(dgeResults[[i]])[dgeResults[[i]]$pvalue <= seqc_pvalue &
                                        !is.na(dgeResults[[i]]$padj)&
                                        dgeResults[[i]]$log2FoldChange > 1]
  
  if (length(seqcUP) > 0) {
    # print heatmap
    annotation_column <- metadata_s
    row.names(annotation_column) <- row.names(metadata_s)
    options(repr.plot.width=12, repr.plot.height=10)
    crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
    colors = crp(255)
    head(cpm(counts(dga)))
    HP <- pheatmap::pheatmap(cpm(counts(dga))[seqcUP,],
                             scale = 'row',
                             cluster_rows = T,
                             cluster_cols = T,
                             show_rownames = F,
                             cutree_cols = 2,
                             fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                             display_numbers = F,
                             col=colors,
                             filename = paste('Heatmap_seqcUP.pdf'),
                             width = 10, height = 11 )
  }
  
}


f = 'DGE_results_PM_PH'
write.xlsx(as.data.frame(dgeResults), 'DGE_results_PM_PH.xlsx', rowNames = T)

dir.create(f, showWarnings=TRUE, recursive=TRUE)

lapply(
  names(dgeResults),
  function(x) write.table(
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname),
    file.path(f, paste(x, ".tsv", sep="")),
    append=F,
    row.names=F,
    col.names=T,
    quote=F,
    sep="\t"))

dgeResults_table = list()
dgeResults_table = lapply(
  names(dgeResults),
  function(x)
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname))

names(dgeResults_table) = names(dgeResults)

write.xlsx(dgeResults_table,
           file = paste(f,'/DGE_results_PM_PH.xlsx', sep=''),
           rowNames = F,
           asTable = T,
           startRow = 1,
           sheetName = str_sub(names(dgeResults),1,31))


n.label = 20
FDR = T
pvalue = 0.01
for (Condition in names(dgeResults)) {
  results = as.data.frame(dgeResults[[Condition]])
  results$DE = 'unm'
  if (!FDR) {
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'up'}
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < .1,]$DE = 'down' }
  } else {
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'up'}
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'down' }
  }
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'SEQCup'}
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]$DE = 'SEQCdown' }
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'FDRup'}
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'FDRdown' }
  
  results$DE <- factor(x = results$DE, levels = c("unm", "FDRdown","FDRup", 'SEQCdown','SEQCup'))
  mycolors = c('grey','dodgerblue4','darkred','dodgerblue2','coral'); names(mycolors) = levels(results$DE)
  results$DE2 = 'unm'; results[results$DE!='unm',]$DE2 = 'mod'
  results$DE2 <- factor(x = results$DE2, levels = c("mod","unm"))
  mysize = c(3,2); names(mysize) = levels(results$DE2)
  myalpha = c(1,0.2); names(mysize) = unique(results$DE2)
  
  # label N genes
  N = min(n.label, length(rownames(results[results$DE == 'FDRup',])))
  up_label = rownames(results[results$DE == 'FDRup',])[1:N]
  N = min(n.label, length(rownames(results[results$DE == 'FDRdown',])))
  down_label = rownames(results[results$DE == 'FDRdown',])[1:N]
  
  NSEQC = min(n.label, length(rownames(results[results$DE == 'SEQCup',])))
  SEQCup_label = rownames(results[results$DE == 'SEQCup',])[1:N]
  NSEQC = min(n.label, length(rownames(results[results$DE == 'SEQdown',])))
  SEQCdown_label = rownames(results[results$DE == 'SEQdown',])[1:N]
  
  
  MAplot = ggplot(results) +
    geom_point(aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    xlim(c(0,1.e5)) +
    scale_x_continuous(trans='log10') +
    ggtitle(paste("MAPlot,", Condition)) +
    scale_color_manual(values = mycolors) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "Mean expression", y = "log2 fold change")
  
  print(MAplot)
  pdf(paste('MAplot_',Condition,'.pdf',sep=''),width=8, height=6.5)
  plot(MAplot)
  dev.off()
  
  # Vulcano plot
  VP = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Vulcano Plot,", Condition)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results["Condition",], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results["Condition",]), size = 3,
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")

  pdf(paste(f,'/','VolcanoPlot',Condition,'.pdf',sep=''),width=8, height=6.5)
  plot(VP)
  dev.off()
  
  options(repr.plot.width=14, repr.plot.height=6.5)
  p1 = MAplot ; p2 = VP;
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +
          plot_layout(guides = "collect"))
  
  pdf(paste(f,'/','MA_VP_',Condition,'.pdf',sep=''),width=16, height=6.5)
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +
          plot_layout(guides = "collect"))
  dev.off()
  
}

DEGs <- as.data.frame(dgeResults_table$Condition)
rownames(DEGs) <- DEGs$Geneid
colnames(DEGs) <- c("GeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue","padj")

df <- DEGs
rownames(df) <- df$Geneid 

# Set thresholds
logFC_threshold <- 0.5
pval_threshold <- 0.05

# Create -log10(padj) for the plot
df$neg_log10_padj <- -log10(df$padj)

# Define color based on thresholds
df$color <- "grey"
df$color[df$log2FoldChange > logFC_threshold & df$padj < pval_threshold] <- "red"
df$color[df$log2FoldChange < -logFC_threshold & df$padj < pval_threshold] <- "blue"


pdf("PM_vs_PH_volcano_l2fc0.5_labels.pdf", )
ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = color), alpha = 0.75) +
  scale_color_identity() +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  geom_text_repel(aes(label = ifelse(abs(log2FoldChange) > logFC_threshold & padj < pval_threshold, GeneID, "")), 
                  box.padding = 0.3, max.overlaps = 10)+
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", title = "Volcano Plot")
dev.off()


logFC_threshold <- 1
pval_threshold <- 0.05

df$neg_log10_padj <- -log10(df$padj)

# Define color based on thresholds
df$color <- "grey"
df$color[df$log2FoldChange > logFC_threshold & df$padj < pval_threshold] <- "red"
df$color[df$log2FoldChange < -logFC_threshold & df$padj < pval_threshold] <- "blue"

table(df$color)

pdf("PM_vs_PH_volcano_l2fc1_labels.pdf", )
# Create the volcano plot
ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = color), alpha = 0.75) +
  scale_color_identity() +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  geom_text_repel(aes(label = ifelse(abs(log2FoldChange) > logFC_threshold & padj < pval_threshold, GeneID, "")), 
                  box.padding = 0.3, max.overlaps = 10)+
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", title = "Volcano Plot")
dev.off()


highlighted_genes <- c("Oasl2", "Isg15", "Irf7", "Zbp1", "Rtp4", "Ifit1", "Cxcl9", 
                       "Oas3", "Oas2", "Ifi205", "Cxcl10", "Ligp1", "Trim30c", 
                       "Gm4951", "Gm12250", "Ifi44", "Clic5", "Ifit3", "Ifi47", 
                       "Saa3", "Xaf1", "Tgtp2", "Ifi211", "Cmpk2", "Gbp6", "Trim30a", 
                       "Slfn5", "Irgm1", "Samd9l", "Igtp", "Gbp3", "Oasl1", "Lyi6i", 
                       "Rsad2", "Slc40a1", "Parp14", "Gbp2", "Clec7a", "Emp1", 
                       "Npy", "Cd28", "Ccl7", "Cdkn1c", "Hist1h1d", "Il6", "Il1b","Il12b",
                       "Sqle", "Cyp51","Msmo1", "Mmp8", "Lgals3", "Plaur")

pdf("PM_vs_PH_volcano_l2fc1_labels_new.pdf" )
ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = color), alpha = 0.75) +
  scale_color_identity() +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  geom_text_repel(aes(label = ifelse(GeneID %in% highlighted_genes, GeneID, "")), 
                  box.padding = 0.3, max.overlaps = Inf) +  # Ensures all highlighted genes are labeled
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", title = "Volcano Plot")
table(df$color)
dev.off()

df <- DEGs
# Set thresholds
logFC_threshold <- 0
pval_threshold <- 0.05
df$neg_log10_padj <- -log10(df$padj)

# Define color based on thresholds
df$color <- "grey"
df$color[df$log2FoldChange > logFC_threshold & df$padj < pval_threshold] <- "red"
df$color[df$log2FoldChange < -logFC_threshold & df$padj < pval_threshold] <- "blue"

table(df$color)

# Define the genes to label (from the provided volcano plot)

highlighted_genes <- c("Oasl2", "Isg15", "Irf7", "Zbp1", "Rtp4", "Ifit1", "Cxcl9", 
                       "Oas3", "Oas2", "Ifi205", "Cxcl10", "Ligp1", "Trim30c", 
                       "Gm4951", "Gm12250", "Ifi44", "Clic5", "Ifit3", "Ifi47", 
                       "Saa3", "Xaf1", "Tgtp2", "Ifi211", "Cmpk2", "Gbp6", "Trim30a", 
                       "Slfn5", "Irgm1", "Samd9l", "Igtp", "Gbp3", "Oasl1", "Lyi6i", 
                       "Rsad2", "Slc40a1", "Parp14", "Gbp2", "Clec7a", "Emp1", 
                       "Npy", "Cd28", "Ccl7", "Cdkn1c", "Hist1h1d", "Il6", "Il1b","Il12b",
                       "Sqle", "Cyp51","Msmo1")

pdf("PM_vs_PH_volcano_l2fc0_labels_neww.pdf" )
ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = color), alpha = 0.75) +
  scale_color_identity() +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  geom_text_repel(aes(label = ifelse(GeneID %in% highlighted_genes, GeneID, "")), 
                  box.padding = 0.3, max.overlaps = Inf) +  # Ensures all highlighted genes are labeled
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", title = "Volcano Plot")
table(df$color)
dev.off()



significant_genes <- df[df$padj < 0.05, ]
significant_up <- significant_genes[significant_genes$log2FoldChange > 0, ]
nrow(significant_up)


significant_genes <- df[df$padj < 0.05, ]
significant_dw <- significant_genes[significant_genes$log2FoldChange < 0, ]
nrow(significant_dw)


significant <- rbind(significant_up, significant_dw)
significant_genes <- significant$GeneID

ccc <- fCountsRPKM[rownames(fCountsRPKM) %in% significant_genes,]

annotation_column_ <- metadata

annotation_column_$Sample <- as.character(annotation_column_$Sample)
annotation_column_$Sample.Name <- NULL

annotation_column_$Sample <- as.character(annotation_column_$Sample)

ann_colors = list(
  Sample = c("1" = "forestgreen",
             "3" = "brown",
             "4" = "#999900"),
  Condition = c("PH"="orange", "PM" = "pink"))


HP <- pheatmap::pheatmap(ccc,
                         scale = 'row',
                         annotation_col = annotation_column_,
                         cluster_rows = T, 
                         cluster_cols = T, 
                         annotation_colors = ann_colors, 
                         show_rownames=T,
                         show_colnames=T,
                         main = "DEGs",
                         treeheight_row  = 0,
                         display_numbers = F, 
                         border_color=NA,
                         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                                               "RdYlBu")))(100),
                         filename = 'DEGs_heatmap.pdf',
                         width = 7, height = 120)

HP <- pheatmap::pheatmap(ccc,
                         scale = 'row',
                         annotation_col = annotation_column_,
                         cluster_rows = T, 
                         cluster_cols = T, 
                         annotation_colors = ann_colors, 
                         show_rownames=F,
                         show_colnames=F,
                         main = "DEGs",
                         treeheight_row  = 0,
                         display_numbers = F, 
                         border_color=NA,
                         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                                               "RdYlBu")))(100),
                         filename = 'DEGs_heatmap_small.pdf')

#### Perform GSEA analysis ####

original_gene_list <- DEGs$log2FoldChange
names(original_gene_list) <- DEGs$GeneID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

names(gene_list)
msigdbr_df <- msigdbr(species = "Mus musculus", category = "C2")

df <- msigdbr_df[,c("gs_name","gene_symbol")]
msigdbr_df <- df 

gsea_results_C2 <- GSEA(gene_list, 
                        TERM2GENE = msigdbr_df,  # Replace 'C1' with the actual gene set data frame
                        verbose = FALSE, 
                        pvalueCutoff = 0.05, 
                        maxGSSize = 5000, 
                        minGSSize = 5, 
                        pAdjustMethod = "BH", 
                        eps = 0)


write.xlsx(gsea_results_C2[gsea_results_C2$Description], "C2_gsea_CLL.xlsx")



msigdbr_df <- msigdbr(species = "Mus musculus", category = "C5")
df <- msigdbr_df[,c("gs_name","gene_symbol")]
msigdbr_df <- df 


gsea_results_C5 <- GSEA(gene_list, 
                        TERM2GENE = msigdbr_df,  # Replace 'C1' with the actual gene set data frame
                        verbose = FALSE, 
                        pvalueCutoff = 0.05, 
                        maxGSSize = 5000, 
                        minGSSize = 5, 
                        pAdjustMethod = "BH", 
                        eps = 0)


write.xlsx(gsea_results_C5[gsea_results_C5$Description], "C5_gsea.xlsx")



msigdbr_df <- msigdbr(species = "Mus musculus", category = "C7")
df <- msigdbr_df[,c("gs_name","gene_symbol")]
msigdbr_df <- df 


gsea_results_C7 <- GSEA(gene_list, 
                        TERM2GENE = msigdbr_df,  # Replace 'C1' with the actual gene set data frame
                        verbose = FALSE, 
                        pvalueCutoff = 0.05, 
                        maxGSSize = 5000, 
                        minGSSize = 5, 
                        pAdjustMethod = "BH", 
                        eps = 0)

write.xlsx(gsea_results_C7[gsea_results_C7$Description], "C7_gsea.xlsx")




pdf('C2_KEGG_STEROID_BIOSYNTHESIS.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C2, geneSetID = which(gsea_results_C2$Description == "KEGG_STEROID_BIOSYNTHESIS"), 
          title = "KEGG_STEROID_BIOSYNTHESIS")
dev.off()


pdf('C2_WP_TOLLLIKE_RECEPTOR_SIGNALING_PATHWAY.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C2, geneSetID = which(gsea_results_C2$Description == "WP_TOLLLIKE_RECEPTOR_SIGNALING_PATHWAY"), 
          title = "WP_TOLLLIKE_RECEPTOR_SIGNALING_PATHWAY")
dev.off()


pdf('C5_GOBP_INNATE_IMMUNE_RESPONSE.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C5, geneSetID = which(gsea_results_C5$Description == "GOBP_INNATE_IMMUNE_RESPONSE"), 
          title = "GOBP_INNATE_IMMUNE_RESPONSE")
dev.off()

pdf('C5_GOBP_STEROL_BIOSYNTHETIC_PROCESS.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C5, geneSetID = which(gsea_results_C5$Description == "GOBP_STEROL_BIOSYNTHETIC_PROCESS"), 
          title = "GOBP_STEROL_BIOSYNTHETIC_PROCESS")
dev.off()

pdf('C7_GSE14308_TH17_VS_NAIVE_CD4_TCELL_UP.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C7, geneSetID = which(gsea_results_C7$Description == "GSE14308_TH17_VS_NAIVE_CD4_TCELL_UP"), 
          title = "GSE14308_TH17_VS_NAIVE_CD4_TCELL_UP")
dev.off()

pdf('C2_REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C2, geneSetID = which(gsea_results_C2$Description == "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"), 
          title = "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM")
dev.off()

#### Perform Enrichment ####

significant_genes <- DEGs[DEGs$padj < 0.05, ]
significant_up <- significant_genes[significant_genes$log2FoldChange > 0, ]
up <- rownames(significant_up)
length(up)

significant_genes <- DEGs[DEGs$padj < 0.05, ]
significant_dw <- significant_genes[significant_genes$log2FoldChange < 0, ]
dw <- rownames(significant_dw)
length(dw)

upregulated_genes_entrez <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

go_enrichment_BP_up <- enrichGO(gene = upregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "BP",  
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)

go_enrichment_CC_up <- enrichGO(gene = upregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "CC",  
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)

go_enrichment_MF_up <- enrichGO(gene = upregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "MF",  
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)



pdf("Enrichment_BP_up_barplot.pdf")
barplot(go_enrichment_BP_up, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_BP_up_dotplot.pdf")
dotplot(go_enrichment_BP_up, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_CC_up_barplot.pdf")
barplot(go_enrichment_CC_up, showCategory = 20)
dev.off()

pdf("Enrichment_CC_up_dotplot.pdf")
dotplot(go_enrichment_CC_up, showCategory = 20)
dev.off()

pdf("Enrichment_MF_up_barplot.pdf",12,10)
barplot(go_enrichment_MF_up, showCategory = 20)  # Show top 10 enriched terms
dev.off()


pdf("Enrichment_MF_up_dotplot.pdf", 12,10)
dotplot(go_enrichment_MF_up, showCategory = 20)  # Show top 10 enriched terms
dev.off()




go_enrichment_BP_up_sub <- go_enrichment_BP_up

go_enrichment_BP_up_sub@result <- go_enrichment_BP_up_sub@result[go_enrichment_BP_up_sub@result$ID %in% c("GO:0006695",
                                                                                                          "GO:1902653",
                                                                                                          "GO:0050900",
                                                                                                          "GO:0032103",
                                                                                                          "GO:0008202",
                                                                                                          "GO:1901615",
                                                                                                          "GO:0070371",
                                                                                                          "GO:0050727"),]



pdf("Enrichment_BP_up_dotplot_interesting_new8.pdf", 8,6)
dotplot(go_enrichment_BP_up_sub, showCategory = 30)
dev.off()



split_gene_ids <- strsplit(go_enrichment_BP_up$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_BP_up <- as.data.frame(go_enrichment_BP_up)

go_enrichment_BP_up$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_BP_up, file = "go_enrichment_BP_up_results.xlsx")




split_gene_ids <- strsplit(go_enrichment_CC_up$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_CC_up <- as.data.frame(go_enrichment_CC_up)

go_enrichment_CC_up$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_CC_up, file = "go_enrichment_CC_up_results.xlsx")



split_gene_ids <- strsplit(go_enrichment_MF_up$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_MF_up <- as.data.frame(go_enrichment_MF_up)

go_enrichment_MF_up$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_MF_up, file = "go_enrichment_MF_up_results.xlsx")





dwregulated_genes_entrez <- bitr(dw, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

go_enrichment_BP_dw <- enrichGO(gene = dwregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "BP",  # Choose "BP", "CC", or "MF"
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)

go_enrichment_CC_dw <- enrichGO(gene = dwregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "CC",  # Choose "BP", "CC", or "MF"
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)

go_enrichment_MF_dw <- enrichGO(gene = dwregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "MF",  # Choose "BP", "CC", or "MF"
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)



pdf("Enrichment_BP_dw_dotplot.pdf")
dotplot(go_enrichment_BP_dw, showCategory = 20)
dev.off()

pdf("Enrichment_BP_dw_barplot.pdf")
barplot(go_enrichment_BP_dw, showCategory = 20)
dev.off()

pdf("Enrichment_MF_dw_barplot.pdf",10,10)
barplot(go_enrichment_MF_dw, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_MF_dw_dotplot.pdf", 10,10)
dotplot(go_enrichment_MF_dw, showCategory = 20)
dev.off()

pdf("Enrichment_CC_dw_barplot.pdf")
barplot(go_enrichment_CC_dw, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_CC_dw_dotplot.pdf")
dotplot(go_enrichment_CC_dw, showCategory = 20)
dev.off()


go_enrichment_BP_dw_sub <- go_enrichment_BP_dw


go_enrichment_BP_dw_sub@result <- go_enrichment_BP_dw_sub@result[go_enrichment_BP_dw_sub@result$ID %in% c("GO:0051607",
                                                                                                          "GO:0032103",
                                                                                                          "GO:0035456",
                                                                                                          "GO:0001819",
                                                                                                          "GO:0046651",
                                                                                                          "GO:0042098",
                                                                                                          "GO:0034341",
                                                                                                          "GO:0019882"),]



pdf("Enrichment_BP_dw_dotplot_interesting_new8.pdf", 8,6)
dotplot(go_enrichment_BP_dw_sub, showCategory = 30)
dev.off()





split_gene_ids <- strsplit(go_enrichment_BP_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_BP_dw <- as.data.frame(go_enrichment_BP_dw)

go_enrichment_BP_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_BP_dw, file = "go_enrichment_BP_dw_results.xlsx")



split_gene_ids <- strsplit(go_enrichment_MF_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_MF_dw <- as.data.frame(go_enrichment_MF_dw)

go_enrichment_MF_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_MF_dw, file = "go_enrichment_MF_dw_results.xlsx")




split_gene_ids <- strsplit(go_enrichment_CC_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_CC_dw <- as.data.frame(go_enrichment_CC_dw)

go_enrichment_CC_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_CC_dw, file = "go_enrichment_CC_dw_results.xlsx")


#### Jaccard Plots ####


# Output directory
dir = paste("./JaccardPlots/", sep = '')
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

# Set p-value threshold
p.value.thr = 0.01

# Define breaks and colors for heatmap
breaksList = seq(0, 1, by = 0.001)
myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList))

# List all .xlsx files in the enrichment results directory
lf = "go_enrichment_BP_up_results.xlsx"

getwd
# Loop through each file
for (file in lf) {
  enrich.file = paste("./GO/", file, sep = '')
  
  # Read the Excel file
  sheets = openxlsx::getSheetNames(enrich.file)
  Pathways.Table = data.frame()
  
  # Read each sheet and combine into a single table
  for (dat in sheets) {
    Table <- read.xlsx(xlsxFile = enrich.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = FALSE,  # clusterProfiler results don’t have row names
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  # Filter pathways based on adjusted p-value
  Pathways.Table = Pathways.Table[Pathways.Table$p.adjust < p.value.thr,]
  pathways = Pathways.Table$Description
  
  # Extract gene lists for each pathway
  gene.list = list()
  gene.all = character()
  
  for (i in 1:nrow(Pathways.Table)) {
    genes = unlist(strsplit(Pathways.Table$geneID[i], '/'))  # clusterProfiler separates genes with "/"
    gene.list[[Pathways.Table$Description[i]]] = genes
    gene.all = c(gene.all, genes)
  }
  gene.all = unique(gene.all)
  
  # Create binary matrix for Jaccard calculation
  if (length(gene.all) != 0 & length(pathways) > 1) {
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all
    
    for (pat in pathways) {
      for (gene in gene.list[[pat]]) {
        if (gene %in% gene.all) {
          M[pat, gene] <- 1
        }
      }
    }
    
    # Compute Jaccard distance
    Jacard.Matrix <- distance(M, method = "jaccard")
    
    if (length(pathways) == 2) {
      Jacard.Matrix_new = as.matrix(rbind(c(0, Jacard.Matrix), c(Jacard.Matrix, 0)))
      Jacard.Matrix = Jacard.Matrix_new
    }
    
    row.names(Jacard.Matrix) <- pathways
    colnames(Jacard.Matrix) <- pathways
    
    # Heatmap parameters
    w = 30
    h = 80
    fs = 5
    cutree_rows_N = 9
    
    # Generate heatmap
    pheatmap(Jacard.Matrix,
             border_color = 'darkgrey',
             color = myc,
             breaks = breaksList,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = cutree_rows_N,
             show_colnames = FALSE,
             main = paste(file, '- Jaccard distance heatmap'),
             fontsize = 18,
             fontsize_row = fs,
             filename = paste(dir, file, '_JaccardDist.pdf', sep = ''),
             width = w, height = h)
  }
}




# Set p-value threshold
p.value.thr = 0.01

# Define breaks and colors for heatmap
breaksList = seq(0, 1, by = 0.001)
myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList))

# List all .xlsx files in the enrichment results directory
lf = "go_enrichment_BP_up_results.xlsx"

getwd
# Loop through each file
for (file in lf) {
  enrich.file = paste("./GO/", file, sep = '')
  
  # Read the Excel file
  sheets = openxlsx::getSheetNames(enrich.file)
  Pathways.Table = data.frame()
  
  # Read each sheet and combine into a single table
  for (dat in sheets) {
    Table <- read.xlsx(xlsxFile = enrich.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = FALSE,  # clusterProfiler results don’t have row names
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  # Filter pathways based on adjusted p-value
  Pathways.Table = Pathways.Table[Pathways.Table$p.adjust < p.value.thr,]
  pathways = Pathways.Table$Description
  
  # Extract gene lists for each pathway
  gene.list = list()
  gene.all = character()
  
  for (i in 1:nrow(Pathways.Table)) {
    genes = unlist(strsplit(Pathways.Table$geneID[i], '/'))  # clusterProfiler separates genes with "/"
    gene.list[[Pathways.Table$Description[i]]] = genes
    gene.all = c(gene.all, genes)
  }
  gene.all = unique(gene.all)
  
  # Create binary matrix for Jaccard calculation
  if (length(gene.all) != 0 & length(pathways) > 1) {
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all
    
    for (pat in pathways) {
      for (gene in gene.list[[pat]]) {
        if (gene %in% gene.all) {
          M[pat, gene] <- 1
        }
      }
    }
    
    # Compute Jaccard distance
    Jacard.Matrix <- distance(M, method = "jaccard")
    
    if (length(pathways) == 2) {
      Jacard.Matrix_new = as.matrix(rbind(c(0, Jacard.Matrix), c(Jacard.Matrix, 0)))
      Jacard.Matrix = Jacard.Matrix_new
    }
    
    row.names(Jacard.Matrix) <- pathways
    colnames(Jacard.Matrix) <- pathways
    
    # Heatmap parameters
    w = 10
    h = 10
    fs = 5
    cutree_rows_N = 9
    
    # Generate heatmap
    pheatmap(Jacard.Matrix,
             border_color = 'darkgrey',
             color = myc,
             breaks = breaksList,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = cutree_rows_N,
             show_colnames = FALSE,
             main = paste(file, '- Jaccard distance heatmap'),
             #fontsize = 12,
             #fontsize_row = fs,
             show_rownames = FALSE,
             filename = paste(dir, file, '_JaccardDist_small.pdf', sep = ''),
             width = w, height = h)
  }
}




# Output directory
dir = paste("./JaccardPlots/", sep = '')
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

# Set p-value threshold
p.value.thr = 0.01

# Define breaks and colors for heatmap
breaksList = seq(0, 1, by = 0.001)
myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList))

# List all .xlsx files in the enrichment results directory
lf = "go_enrichment_BP_dw_results.xlsx"

getwd
# Loop through each file
for (file in lf) {
  enrich.file = paste("./GO/", file, sep = '')
  
  # Read the Excel file
  sheets = openxlsx::getSheetNames(enrich.file)
  Pathways.Table = data.frame()
  
  # Read each sheet and combine into a single table
  for (dat in sheets) {
    Table <- read.xlsx(xlsxFile = enrich.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = FALSE,  # clusterProfiler results don’t have row names
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  # Filter pathways based on adjusted p-value
  Pathways.Table = Pathways.Table[Pathways.Table$p.adjust < p.value.thr,]
  pathways = Pathways.Table$Description
  
  # Extract gene lists for each pathway
  gene.list = list()
  gene.all = character()
  
  for (i in 1:nrow(Pathways.Table)) {
    genes = unlist(strsplit(Pathways.Table$geneID[i], '/'))  # clusterProfiler separates genes with "/"
    gene.list[[Pathways.Table$Description[i]]] = genes
    gene.all = c(gene.all, genes)
  }
  gene.all = unique(gene.all)
  
  # Create binary matrix for Jaccard calculation
  if (length(gene.all) != 0 & length(pathways) > 1) {
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all
    
    for (pat in pathways) {
      for (gene in gene.list[[pat]]) {
        if (gene %in% gene.all) {
          M[pat, gene] <- 1
        }
      }
    }
    
    # Compute Jaccard distance
    Jacard.Matrix <- distance(M, method = "jaccard")
    
    if (length(pathways) == 2) {
      Jacard.Matrix_new = as.matrix(rbind(c(0, Jacard.Matrix), c(Jacard.Matrix, 0)))
      Jacard.Matrix = Jacard.Matrix_new
    }
    
    row.names(Jacard.Matrix) <- pathways
    colnames(Jacard.Matrix) <- pathways
    
    # Heatmap parameters
    w = 30
    h = 80
    fs = 5
    cutree_rows_N = 9
    
    # Generate heatmap
    pheatmap(Jacard.Matrix,
             border_color = 'darkgrey',
             color = myc,
             breaks = breaksList,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = cutree_rows_N,
             show_colnames = FALSE,
             main = paste(file, '- Jaccard distance heatmap'),
             fontsize = 18,
             fontsize_row = fs,
             filename = paste(dir, file, '_JaccardDist.pdf', sep = ''),
             width = w, height = h)
  }
}






# Set p-value threshold
p.value.thr = 0.01

# Define breaks and colors for heatmap
breaksList = seq(0, 1, by = 0.001)
myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList))

# List all .xlsx files in the enrichment results directory
lf = "go_enrichment_BP_dw_results.xlsx"

getwd
# Loop through each file
for (file in lf) {
  enrich.file = paste("./GO/", file, sep = '')
  
  # Read the Excel file
  sheets = openxlsx::getSheetNames(enrich.file)
  Pathways.Table = data.frame()
  
  # Read each sheet and combine into a single table
  for (dat in sheets) {
    Table <- read.xlsx(xlsxFile = enrich.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = FALSE,  # clusterProfiler results don’t have row names
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  # Filter pathways based on adjusted p-value
  Pathways.Table = Pathways.Table[Pathways.Table$p.adjust < p.value.thr,]
  pathways = Pathways.Table$Description
  
  # Extract gene lists for each pathway
  gene.list = list()
  gene.all = character()
  
  for (i in 1:nrow(Pathways.Table)) {
    genes = unlist(strsplit(Pathways.Table$geneID[i], '/'))  # clusterProfiler separates genes with "/"
    gene.list[[Pathways.Table$Description[i]]] = genes
    gene.all = c(gene.all, genes)
  }
  gene.all = unique(gene.all)
  
  # Create binary matrix for Jaccard calculation
  if (length(gene.all) != 0 & length(pathways) > 1) {
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all
    
    for (pat in pathways) {
      for (gene in gene.list[[pat]]) {
        if (gene %in% gene.all) {
          M[pat, gene] <- 1
        }
      }
    }
    
    # Compute Jaccard distance
    Jacard.Matrix <- distance(M, method = "jaccard")
    
    if (length(pathways) == 2) {
      Jacard.Matrix_new = as.matrix(rbind(c(0, Jacard.Matrix), c(Jacard.Matrix, 0)))
      Jacard.Matrix = Jacard.Matrix_new
    }
    
    row.names(Jacard.Matrix) <- pathways
    colnames(Jacard.Matrix) <- pathways
    
    # Heatmap parameters
    w = 10
    h = 10
    fs = 5
    cutree_rows_N = 9
    
    # Generate heatmap
    pheatmap(Jacard.Matrix,
             border_color = 'darkgrey',
             color = myc,
             breaks = breaksList,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = cutree_rows_N,
             show_colnames = FALSE,
             main = paste(file, '- Jaccard distance heatmap'),
             #fontsize = 12,
             #fontsize_row = fs,
             show_rownames = FALSE,
             filename = paste(dir, file, '_JaccardDist_small.pdf', sep = ''),
             width = w, height = h)
  }
}


#### Try performing the enrichment changing the universe ####

expressed_genes <- rownames(DEGs)

expressed_genes_entrez <- bitr(expressed_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
upregulated_genes_entrez <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

head(expressed_genes_entrez)

go_enrichment_BP_up <- enrichGO(
  gene = upregulated_genes_entrez$ENTREZID, 
  OrgDb = org.Mm.eg.db, 
  keyType = "ENTREZID", 
  ont = "BP",  
  pAdjustMethod = "BH", 
  qvalueCutoff = 0.05,
  universe = expressed_genes_entrez$ENTREZID  # Qui metti la tua lista di background genes
)


pdf("Enrichment_BP_up_barplot_universe_exp.pdf", 15,10)
barplot(go_enrichment_BP_up, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_BP_up_dotplot_universe_exp.pdf",15,10)
dotplot(go_enrichment_BP_up, showCategory = 20)  # Show top 10 enriched terms
dev.off()



go_enrichment_CC_up <- enrichGO(gene = upregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "CC",  # Choose "BP", "CC", or "MF"
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                universe = expressed_genes_entrez$ENTREZID )


pdf("Enrichment_CC_up_barplot_universe_exp.pdf", 15,10)
barplot(go_enrichment_CC_up, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_CC_up_dotplot_universe_exp.pdf", 15,10)
dotplot(go_enrichment_CC_up, showCategory = 20)
dev.off()

go_enrichment_MF_up <- enrichGO(gene = upregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "MF",  # Choose "BP", "CC", or "MF"
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05,
                                universe = expressed_genes_entrez$ENTREZID )

pdf("Enrichment_MF_up_barplot_universe_exp.pdf", 15,10)
barplot(go_enrichment_MF_up, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_MF_up_dotplot_universe_exp.pdf", 15,10)
dotplot(go_enrichment_MF_up, showCategory = 20)
dev.off()



split_gene_ids <- strsplit(go_enrichment_BP_up$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_BP_up <- as.data.frame(go_enrichment_BP_up)

go_enrichment_BP_up$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_BP_up, file = "go_enrichment_BP_up_results_universe_exp.xlsx")




split_gene_ids <- strsplit(go_enrichment_MF_up$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_MF_up <- as.data.frame(go_enrichment_MF_up)

go_enrichment_MF_up$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_MF_up, file = "go_enrichment_MF_up_results_universe_exp.xlsx")




split_gene_ids <- strsplit(go_enrichment_CC_up$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_CC_up <- as.data.frame(go_enrichment_CC_up)

go_enrichment_CC_up$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_CC_up, file = "go_enrichment_CC_up_results_universe_exp.xlsx")





expressed_genes <- rownames(DEGs)
expressed_genes_entrez <- bitr(expressed_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
dwregulated_genes_entrez <- bitr(dw, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)


go_enrichment_BP_dw <- enrichGO(
  gene = dwregulated_genes_entrez$ENTREZID, 
  OrgDb = org.Mm.eg.db, 
  keyType = "ENTREZID", 
  ont = "BP",  
  pAdjustMethod = "BH", 
  qvalueCutoff = 0.05,
  universe = expressed_genes_entrez$ENTREZID  # Qui metti la tua lista di background genes
)



pdf("Enrichment_BP_dw_barplot_universe_exp.pdf", 15,10)
barplot(go_enrichment_BP_dw, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_BP_dw_dotplot_universe_exp.pdf",15,10)
dotplot(go_enrichment_BP_dw, showCategory = 20)  # Show top 10 enriched terms
dev.off()



go_enrichment_CC_dw <- enrichGO(gene = dwregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "CC",  # Choose "BP", "CC", or "MF"
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                universe = expressed_genes_entrez$ENTREZID )


pdf("Enrichment_CC_dw_barplot_universe_exp.pdf", 15,10)
barplot(go_enrichment_CC_dw, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_CC_dw_dotplot_universe_exp.pdf", 15,10)
dotplot(go_enrichment_CC_dw, showCategory = 20)
dev.off()

go_enrichment_MF_dw <- enrichGO(gene = dwregulated_genes_entrez$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                keyType = "ENTREZID", 
                                ont = "MF",  # Choose "BP", "CC", or "MF"
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05,
                                universe = expressed_genes_entrez$ENTREZID )

pdf("Enrichment_MF_dw_barplot_universe_exp.pdf", 15,10)
barplot(go_enrichment_MF_dw, showCategory = 20)  # Show top 10 enriched terms
dev.off()

pdf("Enrichment_MF_dw_dotplot_universe_exp.pdf", 15,10)
dotplot(go_enrichment_MF_dw, showCategory = 20)
dev.off()



split_gene_ids <- strsplit(go_enrichment_BP_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_BP_dw <- as.data.frame(go_enrichment_BP_dw)

go_enrichment_BP_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_BP_dw, file = "go_enrichment_BP_dw_results_universe_exp.xlsx")




split_gene_ids <- strsplit(go_enrichment_MF_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_MF_dw <- as.data.frame(go_enrichment_MF_dw)

go_enrichment_MF_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_MF_dw, file = "go_enrichment_MF_dw_results_universe_exp.xlsx")



split_gene_ids <- strsplit(go_enrichment_CC_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Mm.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_CC_dw <- as.data.frame(go_enrichment_CC_dw)

go_enrichment_CC_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_CC_dw, file = "go_enrichment_CC_dw_results_universe_exp.xlsx")

# Output directory
dir = paste("./JaccardPlots_universe_expressed_genes/", sep = '')
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

# Set p-value threshold
p.value.thr = 0.01

# Define breaks and colors for heatmap
breaksList = seq(0, 1, by = 0.001)
myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList))

# List all .xlsx files in the enrichment results directory
lf = "go_enrichment_BP_up_results_universe_exp.xlsx"

getwd
# Loop through each file
for (file in lf) {
  enrich.file = paste("./GO_universe_expressed_genes/", file, sep = '')
  
  # Read the Excel file
  sheets = openxlsx::getSheetNames(enrich.file)
  Pathways.Table = data.frame()
  
  # Read each sheet and combine into a single table
  for (dat in sheets) {
    Table <- read.xlsx(xlsxFile = enrich.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = FALSE,  # clusterProfiler results don’t have row names
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  # Filter pathways based on adjusted p-value
  Pathways.Table = Pathways.Table[Pathways.Table$p.adjust < p.value.thr,]
  pathways = Pathways.Table$Description
  
  # Extract gene lists for each pathway
  gene.list = list()
  gene.all = character()
  
  for (i in 1:nrow(Pathways.Table)) {
    genes = unlist(strsplit(Pathways.Table$geneID[i], '/'))  # clusterProfiler separates genes with "/"
    gene.list[[Pathways.Table$Description[i]]] = genes
    gene.all = c(gene.all, genes)
  }
  gene.all = unique(gene.all)
  
  # Create binary matrix for Jaccard calculation
  if (length(gene.all) != 0 & length(pathways) > 1) {
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all
    
    for (pat in pathways) {
      for (gene in gene.list[[pat]]) {
        if (gene %in% gene.all) {
          M[pat, gene] <- 1
        }
      }
    }
    
    # Compute Jaccard distance
    Jacard.Matrix <- distance(M, method = "jaccard")
    
    if (length(pathways) == 2) {
      Jacard.Matrix_new = as.matrix(rbind(c(0, Jacard.Matrix), c(Jacard.Matrix, 0)))
      Jacard.Matrix = Jacard.Matrix_new
    }
    
    row.names(Jacard.Matrix) <- pathways
    colnames(Jacard.Matrix) <- pathways
    
    # Heatmap parameters
    w = 30
    h = 80
    fs = 5
    cutree_rows_N = 9
    
    # Generate heatmap
    pheatmap(Jacard.Matrix,
             border_color = 'darkgrey',
             color = myc,
             breaks = breaksList,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = cutree_rows_N,
             show_colnames = FALSE,
             main = paste(file, '- Jaccard distance heatmap'),
             fontsize = 18,
             fontsize_row = fs,
             filename = paste(dir, file, '_JaccardDist.pdf', sep = ''),
             width = w, height = h)
  }
}




# Set p-value threshold
p.value.thr = 0.01

# Define breaks and colors for heatmap
breaksList = seq(0, 1, by = 0.001)
myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList))

# List all .xlsx files in the enrichment results directory
lf = "go_enrichment_BP_up_results_universe_exp.xlsx"

getwd
# Loop through each file
for (file in lf) {
  enrich.file = paste("./GO_universe_expressed_genes/", file, sep = '')
  
  # Read the Excel file
  sheets = openxlsx::getSheetNames(enrich.file)
  Pathways.Table = data.frame()
  
  # Read each sheet and combine into a single table
  for (dat in sheets) {
    Table <- read.xlsx(xlsxFile = enrich.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = FALSE,  
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  # Filter pathways based on adjusted p-value
  Pathways.Table = Pathways.Table[Pathways.Table$p.adjust < p.value.thr,]
  pathways = Pathways.Table$Description
  
  # Extract gene lists for each pathway
  gene.list = list()
  gene.all = character()
  
  for (i in 1:nrow(Pathways.Table)) {
    genes = unlist(strsplit(Pathways.Table$geneID[i], '/'))  # clusterProfiler separates genes with "/"
    gene.list[[Pathways.Table$Description[i]]] = genes
    gene.all = c(gene.all, genes)
  }
  gene.all = unique(gene.all)
  
  # Create binary matrix for Jaccard calculation
  if (length(gene.all) != 0 & length(pathways) > 1) {
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all
    
    for (pat in pathways) {
      for (gene in gene.list[[pat]]) {
        if (gene %in% gene.all) {
          M[pat, gene] <- 1
        }
      }
    }
    
    # Compute Jaccard distance
    Jacard.Matrix <- distance(M, method = "jaccard")
    
    if (length(pathways) == 2) {
      Jacard.Matrix_new = as.matrix(rbind(c(0, Jacard.Matrix), c(Jacard.Matrix, 0)))
      Jacard.Matrix = Jacard.Matrix_new
    }
    
    row.names(Jacard.Matrix) <- pathways
    colnames(Jacard.Matrix) <- pathways
    
    # Heatmap parameters
    w = 10
    h = 10
    fs = 5
    cutree_rows_N = 9
    
    # Generate heatmap
    pheatmap(Jacard.Matrix,
             border_color = 'darkgrey',
             color = myc,
             breaks = breaksList,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = cutree_rows_N,
             show_colnames = FALSE,
             main = paste(file, '- Jaccard distance heatmap'),
             #fontsize = 12,
             #fontsize_row = fs,
             show_rownames = FALSE,
             filename = paste(dir, file, '_JaccardDist_small.pdf', sep = ''),
             width = w, height = h)
  }
}



# Output directory
dir = paste("./JaccardPlots/", sep = '')
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

# Set p-value threshold
p.value.thr = 0.01

# Define breaks and colors for heatmap
breaksList = seq(0, 1, by = 0.001)
myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList))

# List all .xlsx files in the enrichment results directory
lf = "go_enrichment_BP_dw_results_universe_exp.xlsx"

getwd
# Loop through each file
for (file in lf) {
  enrich.file = paste("./GO_universe_expressed_genes/", file, sep = '')
  
  # Read the Excel file
  sheets = openxlsx::getSheetNames(enrich.file)
  Pathways.Table = data.frame()
  
  # Read each sheet and combine into a single table
  for (dat in sheets) {
    Table <- read.xlsx(xlsxFile = enrich.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = FALSE,  # clusterProfiler results don’t have row names
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  # Filter pathways based on adjusted p-value
  Pathways.Table = Pathways.Table[Pathways.Table$p.adjust < p.value.thr,]
  pathways = Pathways.Table$Description
  
  # Extract gene lists for each pathway
  gene.list = list()
  gene.all = character()
  
  for (i in 1:nrow(Pathways.Table)) {
    genes = unlist(strsplit(Pathways.Table$geneID[i], '/'))  # clusterProfiler separates genes with "/"
    gene.list[[Pathways.Table$Description[i]]] = genes
    gene.all = c(gene.all, genes)
  }
  gene.all = unique(gene.all)
  
  # Create binary matrix for Jaccard calculation
  if (length(gene.all) != 0 & length(pathways) > 1) {
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all
    
    for (pat in pathways) {
      for (gene in gene.list[[pat]]) {
        if (gene %in% gene.all) {
          M[pat, gene] <- 1
        }
      }
    }
    
    # Compute Jaccard distance
    Jacard.Matrix <- distance(M, method = "jaccard")
    
    if (length(pathways) == 2) {
      Jacard.Matrix_new = as.matrix(rbind(c(0, Jacard.Matrix), c(Jacard.Matrix, 0)))
      Jacard.Matrix = Jacard.Matrix_new
    }
    
    row.names(Jacard.Matrix) <- pathways
    colnames(Jacard.Matrix) <- pathways
    
    # Heatmap parameters
    w = 30
    h = 80
    fs = 5
    cutree_rows_N = 9
    
    # Generate heatmap
    pheatmap(Jacard.Matrix,
             border_color = 'darkgrey',
             color = myc,
             breaks = breaksList,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = cutree_rows_N,
             show_colnames = FALSE,
             main = paste(file, '- Jaccard distance heatmap'),
             fontsize = 18,
             fontsize_row = fs,
             filename = paste(dir, file, '_JaccardDist.pdf', sep = ''),
             width = w, height = h)
  }
}


# Set p-value threshold
p.value.thr = 0.01

# Define breaks and colors for heatmap
breaksList = seq(0, 1, by = 0.001)
myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(breaksList))

# List all .xlsx files in the enrichment results directory
lf = "go_enrichment_BP_dw_results_universe_exp.xlsx"

getwd
# Loop through each file
for (file in lf) {
  enrich.file = paste("./GO_universe_expressed_genes/", file, sep = '')
  
  # Read the Excel file
  sheets = openxlsx::getSheetNames(enrich.file)
  Pathways.Table = data.frame()
  
  # Read each sheet and combine into a single table
  for (dat in sheets) {
    Table <- read.xlsx(xlsxFile = enrich.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = FALSE,  # clusterProfiler results don’t have row names
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  # Filter pathways based on adjusted p-value
  Pathways.Table = Pathways.Table[Pathways.Table$p.adjust < p.value.thr,]
  pathways = Pathways.Table$Description
  
  # Extract gene lists for each pathway
  gene.list = list()
  gene.all = character()
  
  for (i in 1:nrow(Pathways.Table)) {
    genes = unlist(strsplit(Pathways.Table$geneID[i], '/'))  # clusterProfiler separates genes with "/"
    gene.list[[Pathways.Table$Description[i]]] = genes
    gene.all = c(gene.all, genes)
  }
  gene.all = unique(gene.all)
  
  # Create binary matrix for Jaccard calculation
  if (length(gene.all) != 0 & length(pathways) > 1) {
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all
    
    for (pat in pathways) {
      for (gene in gene.list[[pat]]) {
        if (gene %in% gene.all) {
          M[pat, gene] <- 1
        }
      }
    }
    
    # Compute Jaccard distance
    Jacard.Matrix <- distance(M, method = "jaccard")
    
    if (length(pathways) == 2) {
      Jacard.Matrix_new = as.matrix(rbind(c(0, Jacard.Matrix), c(Jacard.Matrix, 0)))
      Jacard.Matrix = Jacard.Matrix_new
    }
    
    row.names(Jacard.Matrix) <- pathways
    colnames(Jacard.Matrix) <- pathways
    
    # Heatmap parameters
    w = 10
    h = 10
    fs = 5
    cutree_rows_N = 9
    
    # Generate heatmap
    pheatmap(Jacard.Matrix,
             border_color = 'darkgrey',
             color = myc,
             breaks = breaksList,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = cutree_rows_N,
             show_colnames = FALSE,
             main = paste(file, '- Jaccard distance heatmap'),
             #fontsize = 12,
             #fontsize_row = fs,
             show_rownames = FALSE,
             filename = paste(dir, file, '_JaccardDist_small.pdf', sep = ''),
             width = w, height = h)
  }
}

