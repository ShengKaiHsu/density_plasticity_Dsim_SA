# Title     : General script
# Objective : Implementation of whole analysis
# Created by: bunibal
# Created on: 14.09.21
#modified by SKH 17.09.21

rm(list=ls())
# DATASOURCE <- "/home/bunibal/PycharmProjects/bachelor/datadensity/"
# GRAPHDIRECTORY <- "/home/bunibal/PycharmProjects/bachelor/graphs/"
# TABLEDIRECTORY <- "/home/bunibal/PycharmProjects/bachelor/tables/"

DATASOURCE <- "/Users/shengkai/Dropbox (PopGen)/popgen_vienna/Stephan Buchner/Pt_density_plasticity/"
GRAPHDIRECTORY <- "/Users/shengkai/Dropbox (PopGen)/popgen_vienna/Stephan Buchner/Pt_density_plasticity/graph/"
TABLEDIRECTORY <- "/Users/shengkai/Dropbox (PopGen)/popgen_vienna/Stephan Buchner/Pt_density_plasticity/table/"

library(edgeR)
library(topGO)
library(WGCNA)
library(EnhancedVolcano)
library(VennDiagram)

####cumstomized function####
makeContrast <- function(matrix, my_contrast, contrast) {
  lrt_res <- glmLRT(matrix, contrast = my_contrast[, contrast])
  lrt_res$table$padj <- p.adjust(lrt_res$table$PValue, method = "BH") # Benjamini Hochbergs correction
  res_table <- lrt_res$table
}

filtering <- function(DGE, res_table_name_list, nameslist, pvalue = 0.05, minchange = log2(1.25)) {
  gene_list <- vector("list", length(res_table_name_list) * 2)
  nam <- vector("list", length(res_table_name_list) * 2)
  i <- 1
  for (y in nameslist) {
    nam[[i]] <- paste0(y, "_up")
    i <- i + 1
    nam[[i]] <- paste0(y, "_dn")
    i <- i + 1
  }
  i <- 1
  for (res in res_table_name_list) {
    gene_list[[i]] <- rownames(DGE)[res$padj < pvalue & res$logFC > minchange]
    i <- i + 1
    gene_list[[i]] <- rownames(DGE)[res$padj < pvalue & res$logFC < -minchange]
    i <- i + 1
  }
  names(gene_list) <- nam
  gene_list
}

create_query_ID_cat <- function(query, index_1, background_1, index_2, background_2, DGE, res_tab1, res_tab2, pvalue = 0.05, minchange = log2(1.25)) {
  thislist <- list(a_up = setdiff(query[[index_1]], unlist(background_1[-index_1])),
                   a_down = setdiff(query[[index_1 + 1]], unlist(background_1[-(index_1 + 1)])),
                   b_up = setdiff(query[[index_2]], unlist(background_2[-index_2])),
                   b_dn = setdiff(query[[index_2 + 1]], unlist(background_2[-(index_2 + 1)])),
                   exa_up = rownames(DGE)[res_tab2$padj < pvalue &
                                            res_tab1$padj < pvalue &
                                            res_tab1$logFC > minchange &
                                            res_tab2$logFC > minchange],
                   exa_dn = rownames(DGE)[res_tab2$padj < pvalue &
                                            res_tab1$padj < pvalue &
                                            res_tab1$logFC < -minchange &
                                            res_tab2$logFC < -minchange],
                   di_up = rownames(DGE)[res_tab2$padj < pvalue &
                                           res_tab1$padj < pvalue &
                                           res_tab1$logFC > minchange &
                                           res_tab2$logFC < -minchange],
                   di_dn = rownames(DGE)[res_tab2$padj < pvalue &
                                           res_tab1$padj < pvalue &
                                           res_tab1$logFC < -minchange &
                                           res_tab2$logFC > minchange])
  return(thislist)
}

gogoanalysis <- function(input) {
  GO_res_table <- list()
  for (i in names(input)) {
    tmp <- as.integer(rownames(y) %in% unlist(input[[i]]))
    names(tmp) <- rownames(y)
    if (sum(tmp) == 0) {
      GO_res_table[[i]] <- FALSE
      next
    }
    tmp <- factor(tmp)
    tgd <- new("topGOdata", ontology = "BP", allGenes = tmp, nodeSize = 5, annot = annFUN.org, mapping = "org.Dm.eg.db", ID = "ensembl")
    resTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")
    resTopGO.weight01 <- runTest(tgd, algorithm = "weight01", statistic = "Fisher")
    tmp_res <- GenTable(tgd, Fisher.classic = resTopGO.classic, Fisher.weight01 = resTopGO.weight01, orderBy = "Fisher.weight01", ranksOf = "Fisher.classic", topNodes = length(resTopGO.classic@score), numChar = 100)
    GO_res_table[[i]] <- tmp_res
  }
  GO_res_table
}

filtering_GO <- function(input) { lapply(input, function(x) {
  if (x == FALSE) {
    return(x)
  }
  x$Fisher.weight01[x$Fisher.weight01 == "< 1e-30"] <- 1e-30
  return(x) })
}

####data_input####
# Input the acclimation-data

dat_use <- read.csv(paste0(DATASOURCE, "countedreads_togenes.csv"), stringsAsFactors = F, row.names = 1)
counts_use <- dat_use[apply(cpm(dat_use), 1, function(x) !sum(x < 1) >= 1),] # filtered for lowly expressed genes
density <- strsplit2(colnames(dat_use), "_")[, 2] # translate labels to Density : HighD, lowD
harv_time <- strsplit2(colnames(dat_use), "_")[, 3] # harvest time : lowD 0 , HighD 24, 48, 72 h
harv_time[!harv_time %in% c("I", "II", "III")] <- "0"
y <- DGEList(counts = counts_use, group = density)
y <- calcNormFactors(y)

# implement a PCA
pca <- prcomp(t(log(cpm(y))))
ve <- round(pca$sdev^2 / sum(pca$sdev^2), 4)

####DE analysis####
# DE - analysis
group <- paste0(density, harv_time)
ModelDesign <- model.matrix(~0 + group)
DGE <- estimateDisp(y, design = ModelDesign, robust = T)
GLM <- glmFit(DGE, design = ModelDesign)
mycontrast <- makeContrasts("density" = ((groupHighDI + groupHighDII + groupHighDIII) / 3) - grouplowD0, "low~I" = grouplowD0 - groupHighDI, "low~II" = grouplowD0 - groupHighDII, "low~III" = grouplowD0 - groupHighDIII,
                            "I~II" = groupHighDI - groupHighDII, "II~III" = groupHighDII - groupHighDIII, "I~III" = groupHighDI - groupHighDIII, "lowD" = grouplowD0, "highD" = (groupHighDI + groupHighDII + groupHighDIII) / 3, levels = ModelDesign)

res_list <- list(
  res_table_density = makeContrast(GLM, mycontrast, "density"),
  res_table_lowI = makeContrast(GLM, mycontrast, "low~I"),
  res_table_lowII = makeContrast(GLM, mycontrast, "low~II"),
  res_table_lowIII = makeContrast(GLM, mycontrast, "low~III"),
  res_table_I_II = makeContrast(GLM, mycontrast, "I~II"),
  res_table_I_III = makeContrast(GLM, mycontrast, "I~III"),
  res_table_II_III = makeContrast(GLM, mycontrast, "II~III"),
  res_table_lowD = makeContrast(GLM, mycontrast, "lowD"),
  res_table_highD = makeContrast(GLM, mycontrast, "highD")
)

write.csv(res_list$res_table_density,paste0(TABLEDIRECTORY,"/Table_S1.csv"),quote = F)
# seperate up and dn regulated genes
query_ID_den <- filtering(y, res_list, c("density", "low~I", "low~II", "low~III", "I~II", "I~III", "II~III"))


#### WGCNA ####
# inverting the dataframe...
# nam <- colnames(dat_use)
# counts <- data.frame(t(dat_use), row.names = nam)
# 
# y <- DGEList(counts = t(counts), group = paste(density, harv_time))
# y <- calcNormFactors(y)

expr_dat <- t(log10(cpm(y)))

net_dat <- WGCNA::blockwiseModules(expr_dat, power = 6, deepSplit = 2,
                                   TOMType = "signed", networkType = "signed", minModuleSize = 100,
                                   reassignThreshold = 0.0001, mergeCutHeight = 0.25,
                                   numericLabels = TRUE, pamRespectsDendro = T,
                                   saveTOMs = F,
                                   verbose = 3)


####GO analysis####
# GO - analysis
all <- gogoanalysis(query_ID_den)
all <- filtering_GO(all)
for (i in names(query_ID_den)) {
  write.table(all[[i]], paste0(TABLEDIRECTORY,"GOall_", i, ".txt"), sep = "\t", quote = F, row.names = F)
}

# GO - analysis of WGCNA

mod_dat_GO_res <- list()
for (i in sort(unique(net_dat$colors))) {
  idx <- net_dat$colors == i
  tmp <- factor(as.integer(idx))
  names(tmp) <- rownames(y) #genelist#
  tgd1 <- new("topGOdata", ontology = "BP", allGenes = tmp, nodeSize = 5, annot = annFUN.org, mapping = "org.Dm.eg.db", ID = "ensembl") #data preparation#
  resTopGO.classic <- runTest(tgd1, algorithm = "classic", statistic = "Fisher") #enrichment test#
  resTopGO.weight01 <- runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res <- GenTable(tgd1, Fisher.classic = resTopGO.classic, Fisher.weight01 = resTopGO.weight01, orderBy = "Fisher.weight01", ranksOf = "Fisher.classic", topNodes = length(resTopGO.classic@score), numChar = 100) #analysis of results#
  mod_dat_GO_res[[i + 1]] <- tmp_res
}

names(mod_dat_GO_res) <- paste0("Module", 0:12)
lapply(mod_dat_GO_res, function(x) x[1:3, c(1:2, 8)])

#output the results
for (i in seq_along(mod_dat_GO_res)) {
  write.table(mod_dat_GO_res[[i]], paste0(TABLEDIRECTORY, names(mod_dat_GO_res)[i], "_GO_res.txt"), quote = F,
              row.names = F, sep = "\t")
}

#### VISUALISATION####
# PCA
mycol_pca <- c()
mycol_pca[density %in% "lowD"] <- "dark magenta"
mycol_pca[density %in% "HighD" & harv_time == "I"] <- "mediumpurple"
mycol_pca[density %in% "HighD" & harv_time == "II"] <- "tan3"
mycol_pca[density %in% "HighD" & harv_time == "III"] <- "steelblue"

mypch_pca <- c()
mypch_pca[density %in% "lowD"] <- 15
mypch_pca[density %in% "HighD" & harv_time == "I"] <- 16
mypch_pca[density %in% "HighD" & harv_time == "II"] <- 18
mypch_pca[density %in% "HighD" & harv_time == "III"] <- 17

# create the graphic
png(paste0(GRAPHDIRECTORY, "PCA.png"), width = 10, height = 10, units = "cm", res = 1080,pointsize = 8)
par(mar = c(2.5, 2.5, 0.5, 0.2), mgp = c(1.1, 0.3, 0), tcl = -0.15,
    las = 1, lwd = 3)
plot(pca$x, xlab = paste0("PC1 (", ve[1] * 100, "%)"), ylab = paste0("PC2 (", ve[2] * 100, "%)"),
     col = mycol_pca, asp = 1, pch = mypch_pca, cex = 1.5)
legend("topleft", pch = mypch_pca, col = mycol_pca, legend = c("lowD", "highD I", "highD II", "highD III"), cex = 1.5)
dev.off()

# create the Volcano plot
volcanoData <- as.data.frame(cbind(res_list[[1]]$logFC, res_list[[1]]$padj))
rownames(volcanoData) <- rownames(res_list[[1]])
colnames(volcanoData) <- c("logFC", "p.adjust")

png(paste0(GRAPHDIRECTORY, "Volcano.png"), width = 10, height = 10, units = "cm", res = 1080,pointsize = 8)
EnhancedVolcano(volcanoData, lab = rownames(volcanoData),
                x = 'logFC',
                y = 'p.adjust',
                pCutoff = 0.05,
                FCcutoff = log2(1.25),
                labSize = 0.0,
                title = "",
                subtitle = "",
                legendPosition = "top",
                borderWidth = 0.2,
                axisLabSize = 12,
                legendLabels = c(paste0("not significant (", nrow(volcanoData[volcanoData$logFC < log2(1.25) & volcanoData$p.adjust > 0.05,]), ")"),
                                 paste0("FC (", nrow(volcanoData[volcanoData$logFC > log2(1.25) & volcanoData$p.adjust > 0.05,]), ")"),
                                 paste0("significant (", nrow(volcanoData[volcanoData$logFC < log2(1.25) & volcanoData$p.adjust < 0.05,]), ")"),
                                 paste0("significant & FC (", nrow(volcanoData[abs(volcanoData$logFC) > log2(1.25) & volcanoData$p.adjust < 0.05,]), ")")),
                legendLabSize = 6
)
dev.off()

# create Venndiagram

query_ID_density_I_II <- create_query_ID_cat(query_ID_den, 1, query_ID_den[9:10], 9, query_ID_den[1:2], y, res_list$res_table_density, res_list$res_table_I_II)
query_ID_density_I_III <- create_query_ID_cat(query_ID_den, 1, query_ID_den[11:12], 11, query_ID_den[1:2], y, res_list$res_table_density, res_list$res_table_I_III)
query_ID_density_II_III <- create_query_ID_cat(query_ID_den, 1, query_ID_den[13:14], 13, query_ID_den[1:2], y, res_list$res_table_density, res_list$res_table_II_III)

venn.diagram(list(unlist(query_ID[3:4]), unlist(query_ID[5:6]), unlist(query_ID[7:8])),
             paste0(GRAPHDIRECTORY, "venndiagram.png"),
             width = 12, height = 12, units = "cm", res = 600, imagetype = "png", cat.cex = 1.25, cex = 1.5,
             category.names = c("highD I", "highD II", "highD III"), col = c("chartreuse3", "darkolivegreen4", "aquamarine4"),
             fill = c("chartreuse4", "darkolivegreen4", "aquamarine4"), lty = "blank")

# Expression of WGCNA modules/DE genes

scl_expr_dat = apply(expr_dat, 2, scale)
scl_avg_expr_dat = apply(scl_expr_dat, 2, function(x) tapply(x, paste(density, harv_time), mean))
scl_avg_expr_dat <- rbind(scl_avg_expr_dat[4,], scl_avg_expr_dat[1:3,])

png(paste0(GRAPHDIRECTORY, "expressionOFmodules.png"), width = 17.4, height = 8.7*1.25, units = "cm", res = 600, pointsize = 8)
par(mfrow = c(2, 3), mar = c(2, 2.5, 3, 2) + 0.1, oma = c(4,3,0,0) + 0.1, las = 1, lwd = 1)
for (i in c(1:5,12)) {
  idx <- net_dat$colors == i
  boxplot(t(scl_avg_expr_dat[, idx]), xlim = c(0.5, 4.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", xaxt = "n", boxlwd = 1.5,
          main = paste0("Module ", i, " - genes: ", sum(idx)), cex.main = 1.5, border = c("dark magenta", "mediumpurple", "tan3", "steelblue"), cex.axis = 1.1, col = "white")
  axis(1, at = 1:4, labels = c("lowD", "HighD I", "HighD II", "HighD III"), cex.axis = 1.1)
}
title(ylab = "normalized mean expression", outer = T, line = 0.5, cex.lab = 2)
dev.off()
png(paste0(GRAPHDIRECTORY, "expressionOFmodules_all.png"), width = 24, height = 16, units = "cm", res = 600, pointsize = 6)
par(mfrow = c(3, 4), mar = c(3, 2.5, 3, 2) + 0.1, oma = c(4,3,0,0) + 0.1, las = 1, lwd = 1)
for (i in 1:12) {
  idx <- net_dat$colors == i
  boxplot(t(scl_avg_expr_dat[, idx]), xlim = c(0.5, 4.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", xaxt = "n", boxlwd = 2.5,
          main = paste0("Module ", i, " - genes: ", sum(idx)), cex.main = 2.5, border = c("dark magenta", "mediumpurple", "tan3", "steelblue"), cex.axis = 1.25, col = "white")
  axis(1, at = 1:4, labels = c("lowD", "HighD I", "HighD II", "HighD III"), cex.axis = 1.5)
}
title(ylab = "normalized mean expression", outer = T, line = 0.5, cex.lab = 2.5)
dev.off()

#figure 2
png(paste0(GRAPHDIRECTORY,"fig2.png"),width = 4*3.5,height = 12,units = 'cm',pointsize = 9,res = 600)
layout(matrix(c(1,1,2,1,1,3),2,3,byrow = T),widths = c(1,1,1.5))
layout.show(3)
par(mar=c(5,5,4,1),las=1)
plot(pca$x, xlab = paste0("PC1 (", ve[1] * 100, "%)"), ylab = paste0("PC2 (", ve[2] * 100, "%)"),
     col = mycol_pca, asp = 1, pch = mypch_pca, cex = 1.5,cex.axis=1.25,cex.lab=1.25)
legend("topleft", pch = mypch_pca, col = mycol_pca, legend = c("lowD", "highD I", "highD II", "highD III"), cex = 1.5,bty = 'n')
par(mar=c(5,5,4,2),las=1)
boxplot(t(scl_avg_expr_dat[, query_ID_den$density_up]), xlim = c(0.5, 4.5), ylim = c(-1.5, 1.5), xlab = "", xaxt = "n", boxlwd = 1.5,
        main = paste0("Plastic up-regulated genes: ", length(query_ID_den$density_up)), border = c("dark magenta", "mediumpurple", "tan3", "steelblue"), 
        col = "white",ylab="Normalized expression",cex.lab=1.25)
axis(1, at = 1:4, labels = c("lowD", "HighD I", "HighD II", "HighD III"))
boxplot(t(scl_avg_expr_dat[, query_ID_den$density_dn]), xlim = c(0.5, 4.5), ylim = c(-1.5, 1.5), xlab = "", xaxt = "n", boxlwd = 1.5,
        main = paste0("Plastic down-regulated genes: ", length(query_ID_den$density_dn)), border = c("dark magenta", "mediumpurple", "tan3", "steelblue"),
        col = "white",ylab="Normalized expression",cex.lab=1.25)
axis(1, at = 1:4, labels = c("lowD", "HighD I", "HighD II", "HighD III"))
dev.off()
