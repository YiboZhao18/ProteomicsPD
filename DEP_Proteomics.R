# DEP interpretation and visualisation
# 1. general visualisation
## 1.1 PCA
PDMerge <- read.delim("Proteomics_input.txt", stringsAsFactors = F)
PDMerge_values <- PDMerge[,-1]
row.names(PDMerge_values) <- PDMerge$NPID
PID <- PDMerge$NPID
pheno <- matrix(unlist(strsplit(colnames(PDMerge_values), split = "_")), ncol = 2, byrow = T)[,2]

# install.packages("factoextra")
# library(factoextra)
# install.packages("FactoMineR")
# library(FactoMineR)
# library(cluster)
# install.packages("ggfortify")
library(ggfortify)
# PCA need subjects as rows and probes as columns
tPDMerge_values <- as.data.frame(t(PDMerge_values))
colnames(tPDMerge_values) <- PID
mean <- sapply(tPDMerge_values, function(x) mean(x))
# head(mean)
Probe_sort <- data.frame(PID, mean)
Probe_sort <- Probe_sort[order(Probe_sort$mean, decreasing = T),]
top500Probe <- Probe_sort$PID[1:500] # PCA on top 500 probes
tPDMerge_values_top500 <- tPDMerge_values[,top500Probe]
tPDMerge_values_top500$pheno <- pheno
PDMerge.pca <- prcomp(tPDMerge_values_top500[,1:500], scale. = T)
autoplot(PDMerge.pca, data = tPDMerge_values_top500, colour = "pheno")
# with top 500 is unable to classify the phenotypes, try again later with DEPs.

## 1.2 Correlation matrix across samples
# install.packages("ggplot2")
library(ggplot2)
# install.packages("corrplot")
library(corrplot)
cor_matrix <- abs(cor(PDMerge_values))
corrplot(cor_matrix, method = "color", 
         col = colorRampPalette(c("lightpink", "darkred"))(200),
         tl.cex = 0.8,  # Text label size
         addCoef.col = "black",  # Add correlation coefficients on plot
         number.cex = 0.7,  # Coefficient text size
         tl.srt = 45)  # Text label rotation
# every sample is very similar to each other, hard to differntiate in this way. shows no batch effect at least...


## 2. DEP-related graphs

get_DEP_res <- function(group1, group2){
        file.name <- paste(group1, "vs", group2, ".txt", sep = "")
        res <- read.delim(file.name, stringsAsFactors = F)
        res <- res[order(res$p),]
        p.names <- colnames(res)[2:3]
        colnames(res)[2:3] <- paste(paste(group1, "vs", group2, sep = ""), p.names, sep = "_")
        return(res)
}

get_sig_probe <- function(group1, group2){
        file.name <- paste(group1, "vs", group2, ".txt", sep = "")
        res <- read.delim(file.name, stringsAsFactors = F)
        sig.list <- res$PID[res$adj.p<0.05]
        return(sig.list)
}

# 2.1. heatmap for all significant proteins
group <- sort(unique(pheno))
all_sig_probe <- integer()
for (i in 1:length(group)){
        for (j in i:length(group)){
                if (group[i] == group[j]){
                        flag = "same pheno, move on"
                } else{
                        all_sig_probe <- c(all_sig_probe, get_sig_probe(group[i], group[j]))
                }
        }
}
all_sig_probe <- unique(all_sig_probe)
PDMerge_values_sig <- tPDMerge_values[,all_sig_probe]
res1 <- get_DEP_res("HC", "LPD")
res2 <- get_DEP_res("HC", "LProdromal")
res3 <- get_DEP_res("HC", "sPD")
res4 <- get_DEP_res("LPD", "LProdromal")
res5 <- get_DEP_res("LPD", "sPD")
res6 <- get_DEP_res("LProdromal", "sPD")
res_list <- list(res1, res2, res3, res4, res5, res6)
all_res <- Reduce(function(x,y) merge(x,y, by = "PID", all = T), res_list)
# colnames(all_res)
# keep only mean for each column
row.names(all_res) <- all_res$PID
column_to_remove <- integer()
i = 1
for (i in 1:ncol(all_res)){
        temp <- colnames(all_res)[i]
        flag1 <- unlist(strsplit(temp, split = ".", fixed = T))[2] == "x"
        flag2 <- unlist(strsplit(temp, split = ".", fixed = T))[2] == "y"
        if (is.na(flag1) & is.na(flag2)){
                message <- "this is not a duplicated column"
        } else if (flag1 | flag2){
                message <- "this is a duplicated column"
                column_to_remove <- c(column_to_remove, i)
        }
}
all_res <- all_res[,-column_to_remove]
write.table(all_res, "DEP_results_all.txt", sep = "\t", quote = F, row.names = F)
all_res.sig <- all_res.sig[all_sig_probe,]
# install.packages("pheatmap")
library("pheatmap")
## 2.1.1 heatmap for test values of all pheno
all_res_nop <- all_res.sig[,group]
pheatmap(all_res_nop, 
         cluster_rows = TRUE,  # Cluster the rows
         cluster_cols = FALSE,  # Cluster the columns
         display_numbers = TRUE,  # Show the values of the heatmap
         show_rownames = F,  # Adjust row font size
         cutree_rows = 5,  # Cut the dendrogram into 4 clusters
         fontsize_number = 8,  # Font size for numbers on heatmap cells
         color = colorRampPalette(c("mistyrose", "darkred"))(20)  # Customize the color palette
)

## 2.1.2 make a heatmap with everything compared to HC
all_res_vsHC <- data.frame(LPD = all_res_nop$LPD-all_res$HC, LProdromal = all_res_nop$LProdromal-all_res$HC, sPD = all_res_nop$sPD-all_res$HC)
rownames(all_res_vsHC) <- all_sig_probe
comp_colours <- colorRampPalette(c("forestgreen","white", "firebrick"))(20)
# to make sure that above 0 are red and below are green, we need to make the breaks smaller
comp_breaks <- c(seq(min(all_res_vsHC), -0.05, length.out = 10), 
                 seq(0, 0.05, length.out = 1), 
                 seq(0.05, max(all_res_vsHC), length.out = 10))
all_res_vsHC.heatmap <- pheatmap(all_res_vsHC, 
                                 cluster_rows = TRUE,  # Cluster the rows
                                 cluster_cols = FALSE,  # Cluster the columns
                                 display_numbers = TRUE,  # Show the values of the heatmap
                                 show_rownames = F,  # Adjust row font size
                                 cutree_rows = 13,  # Cut the dendrogram into 4 clusters
                                 fontsize_number = 8,  # Font size for numbers on heatmap cells
                                 breaks = comp_breaks,
                                 color = comp_colours  # Customize the color palette
)

row_clusters <- cutree(all_res_vsHC.heatmap$tree_row, k = 13)
sep_PID.sig <- as.data.frame(matrix(unlist(strsplit(names(row_clusters), split = ":")), ncol = 3, byrow = T))
colnames(sep_PID.sig) <- c("ProbeID", "Ensebl", "GeneSymbol")
PDMerge_cluster <- data.frame(sep_PID.sig, cluster = row_clusters, all_res_vsHC)
# head(PDMerge_Cluster)
PDMerge_cluster <- PDMerge_cluster[order(PDMerge_cluster$cluster),]
write.table(PDMerge_cluster, "cluster_88sig.txt", sep = "\t", quote = F, row.names = F)

# 2.2 PCA plot on significant genes
tPDMerge_values_sig <- tPDMerge_values[,all_sig_probe]
tPDMerge_values_sig.pca <- prcomp(tPDMerge_values_sig, scale. = T)
tPDMerge_values_sig$pheno <- pheno
autoplot(tPDMerge_values_sig.pca, data = tPDMerge_values_sig, colour = "pheno")
# NOTE: still unable to differetiate the 4 pheno types

# 2.3 Volcano plots
get_volcano_plot <- function(group1, group2){
        all_res <- read.delim("DEP_results_all.txt", stringsAsFactors = F)
        M1 <- all_res[, match(group1, colnames(all_res))]
        M2 <- all_res[, match(group2, colnames(all_res))]
        log2FC <- M2-M1
        pair <- paste(group1, "vs", group2, sep = "")
        NegativeLogp <- -log10(all_res[,match(paste(pair,"p", sep = "_"), colnames(all_res))])
        adj.p <- all_res[,match(paste(pair,"adj.p", sep = "_"), colnames(all_res))]
        is_sig <- rep("non-significant", length(log2FC))
        is_sig[adj.p<0.05&log2FC>0] <- "up-regulated"
        is_sig[adj.p<0.05&log2FC<0] <- "down-regulated"
        data <- data.frame(log2FC, NegativeLogp, is_sig, adj.p)
        if (length(unique(is_sig)) == 3){
                message <- "all 3 groups are here"
                plot <- ggplot(data=data, aes(x=log2FC, y=NegativeLogp, colour = is_sig)) + 
                        geom_point() + 
                        theme_minimal() +
                        # Add lines as before...
                        geom_vline(xintercept = 0, col = "grey", linetype = "dashed")+
                        # geom_vline(xintercept=c(-0.1, 0.1), col="grey", linetype = "dashed") +
                        geom_hline(yintercept=-log10(0.05/length(log2FC)), col="grey", linetype = "dashed") +
                        ## Change point color 
                        scale_color_manual(values=c("steelblue", "grey", "indianred")) +
                        theme(legend.title = element_blank())+
                        labs(title = paste(group1, "vs", group2, sep = " "), y = "-log10(p)") 
        } else if (length(unique(is_sig)) == 2){
                message <- "only 2 groups appear"
                if(is.na(match("up-regulated", is_sig))){
                        message <- "no up-regulated"
                        plot <- ggplot(data=data, aes(x=log2FC, y=NegativeLogp, colour = is_sig)) + 
                                geom_point() + 
                                theme_minimal() +
                                # Add lines as before...
                                geom_vline(xintercept = 0, col = "grey", linetype = "dashed")+
                                # geom_vline(xintercept=c(-0.1, 0.1), col="grey", linetype = "dashed") +
                                geom_hline(yintercept=-log10(0.05/length(log2FC)), col="grey", linetype = "dashed") +
                                ## Change point color 
                                scale_color_manual(values=c("grey", "steelblue")) +
                                theme(legend.title = element_blank())+
                                labs(title = paste(group1, "vs", group2, sep = " "), y = "-log10(p)") 
                } else{
                        message <- "no down-regulated"
                        plot <- ggplot(data=data, aes(x=log2FC, y=NegativeLogp, colour = is_sig)) + 
                                geom_point() + 
                                theme_minimal() +
                                # Add lines as before...
                                geom_vline(xintercept = 0, col = "grey", linetype = "dashed")+
                                # geom_vline(xintercept=c(-0.1, 0.1), col="grey", linetype = "dashed") +
                                geom_hline(yintercept=-log10(0.05/length(log2FC)), col="grey", linetype = "dashed") +
                                ## Change point color 
                                scale_color_manual(values=c("grey", "indianred")) +
                                theme(legend.title = element_blank())+
                                labs(title = paste(group1, "vs", group2, sep = " "), y = "-log10(p)") 
                }
        } else{
                message <- "no significant results!"
                # Create a data frame for the eyes
                eyes_data <- data.frame(x = c(-0.3, 0.3), y = c(0.4, 0.4), r = c(0.1, 0.1))   # Radius for the eyes
                # Create a data frame for the mouth
                x_values <- seq(-0.5, 0.5, length.out = 100)
                mouth_data <- data.frame(x = x_values, y = -0.5 * (x_values^2) - 0.5)   # A parabolic shape for a sad mouth
                # Plot the sad face
                plot <- ggplot() +
                        
                        # Draw the eyes (two small circles)
                        geom_point(data = eyes_data, aes(x, y), size = 5, shape = 21, fill = "black") +
                        
                        # Draw the mouth (a sad curve)
                        geom_line(data = mouth_data, aes(x, y), color = "black", size = 1) +
                        
                        # Set equal aspect ratio to make the face circular
                        coord_fixed() +
                        
                        # Set limits and theme
                        xlim(-1.5, 1.5) +
                        ylim(-1.5, 1) +
                        theme_void() +
                        theme(panel.background = element_rect(fill = "white")) + # Set background to white
                        labs(title = "no significant results...")
                
                
        }
        return(plot)
}
p1 <- get_volcano_plot("HC", "LPD")
p2 <- get_volcano_plot("HC", "LProdromal")
p3 <- get_volcano_plot("HC", "sPD")
p4 <- get_volcano_plot("LPD", "LProdromal") # one group is missing
p5 <- get_volcano_plot("LPD", "sPD")
p6 <- get_volcano_plot("LProdromal", "sPD")

p1
p2
p3
p4
p5
p6

