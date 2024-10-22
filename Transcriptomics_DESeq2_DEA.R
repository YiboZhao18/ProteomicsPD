
# DEA on 4 time points 
PhenoData <- read.delim("pheno.txt", stringsAsFactors = F)
pheno <- integer()
pheno[PhenoData$pheno==0] <- "HC"
pheno[PhenoData$pheno==1 & PhenoData$S.G==0] <- "sPD"
pheno[PhenoData$pheno==1 & PhenoData$L==1] <- "LPD"
pheno[PhenoData$pheno==2 & PhenoData$L==1] <- "LProdromal"
pheno_table <- data.frame(patho = PhenoData$ID, pheno = pheno)
pheno_table <- as.data.frame(na.omit(pheno_table))
write.table(pheno_table, "pheno_table.txt", sep = "\t", row.names = F, quote = F)


# load rc matrix
get_DEA <- function(infile, group1, group2){
        library(DESeq2)
        # read count data
        rcdata <- read.delim(infile, stringsAsFactors = F)
        rcdata<- na.omit(rcdata)
        count_matrix <- rcdata[,-c(1)]
        genename <- rcdata$protein
        rownames(count_matrix) <- genename
        ## in the .txt file, colnames have "X" before patno
        
        # read pheno data
        pheno_table <- read.delim("pheno_table.txt")
        patno <- pheno_table$patho[pheno_table$pheno==group1|pheno_table$pheno==group2]
        # extract 
        patno_in <- intersect(paste("X", patno, sep = ""), colnames(count_matrix))
        pheno_in <- pheno_table$pheno[match(patno_in, paste("X", pheno_table$patho, sep = ""))]
        count_matrix_trimmed <- count_matrix[,patno_in]
        condition_matrix <- as.data.frame(cbind(patno_in, pheno_in))
        
        # DESeq2
        sample75 <- 0.75*dim(count_matrix_trimmed)[2]
        lowcount <- row.names(count_matrix_trimmed)[!rowSums(count_matrix_trimmed >=15) >= sample75]
        gene_in <- row.names(count_matrix_trimmed)[rowSums(count_matrix_trimmed >=15) >= sample75]
        count_matrix_trimmed_in  <- count_matrix_trimmed[rowSums(count_matrix_trimmed >=15) >= sample75,]
        colnames(condition_matrix) <- c("patno", "pheno")
        condition_matrix <- as.data.frame(condition_matrix)
        condition_matrix$pheno <- as.factor(condition_matrix$pheno)
        rownames(count_matrix_trimmed_in) <- gene_in
        # generate a normalised count matrix without design for WGCNA
        dds <- DESeqDataSetFromMatrix(count_matrix_trimmed_in, condition_matrix, design = ~pheno)
        dds <- DESeq(dds)
        # dds <- estimateSizeFactors(dds)
        norm.counts <- counts(dds, normalized = T)
        res <- results(dds, alpha = 0.05) # 1 is the reference level so the res is case vs. control
        p <- res$pvalue
        padj <- p.adjust(p, method = "BH")
        # write.table(norm.counts, paste("norm",infile, sep = "_"), sep = "\t", quote = F, row.names = T)
        outfile.name <- paste(paste("DEA", infile, group1, group2, sep = "_"), "txt", sep = ".")
        write.table(cbind(res$log2FoldChange, res$pvalue, padj), outfile.name, sep = "\t", quote = F, row.names = gene_in, col.names = c("log2FC", "pvalue", "padj"))
}

get_DEA("rc_BL.txt", "HC", "sPD")
get_DEA("rc_BL.txt", "HC", "LPD")
get_DEA("rc_BL.txt", "HC", "LProdromal")
get_DEA("rc_BL.txt", "LPD", "LProdromal")
get_DEA("rc_BL.txt", "LProdromal", "sPD")
get_DEA("rc_BL.txt", "LPD", "sPD")

