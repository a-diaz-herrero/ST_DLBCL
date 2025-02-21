###title: "Cell-Eco signature scores"
## Created by: Alba Diaz Herrero

# Inputs: matrix of survival with RNAseq data and clinical from Lacy et al. 2020
# Outputs: Matrix of cell-type proportions per spot



# Load librairies ---------------------------------------------------------

library(patchwork)
library(umap)
library(dplyr) 
library(ggpubr) 
library(Seurat)
library(ggplot2)
library(ggraph)
library("survival")
library("survminer")
library(clustree)
library(writexl)
library(VennDiagram)
library(grid)
set.seed(123)


# Directories-----------------------------------------------------------------

dir_DEG<-"path " # path to folder of differentially expressed genes
dir_bulk_plots<-"path/plots"

# Set up----------------------------------------------------------------
# Upload necessary files:
DEG<- as.data.frame(readRDS("path")) #path to the list of differentially expressed genes between celleco of integrated spatial objects
matrix_sur<-readRDS("path") #being the matrix of survival with gene expression and clinica data
matrix_name<-"Lacy" #or any name
# Palette
pal <- c("darkgreen", "darkorange")

#  Initialize storage for results
genes_mat <- data.frame(CellEco = character(), Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)
survival_plots <- list()
tertiles_CellEco_list<- list()
histograms <- list()
updated_matrix_list <- list()

# Signature score calculation----------------------------------------------------------------

# Loop through each CellEco value (1 to 6)
for (celleco in 1:6) {
  
  # Filter data for the specific CellEco cluster
  topA_pval <- as.data.frame(DEG[DEG$cluster == as.character(celleco), ])
  rownames(topA_pval) <- topA_pval$gene
  
  # Ensure the p-value data frame is not empty
  if (nrow(topA_pval) == 0) {
    cat("No genes significant for CellEco", celleco, "in Lacy et al. 2020")
    next  # Skip to the next iteration
  }
  
  # Keep all significant genes lower than the 5th percentile
  percentile_5 <- quantile(topA_pval$p_val_adj, probs = 0.05)
  topA_pval <- topA_pval[topA_pval[, "p_val_adj"] < percentile_5, ]
  
  # Get the intersection of genes in signature and bulk matrix
  genes <- intersect(topA_pval$gene, colnames(matrix_sur))
  
  # Extract relevant data for the current CellEco
  current_genes <- data.frame(
    CellEco = rep(paste0("CellEco", celleco), length(genes)),  # Assign CellEco label
    Gene = genes,  # Gene names
    P_Value = topA_pval[genes, "p_val_adj"],  # Corresponding p-values
    stringsAsFactors = FALSE
  )
  
  # Combine with the main genes_mat table
  genes_mat <- rbind(genes_mat, current_genes)
  
  # Handle p-values and calculate scores
  topA_pval$p_val_adj <- as.numeric(as.character(topA_pval$p_val_adj))
  non0pval <- topA_pval[which(topA_pval$p_val_adj > 0), ]
  min_p_val <- min(non0pval$p_val_adj)
  
  # Process each gene
  for (gene in genes) {
    p_val <- topA_pval[gene, "p_val_adj"]
    if (p_val == 0) {
      p_val <- min_p_val  # Substitute 0 with the minimum non-zero p-value
    }
    
    if (p_val > 0) {
      matrix_sur[[paste0(gene, "_pval")]] <- -log10(p_val)
    } else {
      matrix_sur[[paste0(gene, "_pval")]] <- NA
    }
    
    matrix_sur[[paste0(gene, "_sign")]] <- ifelse(topA_pval[gene, "diffexpressed"] == "Upregulated", 1, -1)
    
    matrix_sur[[paste0(gene, "_score_CellEco", celleco)]] <- 
      matrix_sur[[paste0(gene, "_pval")]] * 
      matrix_sur[[gene]] * 
      matrix_sur[[paste0(gene, "_sign")]]
  }
  
  # Summarize scores for the CellEco signature
  score_cols <- grep(paste0("_score_CellEco", celleco), colnames(matrix_sur), value = TRUE)
  matrix_sur[paste0("signature_score_CellEco", celleco)] <- rowSums(matrix_sur[, score_cols, drop = FALSE], na.rm = TRUE)
  
  # Annotate patients based on tertiles
  tertiles <- quantile(matrix_sur[, paste0("signature_score_CellEco", celleco)], probs = c(1/3, 2/3))
  tertiles_CellEco <- rep(NA, nrow(matrix_sur))
  
  higher_indices <- which(matrix_sur[[paste0("signature_score_CellEco", celleco)]] > tertiles[2])
  lower_indices <- which(matrix_sur[[paste0("signature_score_CellEco", celleco)]] < tertiles[1])
  tertiles_CellEco[higher_indices] <- "high"
  tertiles_CellEco[lower_indices] <- "low"
  
  # Save tertiles to the list for the current CellEco
  tertiles_CellEco_list[[paste0("CellEco", celleco)]] <- tertiles_CellEco
  
  # Perform survival analysis and plot
  matrix_sur[[paste0("tertiles_CellEco", celleco)]] <- tertiles_CellEco  # Add annotation column
  fit <- survfit(Surv(time, status) ~ matrix_sur[[paste0("tertiles_CellEco", celleco)]], data = matrix_sur)
  # Plots----------------------------------------------------------------
    survival_plot <- ggsurvplot(fit, 
                              pval = TRUE,
                              palette = pal,
                              title = paste0("Survival Analysis for CellEco ", celleco),
                              xlab = "Time (months)",
                              ylab = "Survival Probability",
                              legend="right",
                              legend.title=paste0("signature score CellEco ", celleco),
                              subtitle = paste0(matrix_name, " et al."),
                              legend.labs=c("high","low"))
  
  # Store the plot in the list
  survival_plots[[paste0(matrix_name, "_Survival_CellEco", celleco, "_p5percentile_du")]] <- survival_plot
  # Generate and store histogram for the scores
  hist_plot <- ggplot(matrix_sur, aes(x = .data[[paste0("signature_score_CellEco", celleco)]])) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    geom_vline(xintercept = tertiles, color = "red", linetype = "dashed", size = 1) +  # Add vertical lines for tertiles
    ggtitle(paste("Histogram of Signature Scores for CellEco", celleco)) +
    xlab("Signature Score") +
    ylab("Frequency")+
    labs(subtitle = paste0(matrix_name," et al."))
  histograms[[paste0(matrix_name, "_Histogram_CellEco", celleco)]] <- hist_plot
}

# Save the updated matrix with scores and other results in the new list
updated_matrix_list[[matrix_name]] <- matrix_sur

# Acess to plots:
celleco<-"1"
histo<- paste0(matrix_name, "_Histogram_CellEco", celleco)
survival_plot <- survival_plots[[paste0(matrix_name, "_Survival_CellEco", celleco, "_p5percentile_du")]]
survival_plot

# Cox HR ----------------------------------------------------------------

# Replace "low" with "alow" for HR reference simply to have it alfabetically ordered, and have low as reference:
tertiles_columns <- grep("^tertiles", colnames(matrix_sur), value = TRUE)
for (col in tertiles_columns) {
  matrix_sur[[paste0("hr_", col)]] <- ifelse(matrix_sur[[col]] == "low", "alow", matrix_sur[[col]])
}

# Initialize Cox results list
cox_results <- data.frame(
  HR = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric(),
  p_value = numeric(),
  CellEco = numeric(),
  stringsAsFactors = FALSE
)

# Loop through CellEco clusters for Cox analysis
for (celleco in 1:6) {
  tertiles_col <- paste0("hr_tertiles_CellEco", celleco)
  if (!tertiles_col %in% colnames(matrix_sur)) next
  
  formula <- as.formula(paste("Surv(time, status) ~", tertiles_col))
  res.cox <- tryCatch({ coxph(formula, data = matrix_sur) }, error = function(e) NULL)
  
  if (!is.null(res.cox)) {
    summary_cox <- summary(res.cox)
    coef_name <- paste0(tertiles_col, "high")
    
    if (coef_name %in% rownames(summary_cox$coefficients)) {
      cox_results <- rbind(cox_results, data.frame(
        HR = summary_cox$coefficients[coef_name, "exp(coef)"],
        LowerCI = summary_cox$conf.int[coef_name, "lower .95"],
        UpperCI = summary_cox$conf.int[coef_name, "upper .95"],
        p_value = summary_cox$coefficients[coef_name, "Pr(>|z|)"],
        CellEco = celleco
      ))
    }
  }
}

# Add log2 HR & significance
cox_results$log2HR <- log2(cox_results$HR)
cox_results$Significance <- ifelse(
  cox_results$p_value < 0.001, "***",
  ifelse(cox_results$p_value < 0.01, "**",
         ifelse(cox_results$p_value < 0.05, "*", " ")))

# Plot Cox HR results ----------------------------------------------------------------

ggplot(cox_results, aes(x = log2HR, y = factor(CellEco, levels = rev(unique(CellEco))))) +  
  geom_errorbarh(aes(xmin = log2(LowerCI), xmax = log2(UpperCI)), height = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_point(aes(color = ifelse(p_value < 0.05, "significant", "non-significant")), size = 3) +
  geom_text(aes(label = paste0("p=", format(p_value, digits = 2))), vjust = 2) +
  scale_color_manual(values = c("non-significant" = "gray", "significant" = "black")) +
  labs(
    title = "Hazard Ratios for CellEco Clusters",
    x = "Log2 Hazard Ratio",
    y = "CellEco Cluster"
  ) + 
  theme_minimal()



