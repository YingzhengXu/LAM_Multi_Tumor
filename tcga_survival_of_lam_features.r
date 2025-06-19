# Load libraries
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.rnaseq)
library(survival)
library(survminer)
library(ggpubr)

# Load clinical + RNA-seq data
data(BRCA.clinical)
data(BRCA.rnaseq)

# Prepare survival data
survival_data <- survivalTCGA(BRCA.clinical, extract.cols = c("days_to_death", "vital_status"))
rownames(survival_data) <- survival_data$bcr_patient_barcode
survival_data <- na.omit(survival_data)

# Prepare RNA-seq data
rnaseq_data <- as.data.frame(BRCA.rnaseq)
colnames(rnaseq_data) <- sub("\\|.*", "", colnames(rnaseq_data)) # Remove ENSG suffix
rownames(rnaseq_data) <- substr(rnaseq_data$bcr_patient_barcode, 1, 12)
rnaseq_data <- rnaseq_data[!duplicated(rownames(rnaseq_data)), ]

# Match sample IDs
common_samples <- intersect(rownames(rnaseq_data), rownames(survival_data))
rnaseq_data <- rnaseq_data[common_samples, ]
survival_data <- survival_data[common_samples, ]

# Genes of interest (LAM signature or AUCell top TFs)
genes_list <- c("APOC1", "LIPA", "SPP1", "TREM2", "CD9", "LPL", "APOE", "GPNMB", "PLD3")

# Survival analysis loop
results <- list()
plot_list <- list()

for (gene in genes_list) {
  if (!gene %in% colnames(rnaseq_data)) {
    cat(sprintf("Gene %s not found in data.\n", gene))
    next
  }
  
  # Expression and group
  expr <- as.numeric(rnaseq_data[, gene])
  group <- ifelse(expr > median(expr, na.rm=TRUE), "High", "Low")
  
  survival_data$group <- group
  
  # Cox model / KM
  surv_obj <- Surv(time = survival_data$times, event = survival_data$patient.vital_status)
  fit <- survfit(surv_obj ~ group, data = survival_data)
  
  # Log-rank p-value
  logrank <- survdiff(surv_obj ~ group, data = survival_data)
  pval <- 1 - pchisq(logrank$chisq, df = 1)
  
  # Save result + plot
  results[[gene]] <- list(fit=fit, pval=pval)
  plot_list[[gene]] <- ggsurvplot(fit, data=survival_data, pval=TRUE, 
                                  title=gene, risk.table=TRUE)$plot
}

# Combine plots
combined_plot <- ggarrange(plotlist=plot_list, ncol=3, nrow=3, common.legend=TRUE)
print(combined_plot)

# Save PDF
pdf("TCGA_LAM_gene_survival.pdf", width=12, height=10)
print(combined_plot)
dev.off()

# Summary p-values
for (gene in names(results)) {
  cat(sprintf("Gene: %s, p-value: %.4g\n", gene, results[[gene]]$pval))
}
