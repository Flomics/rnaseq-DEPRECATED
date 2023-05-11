#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(stringr)
library(dplyr)
library(ggplot2)

# load concentration table and tpm matrix
ercc_conc <- read.csv(args[1]) #change this in pipeline with the path to bin
#ercc_conc <- read.csv("/home/ctuni/.nextflow/assets/Flomics/rnaseq/assets/ercc_concentration_table.csv")
ercc_len <- ercc_conc[,-2]
ercc_conc <- ercc_conc[,-3]
tpm_matrix <- read.delim(args[2], sep = "\t")
#tpm_matrix <- read.delim("/home/ctuni/Downloads/salmon.merged.gene_tpm.tsv", sep = "\t")

if (any(grep("ERCC-", tpm_matrix$gene_name, fixed = TRUE))) {
  if (grep("ERCC-", tpm_matrix$gene_name, fixed = TRUE)[1] != 0) {
    
    ##for tpm vs conc
    
    #remove all first column and all rows that are not ERCC
    ercc_tpm_matrix <- tpm_matrix %>%
      filter(str_detect(gene_id, "ERCC")) %>%
      select(-2)
    
    #alphabetically order both matrixes so they can be joined
    ercc_conc <- ercc_conc[order(ercc_conc$ERCC_ID), ]
    ercc_tpm_matrix <- ercc_tpm_matrix[order(ercc_tpm_matrix$gene_id), ]
    
    #crate final object to plot
    ercc_plot <- cbind(ercc_tpm_matrix, ercc_conc$conc_amoles_microl)
    
    colnames(ercc_plot) <- c(colnames(ercc_tpm_matrix), "conc")
    
    #plot in a loop and save as image
    ## remove first column
    ercc_plot <- ercc_plot %>%
      select(-1)
    
    i <- 1
    while (i < ncol(ercc_plot)) {
      ml <- lm(ercc_plot$conc ~ ercc_plot[, i])
      # Extracting R-squared parameter from summary
      r_squared <- summary(ml)$r.squared
      plot_title <- gsub("X", "", colnames(ercc_plot)[i])
      
      png(file = paste(gsub("X", "", colnames(ercc_plot)[i]), "_TPM_correlation_plot.png", sep = ""), width = 700, height = 700) #nolint
      print(ggplot(ercc_plot, aes(x = log10(ercc_plot$conc), y = log10(ercc_plot[, i] + 0.01))) + geom_point() + #nolint
              ylab("Spike-in TPMs (log10+0.01)") +
              xlab("Spike-in concentration [attomoles/microL] (log10)") +
              labs(title = plot_title) +
              geom_smooth(method = "lm", se = TRUE) +
              annotate(geom = "text", x = -1.3, y = max(log10(ercc_plot[, i] + 0.01)), label = paste("R² = ", r_squared, sep = "")) + #nolint
              theme_minimal())
      dev.off()
      i <- i + 1
    }
    
    #obtain correlation coefficients
    pearson_vector <- vector()
    spearman_vector <- vector()
    i <- 1
    while (i < ncol(ercc_plot)) {
      pearson_cor <- cor(ercc_plot$conc, ercc_plot[, i], method = "pearson")
      spearman_cor <- cor(ercc_plot$conc, ercc_plot[, i], method = "spearman")
      pearson_vector <- c(pearson_vector, pearson_cor)
      spearman_vector <- c(spearman_vector, spearman_cor)
      i <- i + 1
    }
    
    
    #create linear model to obtain R squared
    r_squared_vector <- vector()
    i <- 1
    while (i < ncol(ercc_plot)) {
      ml <- lm(ercc_plot$conc ~ ercc_plot[, i])
      # Extracting R-squared parameter from summary
      r_squared <- summary(ml)$r.squared
      r_squared_vector <- c(r_squared_vector, r_squared)
      i <- i + 1
    }
    
    
    #create small table with those three coefficients
    sample_names_vector <- colnames(ercc_plot)[-length(colnames(ercc_plot))]
    final_df <- data.frame(sample_names_vector, pearson_vector, spearman_vector, r_squared_vector) #nolint
    final_df$sample_names_vector <- gsub("X", "", final_df$sample_names_vector)
    colnames(final_df) <- c("Sample", "pearson_coef", "spearman_coef", "r_squared") #nolint
    write.table(final_df, file = "correlation_coefs.tsv", row.names = FALSE, quote = FALSE)
    
    
    ##for tmp vs len
    ercc_len <- ercc_len[order(ercc_len$ERCC_ID), ]
    ercc_plot_len <- cbind(ercc_tpm_matrix, ercc_len$length)
    
    colnames(ercc_plot_len) <- c(colnames(ercc_tpm_matrix), "length")
    
    #plot in a loop and save as image
    ## remove first column
    ercc_plot_len <- ercc_plot_len %>%
      select(-1)
    
    i <- 1
    while (i < ncol(ercc_plot_len)) {
      plot_title <- gsub("X", "", colnames(ercc_plot_len)[i])
      
      png(file = paste(gsub("X", "", colnames(ercc_plot_len)[i]), "_TPM_vs_length_correlation_plot.png", sep = ""), width = 700, height = 700) #nolint
      print(ggplot(ercc_plot_len, aes(x = ercc_plot_len$length, y = ercc_plot_len[, i])) + geom_point() + #nolint
              ylab("Spike-in TPMs") +
              xlab("Spike-in length [nt]") +
              labs(title = plot_title) +
              theme_minimal())
      dev.off()
      i <- i + 1
    }
    
  } else {
    
    ercc_tpm_matrix <- tpm_matrix %>%
      select(c(-1, -2))
    n_samples <- ncol(ercc_tpm_matrix)
    
    final_matrix <- matrix(data = "no_spike_ins", nrow = n_samples, ncol = 3)
    sample_names_vector <- colnames(ercc_tpm_matrix)
    final_df <- data.frame(sample_names_vector, final_matrix)
    final_df$sample_names_vector <- gsub("X", "", final_df$sample_names_vector)
    colnames(final_df) <- c("sample", "pearson_coef", "spearman_coef", "r_squared") #nolint
    write.table(final_df, file = "correlation_coefs.tsv", row.names = FALSE)
    
  }
} else {
  ercc_tpm_matrix <- tpm_matrix %>%
    select(c(-1, -2))
  n_samples <- ncol(ercc_tpm_matrix)
  
  final_matrix <- matrix(data = "no_spike_ins", nrow = n_samples, ncol = 3)
  sample_names_vector <- colnames(ercc_tpm_matrix)
  final_df <- data.frame(sample_names_vector, final_matrix)
  final_df$sample_names_vector <- gsub("X", "", final_df$sample_names_vector)
  colnames(final_df) <- c("sample", "pearson_coef", "spearman_coef", "r_squared") #nolint
  write.table(final_df, file = "correlation_coefs.tsv", row.names = FALSE)
}

