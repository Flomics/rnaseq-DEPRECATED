#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(stringr)
library(dplyr)
library(ggplot2)

# load concentration table and tpm matrix
ercc_conc <- read.csv(args[1]) #change this in pipeline with the path to bin
tpm_matrix <- read.delim(args[2], sep = "\t")
#tpm_matrix <- read.delim("/home/ctuni/Flomics/spike-in_qc/salmon.merged.gene_tpm.tsv", sep = "\t")


#remove all first column and all rows that are not ERCC
ercc_tpm_matrix <- tpm_matrix %>%
   filter(str_detect(gene_id, "ERCC")) %>%
   select(-2)

#alphabetically order both matrixes so they can be joined
ercc_conc <- ercc_conc[order(ercc_conc$ERCC_ID),]
ercc_tpm_matrix <- ercc_tpm_matrix[order(ercc_tpm_matrix$gene_id),]

#crate final object to plot
ercc_plot <- cbind(ercc_tpm_matrix, ercc_conc$conc_amoles_microl)

colnames(ercc_plot) <- c(colnames(ercc_tpm_matrix), "conc")

#plot in a loop and save as image
## remove first column
ercc_plot <- ercc_plot %>%
  select(-1)

i=1
while (i < ncol(ercc_plot)) {
  #print(length(ercc_plot$conc)==length(ercc_plot[,i]))
  png(file=paste(gsub("X", "",colnames(ercc_plot)[i]), "_TPM_correlation_plot.png", sep = ""))
  print(ggplot(ercc_plot, aes(x=ercc_plot$conc, y=ercc_plot[,i])) + geom_point() +
    ylab("Spike-in TPMs") +
    xlab("Spike-in concentration [attomoles/μL]") +
    theme_minimal())
  #plot(x=ercc_plot$conc, y=ercc_plot[,i], main = gsub("X", "",colnames(ercc_plot)[i]), xlab = "Spike-in concentration [attomoles/μL]", ylab = "Spike-in TPMs")
  dev.off()
  i <- i + 1
}

#obtain correlation coeficients
pearson_vector <- vector()
spearman_vector <- vector()
i=1
while (i < ncol(ercc_plot)) {
  pearson_cor <- cor(ercc_plot$conc, ercc_plot[,i], method = "pearson")
  spearman_cor <-cor(ercc_plot$conc, ercc_plot[,i], method = "spearman")
  pearson_vector <- c(pearson_vector, pearson_cor)
  spearman_vector <- c(spearman_vector, spearman_cor)
  i <- i + 1
}


#create linear model to obtain R squared
r_squared_vector <- vector()
i=1
while (i < ncol(ercc_plot)) {
  ml <- lm(ercc_plot$conc~ercc_plot[,i])
  # Extracting R-squared parameter from summary
  r_squared <- summary(ml)$r.squared
  r_squared_vector <- c(r_squared_vector, r_squared)
  i <- i + 1
}


#create small table with those three coefficients
sample_names_vector <- colnames(ercc_plot)[-length(colnames(ercc_plot))]
final_df <- data.frame(sample_names_vector,pearson_vector, spearman_vector, r_squared_vector)
final_df$sample_names_vector <- gsub("X", "", final_df$sample_names_vector)
colnames(final_df) <- c("sample", "pearson_coef", "spearman_coef", "r_squared")
write.table(final_df, file = "correlation_coefs.tsv", row.names = FALSE)
