# This is a code for predicting the five CRCAssigner subtypes by correlating the expression of each sample 
# to the CRCAssigner PAM centroids (786 genes). The sample is assigned into the subtype with the maximun 
# Pearson correlation.
# direc: is the results folder
# file: is the file name of the expression data with the samples to be predicted
# PAM: is the 786 genes CRCAssigner PAM centroids
# und_cutoff: cutoff of identifying low confidence samples (und_cutoff = 0.15)
# mix_cutoff: cutoff of identifying mixed samples (mix_cutoff = 0.06)
 
correl_subtypes<- function(direc, file, PAM, und_cutoff, mix_cutoff)
{
  
  setwd( direc )
  file_s <- strsplit(file, ".txt")
  
  outfile <- paste0( file_s, "_CRCassigner_map.txt", sep ="" )
  outfile_cor <- paste0( file_s, "_correlation.txt", sep ="" )
  outfile_tab <- paste0( file_s, "_subtypes_table.txt", sep ="" )
  outfile_file <- paste0( file_s, "_CRCAssigner_786_data.txt", sep ="" )
  
  # reading data
  data <- read.delim(file, stringsAsFactors = FALSE)
  
  # making a matrix
  data_m <- data.matrix(data[,2:ncol(data)])
  row.names(data_m) <- data[,1]
  
  # median centering data
  med <- apply(data_m, 1, median)
  data_m_med <- data_m - med
  
  # reading PAM centroids
  crc <- read.delim(PAM, stringsAsFactors = FALSE)
  
  # making a matrix
  crc_m <- data.matrix(crc[,2:ncol(crc)])
  row.names(crc_m) <- crc[,1]
  
  # matching genes between data and PAM centroids
  m <- match(rownames(crc_m), rownames(data_m))
  w <- which(!is.na(m))
  
  # performing correlation
  corr <- cor(data_m_med[m[w],], crc_m[w,])
  
  # finding which subtype has maximum correlation, then first, second, third, fourth and fifth subtypes
  corr_max_w <- apply(corr, 1, function(x) which.max(x))
  # corr_sort <- apply(corr, 1, function(x) sort(x,partial=length(x)-1)[length(x)-1])
  corr_second <- apply(corr, 1, function(x) which(x == sort(x,partial=length(x)-1)[length(x)-1]))
  corr_third <- apply(corr, 1, function(x) which(x == sort(x,partial=length(x)-2)[length(x)-2]))
  corr_fourth <- apply(corr, 1, function(x) which(x == sort(x,partial=length(x)-3)[length(x)-3]))
  corr_fifth <- apply(corr, 1, function(x) which(x == sort(x,partial=length(x)-4)[length(x)-4]))
  
  # Low confidence subtypes
  corr_low_conf <-apply(corr, 1, function(x) which(x[which.max(x)] < und_cutoff))
  n_corr_low_conf <- which(corr_low_conf == TRUE)
  names_corr_low_conf <- colnames(corr)[corr_max_w]
  names_corr_low_conf[n_corr_low_conf] <- "Undetermined"
  
  # Mixed subtypes
  corr_mixed <-apply(corr, 1, function(x) which((x[which.max(x)] -  x[which(x == sort(x,partial=length(x)-1)[length(x)-1])]) < mix_cutoff))
  n_mixed <- which(corr_mixed == TRUE)
  names_mixed <- colnames(corr)[corr_max_w]
  names_mixed[n_mixed] <- "mixed"
  
  # making a matrix of subtype identity for samples with PAM centroids
  sub <- cbind(rownames(corr), corr, colnames(corr)[corr_max_w], names_corr_low_conf, names_mixed, colnames(corr)[corr_second], colnames(corr)[corr_third], colnames(corr)[corr_fourth], colnames(corr)[corr_fifth]) 
  colnames(sub)[1] <- "Samples"
  colnames(sub)[7] <- "First subtype"
  colnames(sub)[8] <- "Undetermined subtypes"
  colnames(sub)[9] <- "Mixed subtypes"
  colnames(sub)[10] <- "Second subtype"
  colnames(sub)[11] <- "Third subtype"
  colnames(sub)[12] <- "Fourth subtype"
  colnames(sub)[13] <- "Fifth subtype"
  
  # calculating table of subtypes for statistics
  sub_table <- table(colnames(corr)[corr_max_w])
  
  #writing data
  write.table(sub, outfile, row.names = FALSE, sep = "\t", quote = FALSE)
  
  #writing subtype table
  write.table(sub_table, outfile_tab, row.names = FALSE, sep = "\t", quote = FALSE)
  
  #writing CRCAssigner-786 data
  write.table(data_m_med[m[w],], outfile_file, sep = "\t", quote = FALSE)
  
  
  #writing correlation
  write.table(corr, outfile_cor, sep = "\t", quote = FALSE)
  
  return( sub_table )
  
}
