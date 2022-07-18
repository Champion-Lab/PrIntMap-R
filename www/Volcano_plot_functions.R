read_peptide_csv_PEAKS_volcano <- function(peptide_file, sample_pattern, min_valid_sample = 2,
                                           intensity_metric = "PSM") {
  check_file(peptide_file, "PEAKS")
  peptide_import <- read.csv(peptide_file)
  filetype(peptide_import, "Combined", "PEAKS")
  names(peptide_import)[names(peptide_import) == "X.Spec"] <- "total_spectra"
  peptides <- peptide_import
  peptides$sequence <- str_remove_all(peptides$Peptide, "[a-z1-9()+-:.]")
  peptides$PEPTIDE <- peptides$Peptide
  peptides$protein <- NA
  for(i in 1:nrow(peptides)) {
    protein <- str_split(peptides$Accession[i], ";")[[1]]
    protein_vec <- rep("", length(protein))
    for (j in 1:length(protein)) {
      protein_vec[j] <- str_split(protein[j], "\\|")[[1]][1]
    }
    peptides$protein[i] <- paste(protein_vec, sep = ";", collapse = ";")
  }
  
  if(length(names(peptides)[grepl(sample_pattern, names(peptides))])<=0){
    stop("Sample Pattern not found in file.")
  }
  else {
    if (intensity_metric == "PSM") {
      pattern <- paste0("X\\.Spec", ".*", sample_pattern, ".*")
    } else if (intensity_metric == "Intensity") {
      pattern <- paste0("Intensity", ".*", sample_pattern, ".*")
    } else if (intensity_metric == "Area") {
      pattern <- paste0("Area", ".*", sample_pattern, ".*")
    }
    dataframe <- peptide_import[,grepl(pattern, names(peptide_import))]
    if (intensity_metric == "PSM") {
      dataframe[dataframe == 0] <- NA
    }
    sample_count <- ncol(as.data.frame(dataframe))
    for (i in 1:sample_count) {
      peptides[[paste0("Volcano_intensity_", i)]] <- NA
    }
    for (i in 1:nrow(peptides)) {
      for (j in 1:sample_count) {
        peptides[[paste0("Volcano_intensity_", j)]][i] <- dataframe[i,j]
      }
    }
    return(peptides)
  }}

combine_two_volcano_dfs <- function(df_1, df_2, min_valid_sample = 2, fdr = 0.05,
                                    fold_change_cutoff_plot = 1, fold_change_cutoff_sig = 5,
                                    equal_variance_bool = T, remove_na = T, set_na = 0) {
  combine_df <- full_join(df_1, df_2, by = c("PEPTIDE", "sequence", "protein"))
  
  volcano_df1 <- combine_df[,grepl("Volcano_intensity_.*\\.x", names(combine_df))]
  volcano_df2 <- combine_df[,grepl("Volcano_intensity_.*\\.y", names(combine_df))]
  
  combine_df$count.x <- NA
  combine_df$count.y <- NA
  combine_df$valid.x <- F
  combine_df$valid.y <- F
  combine_df$fold_change_xy <- NA
  combine_df$fold_change_category <- "Normal"
  combine_df$p_val <- NA
  combine_df$avg.x <- NA
  combine_df$avg.y <- NA
  combine_df$l2fc_xy <- NA
  combine_df$neg10logp <- NA
  combine_df$p_val_error <- F
  
  for (i in 1:nrow(combine_df)) {
    combine_df$count.x[i] <- sum(!is.na(as.numeric(volcano_df1[i,])))
    combine_df$count.y[i] <- sum(!is.na(as.numeric(volcano_df2[i,])))
  }
  
 
  combine_df$valid.x[combine_df$count.x >= min_valid_sample] <- T
  combine_df$valid.y[combine_df$count.y >= min_valid_sample] <- T
  combine_df$fold_change_category[combine_df$valid.x == F & combine_df$valid.y == F] <- "Invalid"
  combine_df$fold_change_category[combine_df$valid.x == T & combine_df$count.y == 0] <- "Negative_Infinite"
  combine_df$fold_change_category[combine_df$count.x == 0 & combine_df$valid.y == T] <- "Infinite"
  combine_df$fold_change_category[(combine_df$valid.x == T & combine_df$valid.y == F & combine_df$count.y > 0)] <- "Compromised.x"
  combine_df$fold_change_category[(combine_df$valid.y == T & combine_df$valid.x == F & combine_df$count.x > 0)] <- "Compromised.y"
  
  for (i in 1:nrow(combine_df)) {
    if (combine_df$fold_change_category[i] == "Normal" || combine_df$fold_change_category[i] == "Compromixed.x" || combine_df$fold_change_category[i] == "Compromised.y") {
      if (remove_na == T) {
        combine_df$avg.x[i] <- mean(as.numeric(volcano_df1[i,]), na.rm = T)
        combine_df$avg.y[i] <- mean(as.numeric(volcano_df2[i,]), na.rm = T)
        
        tryCatch({
          combine_df$p_val[i] <- t.test(as.numeric(volcano_df1[i,]),
                                        as.numeric(volcano_df2[i,]),
                                        var.equal = equal_variance_bool)[[3]]
        }, error = function(e){
          combine_df$p_val[i] <- NA
          combine_df$p_val_error[i] <- T})
      } else {
        volcano_df1[is.na(volcano_df1)] <- set_na
        volcano_df2[is.na(volcano_df2)] <- set_na
        combine_df$avg.x[i] <- mean(as.numeric(volcano_df1[i,]), na.rm = F)
        combine_df$avg.y[i] <- mean(as.numeric(volcano_df2[i,]), na.rm = F)
        combine_df$p_val[i] <- t.test(as.numeric(volcano_df1[i,]),
                                      as.numeric(volcano_df2[i,]),
                                      var.equal = equal_variance_bool)[[3]]
        
      }
      combine_df$fold_change_xy[i] <- combine_df$avg.y[i] / combine_df$avg.x[i]
      combine_df$l2fc_xy[i] <- log(combine_df$fold_change_xy[i], base = 2)
      combine_df$neg10logp[i] <- -log(combine_df$p_val[i], base = 10)
    }
  }
  
  
  return(combine_df)
}

create_volcano_plot <- function(df, fdr = 0.05,
                                fold_change_cutoff_plot = 1, fold_change_cutoff_sig = 5,
                                equal_variance_bool = T, intensity_metric = "PSM",
                                sample1 = "Sample_1", sample2 = "Sample_2", BH_correction = F,
                                protein_of_interest = NULL, display_comp_vals = T, display_infinites = T) {
  
  minx <- min(df$l2fc_xy, na.rm = T)
  maxx <- max(df$l2fc_xy, na.rm = T)
  maxy <- max(df$neg10logp, na.rm = T)
  rangex <- maxx - minx
  
  if (display_comp_vals == F) {
    df <- df[df$fold_change_category != "Compromised.x",]
    df <- df[df$fold_change_category != "Compromised.y",]
  }
  
  if (BH_correction == T) {
    df <- df[order(df$p_val),]
    BH_df <- df[!is.na(df$p_val),]
    df$critical_val_BH <- NA
    for (i in 1:nrow(BH_df)) {
      df$critical_val_BH[i] <- (i / nrow(BH_df)) * fdr
    }
    df$BH_cutoff <- df$critical - df$p_val
    sigdf <- df[df$BH_cutoff > 0,]
    sigdf <- sigdf[!is.na(sigdf$BH_cutoff),]
    print(nrow(sigdf))
    if (nrow(sigdf) == 0) {
      y_cutoff <- maxy * 1.1
    } else {
      new_p_val <- max(sigdf$p_val, na.rm = T)
      y_cutoff <- -log(new_p_val, base = 10)
    }
  } else {
    y_cutoff <- -log(fdr, base = 10)
  }
  
 
  df$color <- "Not_Significant"
  df$color[df$l2fc_xy > fold_change_cutoff_plot & df$neg10logp > y_cutoff] <- "Significant"
  df$color[df$l2fc_xy < -fold_change_cutoff_plot & df$neg10logp > y_cutoff] <- "Significant"
  
  if (!is.null(protein_of_interest)) {
    df$color[grepl(protein_of_interest, df$protein)] <- "Prot_of_interest"
  }
  
  if (length(unique(df$color)) == 1) {
    color_list <- "grey65"
  } else if (length(unique(df$color)) == 2) {
    color_list <- c("grey65", "black")
  } else {
    color_list <- c("grey65", "green", "black")
  }
  
  
  minx <- min(df$l2fc_xy, na.rm = T)
  maxx <- max(df$l2fc_xy, na.rm = T)
  maxy <- max(df$neg10logp, na.rm = T)
  rangex <- maxx - minx
  
  pos_range <- maxx - fold_change_cutoff_plot
  neg_range <- minx + fold_change_cutoff_plot
  y_range <- maxy - y_cutoff
  
  infinite_x <- maxx - (0.1*maxx)
  infinite_y <- maxy - (0.2*maxy)
  neg_infinite_x <- minx + (0.1*abs(minx))
  
  
  infinite <- df[df$fold_change_category == "Infinite",]
  neg_infinite <- df[df$fold_change_category == "Negative_Infinite",]
  
  infinite$x <- runif(nrow(infinite), (fold_change_cutoff_plot + (pos_range/2)*1.2), maxx)
  infinite$y <- runif(nrow(infinite), (y_cutoff + (y_range/2)*1.4), maxy)
  neg_infinite$x <- -1 * (runif(nrow(neg_infinite), (fold_change_cutoff_plot + (abs(neg_range)/2)*1.2), abs(minx)))
  neg_infinite$y <- runif(nrow(neg_infinite), (y_cutoff + (y_range/2)*1.4), maxy)
  
  if (is.null(protein_of_interest)) {
    plot_title <- paste0("Peptide Volcano Plot: ", intensity_metric)
  } else {
    plot_title <- paste0("Peptide Volcano Plot: ", intensity_metric, " (", protein_of_interest, ")")
  }
  
  plot <- ggplot() +
    geom_vline(xintercept = fold_change_cutoff_plot, linetype = 2) +
    geom_vline(xintercept = -fold_change_cutoff_plot, linetype = 2) +
    geom_hline(yintercept = y_cutoff, linetype = 2) +
    geom_point(data = df, aes(x = l2fc_xy, y = neg10logp,
                              PEPTIDE = PEPTIDE,
                              sequence = sequence,
                              protein = protein,
                              colour = color,
                              p_val = p_val,
                              fold_change_category = fold_change_category)) +
    
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(x = (paste("Log2 fold-change: (", sample2, "/", sample1, ")",sep="")),
         y = "-Log10 (p-value)",
         title = plot_title) +
    scale_x_continuous(breaks = round(min(df$l2fc_xy, na.rm = T)):round(max(df$l2fc_xy, na.rm = T))) +
    scale_y_continuous(breaks = 0:round(max(df$neg10logp, na.rm = T)))+
    scale_color_manual(values = color_list)
  
  if (display_infinites) {
    plot <- plot + geom_point(data = infinite, 
                               aes(x = x, y = y,
                                   PEPTIDE = PEPTIDE,
                                   sequence = sequence,
                                   protein = protein,
                                   colour = color,
                                   fold_change_category = fold_change_category)) +
      geom_point(data = neg_infinite, 
                  aes(x = x, y = y,
                      PEPTIDE = PEPTIDE,
                      sequence = sequence,
                      protein = protein,
                      colour = color,
                      fold_change_category = fold_change_category))
  }
  
  return(plot)
}

create_volcano_plotly <- function(plot) {
  plotly <- ggplotly(plot)
  return(plotly)
}

