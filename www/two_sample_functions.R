plot_origin_comb <- function(AA_df, protein,
                             intensity_label = "PSM",
                             font_size = 15, 
                             alpha = 1, plot_type = "line") {
  if (plot_type == "bar") {
    plot <- ggplot(data = AA_df) +
      geom_bar(aes(x = AA_index, y = intensity, fill = sample),
               stat = "identity", position = "dodge", alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, AA = AA, origin_pep = origin_pep, color = sample),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label) +
      scale_x_continuous(breaks = pretty_breaks(10))
  } else {
    plot <- ggplot(data = AA_df) +
      geom_line(aes(x = AA_index, y = intensity, color = sample), alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, AA = AA, origin_pep = origin_pep, color = sample),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank(), legend.position = "right",
            legend.text = element_text(size = 8)) +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
      scale_x_continuous(breaks = pretty_breaks(10))
  }
  return(plot)
}

plot_intensity_comb <- function(AA_df, protein, intensity_label = "PSM",
                                font_size = 18,
                                alpha = 1, plot_type = "line") {
  if (plot_type == "bar") {
    plot <- ggplot(data = AA_df) +
      geom_bar(aes(x = AA_index, y = intensity, color = sample),
               stat = "identity", position = "dodge",
               fill = color, alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, AA = AA, sample = sample),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank()) +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
      scale_x_continuous(breaks = pretty_breaks(10))
  } else {
    plot <- ggplot(data = AA_df) +
      geom_line(aes(x = AA_index, y = intensity, color = sample), alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, AA = AA, sample = sample),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank(), legend.position = "right",
            legend.text = element_text(size = 8)) +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
      scale_x_continuous(breaks = pretty_breaks(10))
  }
  return(plot)
}


plot_difference_comb <- function(AA_df, protein, intensity_label = "PSM",
                                 font_size = 18,
                                 alpha = 1, plot_type = "line") {
  plot <- ggplot(data = AA_df) +
    geom_hline(yintercept = 0, color = "red", size = 0.1, linetype = "dashed") +
    geom_line(aes(x = AA_index, y = difference), alpha = alpha) +
    geom_point(aes(x = AA_index, y = difference, AA = AA),
               size = 0, alpha = 0) +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(title = protein, x = "Amino Acid Number", y = paste0("Difference ", intensity_label, " (S2 - S1)"))+
    scale_x_continuous(breaks = pretty_breaks(10))
  
  return(plot)
}

plot_foldchange_comb <- function(AA_df, protein, intensity_label = "PSM",
                                 font_size = 18,
                                 alpha = 1, plot_type = "line") {
  plot <- ggplot(data = AA_df) +
    geom_hline(yintercept = 1, color = "red", size = 0.1, linetype = "dashed") +
    geom_line(aes(x = AA_index, y = fold_change), alpha = alpha) +
    geom_point(aes(x = AA_index, y = fold_change, AA = AA),
               size = 0, alpha = 0) +
    geom_point(data = AA_df[AA_df$fold_change_label == "Zero",], aes(x = AA_index,
                                                                     y = fold_change,
                                                                     AA = AA,
                                                                     fold_change_label = fold_change_label),
               size = 0.5, color = "red") +
    geom_point(data = AA_df[AA_df$fold_change_label == "Infinite",], aes(x = AA_index,
                                                                     y = fold_change,
                                                                     AA = AA,
                                                                     fold_change_label = fold_change_label),
               size = 0.5, color = "green") +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(), legend.position = "right",
          legend.text = element_text(size = 8)) +
    labs(title = protein, x = "Amino Acid Number", y = paste0("Ratio ", intensity_label, " (S2 / S1)")) +
    scale_alpha_manual(values = c(0.5, 0.75, 1))+
    scale_x_continuous(breaks = pretty_breaks(10))
  
  return(plot)
}

#create dataframe for making a 2 sample stacked plot
create_stack_line_df2 <- function(peptide_df1, peptide_df2,
                                  sample1, sample2,
                                  protein, 
                                  intensity1="PSM", intensity2="PSM") {
  stack_line_dataframe1 <- create_stack_line_df(peptide_df1, protein, intensity1)
  stack_line_dataframe1$sample <- sample1
  stack_line_dataframe2 <- create_stack_line_df(peptide_df2, protein, intensity2)
  stack_line_dataframe2$sample <- sample2
  combo_stack_line_dataframe <- bind_rows(stack_line_dataframe1, stack_line_dataframe2)
  combo_stack_line_dataframe <- combo_stack_line_dataframe[order(combo_stack_line_dataframe$start, combo_stack_line_dataframe$sample), ]
  combo_stack_line_dataframe$y_val <- 1:nrow(combo_stack_line_dataframe)
  return(combo_stack_line_dataframe)
}

#create plot for 2 sample stacked plot with y value as y_val
create_stacked_line_plot_yval2 <- function(peptide_dataframe, protein_name, protein_seq){
  plot <- ggplot(data = peptide_dataframe) +
    geom_segment(aes(y = y_val, x = start, xend = end, yend = y_val,
                     length = length, peptide = peptide, color=sample),
                 size = 1, lineend = "round") +
    theme_bw(base_size = 18) +
    theme(panel.grid = element_blank(),
          legend.position = "right", legend.text = element_text(size = 8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Amino Acid Position", y = element_blank(),
         title = protein_name)
  return(plot)
}

#create plot for 2 sample stacked plot with intensity as y_val
create_stacked_line_plot_intensity2 <- function(peptide_dataframe, protein_name, protein_seq,
                                                intensity_label = "Intensity"){
  plot <- ggplot(data = peptide_dataframe) +
    geom_segment(aes(y = intensity_value, x = start, xend = end, yend = intensity_value,
                     length = length, peptide = peptide, color=sample),
                 size = 1, lineend = "round") +
    theme_bw(base_size = 18) +
    theme(panel.grid = element_blank(),
          legend.position = "right", legend.text = element_text(size = 8)) +
    labs(x = "Amino Acid Position", y = intensity_label,
         title = protein_name)
  return(plot)
}