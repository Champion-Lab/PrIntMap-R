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
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
      scale_x_continuous(breaks = pretty_breaks(10))
  }
  return(plot)
}

plot_intensity_comb <- function(AA_df, protein, intensity_label = "PSM",
                                font_size = 15,
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
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
      scale_x_continuous(breaks = pretty_breaks(10))
  }
  return(plot)
}


plot_difference_comb <- function(AA_df, protein, intensity_label = "PSM",
                                 font_size = 15,
                                 alpha = 1, plot_type = "line") {
  plot <- ggplot(data = AA_df) +
    geom_line(aes(x = AA_index, y = difference), alpha = alpha) +
    geom_point(aes(x = AA_index, y = difference, AA = AA),
               size = 0, alpha = 0) +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(title = protein, x = "Amino Acid Number", y = "Difference in Intensity (S2 - S1)")+
    scale_x_continuous(breaks = pretty_breaks(10))
  
  return(plot)
}

plot_foldchange_comb <- function(AA_df, protein, intensity_label = "PSM",
                                 font_size = 15,
                                 alpha = 1, plot_type = "line") {
  plot <- ggplot(data = AA_df) +
    geom_line(aes(x = AA_index, y = fold_change), alpha = alpha) +
    geom_point(aes(x = AA_index, y = fold_change, AA = AA),
               size = 0, alpha = 0) +
    geom_point(data = AA_df[AA_df$fold_change_label == "Zero",], aes(x = AA_index,
                                                                     y = fold_change,
                                                                     AA = AA),
               size = 0.5, color = "red") +
    geom_point(data = AA_df[AA_df$fold_change_label == "Infinite",], aes(x = AA_index,
                                                                     y = fold_change,
                                                                     AA = AA),
               size = 0.5, color = "green") +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(title = protein, x = "Amino Acid Number", y = "Ratio Intensity (S2 / S1)") +
    scale_alpha_manual(values = c(0.5, 0.75, 1))+
    scale_x_continuous(breaks = pretty_breaks(10))
  
  return(plot)
}