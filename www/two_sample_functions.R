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
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)
  } else {
    plot <- ggplot(data = AA_df) +
      geom_line(aes(x = AA_index, y = intensity, color = sample), alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, AA = AA, origin_pep = origin_pep, color = sample),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)
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
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)
  } else {
    plot <- ggplot(data = AA_df) +
      geom_line(aes(x = AA_index, y = intensity, color = sample), alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, AA = AA, sample = sample),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)
  }
  return(plot)
}