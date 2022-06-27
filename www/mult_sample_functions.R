plot_difference_comb_mult <- function(AA_df, protein, intensity_label = "PSM",
                                 font_size = 15,
                                 alpha = 1, plot_type = "line") {
  plot <- ggplot(data = AA_df) +
    geom_hline(yintercept = 0, color = "red", size = 0.1, linetype = "dashed") +
    geom_line(aes(x = AA_index, y = difference, color = sample), alpha = alpha) +
    geom_point(aes(x = AA_index, y = difference, AA = AA, sample = sample),
               size = 0, alpha = 0) +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(), legend.position = "right") +
    labs(title = protein, x = "Amino Acid Number", y = paste0("Difference ", intensity_label, " (Sx - S1)"))+
    scale_x_continuous(breaks = pretty_breaks(10))
  
  return(plot)
}

plot_foldchange_comb_mult <- function(AA_df, protein, intensity_label = "PSM",
                                 font_size = 15,
                                 alpha = 1, plot_type = "line") {
  plot <- ggplot(data = AA_df) +
    geom_hline(yintercept = 1, color = "red", size = 0.1, linetype = "dashed") +
    geom_line(aes(x = AA_index, y = fold_change, color = sample), alpha = alpha) +
    geom_point(aes(x = AA_index, y = fold_change, AA = AA, sample = sample),
               size = 0, alpha = 0) +
    geom_point(data = AA_df[AA_df$fold_change_label == "Zero",], aes(x = AA_index,
                                                                     y = fold_change,
                                                                     AA = AA,
                                                                     fold_change_label = fold_change_label,
                                                                     sample = sample),
               size = 0.5, color = "red") +
    geom_point(data = AA_df[AA_df$fold_change_label == "Infinite",], aes(x = AA_index,
                                                                         y = fold_change,
                                                                         AA = AA,
                                                                         fold_change_label = fold_change_label,
                                                                         sample = sample),
               size = 0.5, color = "green") +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(), legend.position = "right",
          legend.text = element_text(size = 8)) +
    labs(title = protein, x = "Amino Acid Number", y = paste0("Ratio ", intensity_label, " (Sx / S1)")) +
    scale_alpha_manual(values = c(0.5, 0.75, 1))+
    scale_x_continuous(breaks = pretty_breaks(10))
  
  return(plot)
}