#create vector with repeat information
create_repeat_vec <- function(peptide_df,
                              protein,
                              protein_name,
                              database,
                              intensity = "PSM") {
  vector_list_proteome_repeat <- list()
  vector_list_intraprotein_repeat <- list()
  for (i in 1:nrow(peptide_df)) {
    peptide <- peptide_df$sequence[i]
    matches_df <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    proteome_repeat_vec <- rep(0, nchar(protein))
    intraprotein_repeat_vec <- rep(0, nchar(protein))
    if (matches_count > 0) {
      for (j in 1:matches_count) {
        start <- matches_df[[j,1]]
        end <- matches_df[[j,2]]
        if (matches_count > 1) {
          intraprotein_repeat_vec[start:end] <- 1
        }
      }
      for (k in 1:length(database)) {
        protein_matches <- str_count(database[[k]], peptide)
        if (protein_matches > 0) {
          name <- attr(database[[k]], "name")
          if (name != protein_name){
            for (j in 1:matches_count) {
              start <- matches_df[[j,1]]
              end <- matches_df[[j,2]]
              proteome_repeat_vec[start:end] <- 1
            }
          }
        }
      }
    }
    vector_list_proteome_repeat[[i]] <- proteome_repeat_vec
    vector_list_intraprotein_repeat[[i]] <- intraprotein_repeat_vec
  }
  proteome_repeat_df <- as.data.frame(sapply(vector_list_proteome_repeat, unlist))
  proteome_repeat_sum <- rowSums(proteome_repeat_df)
  intraprotein_repeat_df <- as.data.frame(sapply(vector_list_intraprotein_repeat, unlist))
  intraprotein_repeat_sum <- rowSums(intraprotein_repeat_df)
  
  repeat_vec <- rep("unique", nchar(protein))
  for (i in 1:nchar(protein)) {
    if ((intraprotein_repeat_sum[i] > 0) && (proteome_repeat_sum[i] == 0)) {
      repeat_vec[i] <- "unique_repeated"
    }
    if ((intraprotein_repeat_sum[i] == 0) && (proteome_repeat_sum[i] > 0)) {
      repeat_vec[i] <- "non_unique"
    }
    if ((intraprotein_repeat_sum[i] > 0) && (proteome_repeat_sum[i] > 0)) {
      repeat_vec[i] <- "non_unique_repeated"
    }
  }
  return(repeat_vec)
}


create_unique_plot_origin <- function(AA_df, protein,
                               intensity_label = "PSM",
                               font_size = 15) {
  plot <- ggplot(data = AA_df) +
    geom_line(aes(x = AA_index, y = intensity)) +
    geom_bar(aes(x = AA_index, y = intensity, fill = repeated, AA = AA),
             position = "dodge", stat = 'identity')+
    geom_point(aes(x = AA_index, y = intensity, AA = AA, fill = repeated,
                   origin_pep = origin_pep),
               size = 0, alpha = 0) +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(), legend.position = "right",
          legend.text = element_text(size = 8)) +
    labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
    scale_x_continuous(breaks = pretty_breaks(10))
  
  return(plot)
}

create_unique_plot_intensity <- function(AA_df, protein,
                                      intensity_label = "PSM",
                                      font_size = 15) {
  plot <- ggplot(data = AA_df) +
    geom_line(aes(x = AA_index, y = intensity)) +
    geom_bar(aes(x = AA_index, y = intensity, fill = repeated, AA = AA),
             position = "dodge", stat = 'identity')+
    geom_point(aes(x = AA_index, y = intensity, AA = AA, fill = repeated),
               size = 0, alpha = 0) +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(), legend.position = "right",
          legend.text = element_text(size = 8)) +
    labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
    scale_x_continuous(breaks = pretty_breaks(10))
  
  return(plot)
}