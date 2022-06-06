#create dataframe for making stacked plot
create_stack_line_df <- function(peptide_df,
                                 protein,
                                 intensity = "PSM") {
  vector_list <- list()
  for (i in 1:nrow(peptide_df)) {
    peptide <- peptide_df$sequence[i]
    intensity_value <- peptide_df[[intensity]][i]
    if (is.na(intensity_value)) {
      intensity_value <- 0
    }
    matches_df <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    if (matches_count > 0) {
      for (j in 1:matches_count) {
        start <- matches_df[[j,1]]
        end <- matches_df[[j,2]]
        length <- end - start + 1
        if (matches_count > 1) {
          repeated <- T
        } else {
          repeated <- F
        }
        peptide_vector <- c(peptide, as.numeric(intensity_value),
                            as.numeric(length), as.numeric(start),
                            as.numeric(end), as.logical(repeated))
        vector_list[[length(vector_list)+1]] <- peptide_vector
      }
    }
  }
  peptide_df <- as.data.frame(do.call(rbind, vector_list))
  names(peptide_df) <- c("peptide", "intensity_value", "length", "start", "end", "repeated")
  peptide_df$intensity_value <- as.numeric(peptide_df$intensity_value)
  peptide_df$length <- as.numeric(peptide_df$length)
  peptide_df$start <- as.numeric(peptide_df$start)
  peptide_df$end <- as.numeric(peptide_df$end)
  peptide_df$repeated <- as.logical(peptide_df$repeated)
  peptide_df <- peptide_df %>%
    group_by(peptide, length, start, end, repeated) %>%
    dplyr::summarise(intensity_value = sum(intensity_value)) %>%
    as.data.frame()
  peptide_df <- peptide_df[order(peptide_df$start, peptide_df$length),]
  peptide_df$y_val <- 1:nrow(peptide_df)
  peptide_df$intensity_value[peptide_df$intensity_value == 0] <- NA
  return(peptide_df)
}

#create plot for stacked plot
create_stacked_line_plot <- function(peptide_dataframe, protein_name, protein_seq){
  plot <- ggplot(peptide_dataframe) +
    geom_segment(aes(y = y_val, x = start, xend = end, yend = y_val,
                     color = intensity_value, length = length, peptide = peptide),
                 size = 1, lineend = "round") +
    theme_bw(base_size = 15) +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Amino Acid Position", y = element_blank(),
         title = protein_name, color = element_blank()) +
    scale_x_continuous(limits = c(0, nchar(protein_seq)), breaks = pretty_breaks(10)) +
    scale_colour_gradient(low = "turquoise", high = "red", na.value = "black")
  return(plot)
}

#create plotly for stacked plot
create_stack_line_plotly <- function(plot) {
  plotly_obj <- ggplotly(plot, tooltip = c("length", "peptide",
                                           "x", "colour"))
  return(plotly_obj)
}