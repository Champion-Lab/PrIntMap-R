create_annotation_vec <- function(pattern, protein) {
  split_AA <- str_split(protein, "")
  AA_df <- data.frame(AA = split_AA[[1]])
  AA_df$AA_index <- 1:nrow(AA_df)
  motifs <- str_locate_all(protein, pattern)[[1]]
  matches_count <- nrow(motifs)
  motif_vector <- rep(F, nrow(AA_df))
  if (matches_count > 0) {
    for (j in 1:matches_count) {
      start <- motifs[[j,1]]
      end <- motifs[[j,2]]
      motif_vector[start:end] <- T
    }
  }
  return(motif_vector) 
}

#create annotation dataframe
create_annotation_df <- function(AA_df, annotation_vec,
                                 annotation_name = "annotation"){
  AA_df[[annotation_name]] <- annotation_vec
  annotation_df <- AA_df[AA_df[[annotation_name]] == T,]
  return(annotation_df)
}

#take a plot and add verticle lines at the annotation location
add_annotation_layer <- function(plot,
                                 annotation_df,
                                 color = "yellow", annot_title = "annotation") {
  plot <- plot +
    geom_vline(xintercept = annotation_df$AA_index, color = color,
               size = 0.5)
  return(plot)
}