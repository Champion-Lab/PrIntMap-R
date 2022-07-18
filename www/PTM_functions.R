

create_PTM_vec_PEAKS <- function(peptide_df, protein,){
  for(i in 1:nrow(peptide_df)){
    mod_peptide<-peptide_df$Peptide[i]
    peptide <- peptide_df$sequence[i]
    matches_df <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    PTM_vector <- rep("", nchar(protein))
    PTM_df <-str_locate_all(mod_peptide, "\\(.?.?.?.?.?.?\\)")[[1]]
    PTM_count <- nrow(PTM_df)
    vector_list <- list()
    AA_vec <- vector()
    PTM_value <- vector()
    if (matches_count > 0){
      mod_peptide_strsplit <- strsplit(mod_peptide, "")[[1]]
      for(j in 1: length(mod_peptide_strsplit)){
        if (!grepl("[a-z1-9()+-:.]", mod_peptide_strsplit[j])){
          AA_vec <- c(AA_vec, mod_peptide_strsplit[j])
          PTM_value <- c(PTM_value, "")
          for (k in 1:PTM_count) {
            start <- PTM_df[[k,1]]
            end <- PTM_df[[k,2]]
            if(j == start-1){
            PTM_value<- PTM_value[-length(PTM_value)]
            PTM_value<- c(PTM_value, substr(mod_peptide, start, end))
            }
          }}}
        for (j in 1:matches_count) {
          start <- matches_df[[j,1]]
          end <- matches_df[[j,2]]
          PTM_vector[start:end] <- PTM_value
        }
    }
    vector_list[[i]] <- PTM_vector
  }

  PTM_df <- as.data.frame(sapply(vector_list, unlist))
  PTM_vec <- apply(PTM_df, 1, paste, collapse="")
  return(PTM_vec)
  }
          
#remove unwanted PTMs
refine_PTM_vec <- function(AA_df, PTM_pattern_vec, PTM_names_vec){
  for(i in 1:length(PTM_pattern_vec)){
    for(j in 1:nrows(AA_df)){
      store_vec <- vector()
      count_vec <- vector()
      PTM <- AA_df$PTMs[i]
      matches_df <- str_locate_all(PTM, PTM_pattern_vec[i])[[1]]
      matches_count <- nrows(matches_df)
      count_vec <- c(count_vec, matches_count)
      if(matches_count >0){
        store_vec <- c(store_vec, PTM_pattern_vec[[i]])
      }else{
        store_vec <- c(store_vec, "")
      }
    }
    AA_df[[paste0(PTM_names_vec[[i]], " Count")]] <- count_vec
    AA_df$PTM_names_vec[[i]] <- store_vec
  }
  return(AA_df)
}

#combine PTM vector with AA_df
create_PTM_df <- function(AA_df,
                           PTM_vec) {
  AA_df$PTMs <- PTM_vec
  return(AA_df)
}

#add PTM layer to plots
add_PTM_layer <- function(plot,
                          PTM_df, PTM_value,
                          color = "yellow") {
  
  plot <- plot +
    
    geom_vline(xintercept = annotation_df$AA_index, color = color,
               size = 0.5)
  return(plot)
}