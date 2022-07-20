
#create new dataframe with only peptides of interest

create_PTM_df_PEAKS <- function(peptide_df, protein){
  modpeptide_store <- vector()
  peptide_store<- vector()
  for(i in 1:nrow(peptide_df)){
    peptide <- peptide_df$sequence[i]
    matches_df <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    if(matches_count >0){
      peptide_store <- append(peptide_store, peptide)
      modpeptide_store <- append(modpeptide_store, peptide_df$Peptide[i])

    }
  }
  peptide_df <- data.frame(peptide_store, modpeptide_store)
  print(peptide_df)
  return(peptide_df)
}

##function to pull out modifications and add to new col in dataframe
add_mod_PEAKS <- function(peptide_df, regex_pattern){
  PTM_store <- list()
  for(i in 1:nrow(peptide_df)){
    mod_peptide <- peptide_df$modpeptide_store[i]
    peptide <- peptide_df$peptide_store[i]
    PTM_value <- rep("", nchar(peptide))
    PTM_df <-str_locate_all(mod_peptide, regex_pattern)[[1]]
    print(PTM_df)
    PTM_count <- nrow(PTM_df)
    if(PTM_count > 0){
      for(j in 1:PTM_count){
        start <- PTM_df[[j,1]]
        PTM_value[[start-1]] <- regex_pattern
      }
    }
  
    PTM_store[[i]] <- PTM_value
  }
 
  peptide_df[[regex_pattern]] <- PTM_store
  return(peptide_df)
  }


##create new df the length of #AA in protein, including PTMs, intensity, etc.
create_PTM_plot_df_PEAKS <- function(peptide_df, protein, regex_pattern, intensity_vector){
  vector_list <- list()
  split_AA <- str_split(protein, "")
  AA_df <- data.frame(AA = split_AA[[1]])
  AA_df$AA_index <- 1:nrow(AA_df)
  for(i in 1:nrow(peptide_df)){
    PTM_vector <- rep("", nchar(protein))
    PTM <- peptide_df[[regex_pattern]][i]
    print(PTM)
    peptide <- peptide_df$peptide_store[i]
    matches_df <-str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    if(matches_count >0){
    for(j in 1:matches_count){
      start <- matches_df[[j,1]]
      print(start)
      end <- matches_df[[j,2]]
      print(end)
      PTM_vector[start:end] <- PTM[[1]]
  
    }
    }
    vector_list[[i]] <- PTM_vector
    }
  PTM_df <- as.data.frame(sapply(vector_list, unlist))
  PTM_vec <- apply(PTM_df, 1, paste, collapse="")
  AA_df[[regex_pattern]] <- PTM_vec
  AA_df$Intensity <- intensity_vector
  PTM_df2 <- subset(AA_df, AA_df[[regex_pattern]] != "")
  return(PTM_df2)
}

  
  
#add PTM layer to plots 
##add a new layer for eachPTM
add_PTM_layer <- function(plot, PTM_df, color = "red" ){
  print(PTM_df)
  plot <- plot +
    geom_point(data = PTM_df, x= PTM_df$AA_index, y = PTM_df$Intensity,  color = color ,
           size = .5) 

  return(plot)
}