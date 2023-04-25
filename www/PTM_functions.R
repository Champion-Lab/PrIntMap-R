
#create new dataframe with only peptides of interest

create_PTM_df_PEAKS <- function(peptide_df, protein, regex_pattern, intensity = "PSM"){
  modpeptide_store <- vector()
  peptide_store<- vector()
  intensity_store <- vector()
  for(i in 1:nrow(peptide_df)){
    peptide <- peptide_df$sequence[i]
    matches_df <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    if(matches_count >0){
      peptide_store <- append(peptide_store, peptide)
      modpeptide_store <- append(modpeptide_store, peptide_df$Peptide[i])
      intensity_store <- append(intensity_store, peptide_df[[intensity]][i])

    }
  }
  peptide_df <- data.frame(peptide_store, modpeptide_store, intensity_store)
  intensity_vector <- vector()
  PTM_store <- list()
  for(i in 1:nrow(peptide_df)){
    mod_peptide <- peptide_df$modpeptide_store[i]
    mod_peptide_strsplit <- strsplit(mod_peptide, "")[[1]]
    peptide <- peptide_df$peptide_store[i]
    PTM_df <-str_locate_all(mod_peptide, fixed(regex_pattern))[[1]]
    PTM_count <- nrow(PTM_df)
    if(PTM_count > 0){
      if(is.na(peptide_df$intensity_store[[i]])){
        intensity_value <- 0
      }else{      
        intensity_value <- peptide_df$intensity_store[[i]]
}
      for(j in 1:PTM_count){
        start <- PTM_df[[j,1]]
        end <- PTM_df[[j,2]]
        mod_peptide_strsplit[start-1] <- "x"
      }
    }else{
      intensity_value <- 0
    }
    intensity_vector<- c(intensity_vector, intensity_value)
    length <- length(mod_peptide_strsplit)
    for(k in length:1){
      if (!grepl("[a-zA-Z]", mod_peptide_strsplit[k])){
        mod_peptide_strsplit<- mod_peptide_strsplit[-k]
      } else if(!grepl("x", mod_peptide_strsplit[k])){
        mod_peptide_strsplit[[k]] <- ""
      }
      }
    PTM_store[[i]] <- mod_peptide_strsplit
  }
  peptide_df$intensity_store_mod <- intensity_vector
  peptide_df[[regex_pattern]] <- PTM_store
  return(peptide_df)
}

##create df for stacked peptide plot
create_PTM_df_stacked <- function(peptide_df, protein, regex_pattern, stacked_df){
  peptide_df <- peptide_df[peptide_df$intensity_store_mod > 0,]
  stacked_df$intensity_store_mod <-0
  stacked_df$mod_position <- NA
  stacked_df$PTM <- regex_pattern
  for(i in 1:nrow(stacked_df)){
    for(j in 1:nrow(peptide_df)){
      mod_position_vec<-vector()
      mod <- peptide_df[[regex_pattern]][j]
      mod_strsplit <- strsplit(as.character(mod[[1]]), "")
      if(stacked_df$peptide[i] == peptide_df$peptide_store[j]){
        for(k in 1:length(mod_strsplit)){
          if (mod_strsplit[k]=="x"){
            mod_position_vec <- append(mod_position_vec, as.character(stacked_df$start[i] + k - 1))
          }
        }
        stacked_df$intensity_store_mod[i] <- stacked_df$intensity_store_mod[i] + peptide_df$intensity_store_mod[j]
        if(length(mod_position_vec)>1){
          stacked_df$mod_position[i] <-mod_position_vec[1]
          for(l in 2:length(mod_position_vec)){
            stacked_df %>% add_row(mod_position = mod_position_vec[l],
                                   peptide=stacked_df$peptide[i],
                                   intensity_value = stacked_df$intensity_value[i],
                                   y_val = stacked_df$y_val[i], intensity_store_mod = stacked_df$intensity_store_mod[i],
                                   PTM = stacked_df$PTM[i])
          }
        }else{
          stacked_df$mod_position[i] <- mod_position_vec
        }
      }
    }
  }
 stacked_df$mod_position <- as.numeric(stacked_df$mod_position)
 stacked_df <- stacked_df[stacked_df$intensity_value > 0,]

  return(stacked_df)
}


##create new df the length of #AA in protein, including PTMs, intensity, etc.
create_PTM_plot_df_PEAKS <- function(peptide_df, protein, regex_pattern, intensity_vector_input, origin_pep_vector){
  vector_list <- list()
  vector_list_inten <- list()
  split_AA <- str_split(protein, "")
  AA_df <- data.frame(AA = split_AA[[1]])
  AA_df$AA_index <- 1:nrow(AA_df)
  for(i in 1:nrow(peptide_df)){
    PTM_vector <- rep("", nchar(protein))
    PTM <- peptide_df[[regex_pattern]][i]
    peptide <- peptide_df$peptide_store[i]
    intensity_value <-rep(0, nchar(peptide))
    PTM_strsplit <- strsplit(PTM[[1]], "")
    for(k in 1:length(PTM_strsplit)){
      if(PTM_strsplit[k]=="x")
        intensity_value[k]<- peptide_df$intensity_store_mod[i]
    }
    matches_df <-str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    intensity_vector <- rep(0, nchar(protein))
    if(matches_count >0){
    for(j in 1:matches_count){
      start <- matches_df[[j,1]]
      end <- matches_df[[j,2]]
      PTM_vector[start:end] <-PTM[[1]]
      intensity_vector[start:end] <- intensity_value
  
    }
    }
    vector_list[[i]] <- PTM_vector
    vector_list_inten[[i]] <- intensity_vector
    }
  PTM_df <- as.data.frame(sapply(vector_list, unlist))
  PTM_vec <- apply(PTM_df, 1, paste, collapse="")
  intensity_df <- as.data.frame(sapply(vector_list_inten, unlist))
  intensity_sum <- rowSums(intensity_df)
  AA_df$intensitymod <- intensity_sum
  AA_df[[regex_pattern]] <- PTM_vec
  AA_df$intensity <- intensity_vector_input
  AA_df$origin_pep <- origin_pep_vector
  if(all(AA_df[[regex_pattern]]=="")){
    stop("One or more PTM values not found in your data.")
  }
  PTM_df2 <- subset(AA_df, AA_df[[regex_pattern]] != "")
  PTM_df2 <- PTM_df2[-4]
  PTM_df2$PTM <- regex_pattern
  return(PTM_df2)
}


  
#add PTM layer to plots 
##add a new layer for eachPTM
add_PTM_layer_origin <- function(plot, PTM_df_list){
  
   plot <- plot +
   geom_point(data = PTM_df_list, shape = 16, aes(x = AA_index, y = intensity, text = paste0("AA: ",AA), 
                                           pep_label = origin_pep, intensity_label = intensitymod/intensity,
                                      fill = PTM),
              size = 2, stroke = NA) 
  return(plot)
}




add_PTM_layer <- function(plot, PTM_df_list){

plot <- plot +
  geom_point(data = PTM_df_list, shape = 16, aes(x = AA_index, y = intensity, 
                                          text = paste0("AA: ",AA) ,
                                          intensity_label = intensitymod/intensity,
                                          fill = PTM),
             size = 2, stroke = NA)
  return(plot)
}

#
add_PTM_layer_stacked <- function(plot, PTM_df){

    plot <- plot + 
      geom_point(data = PTM_df, shape = 16, aes(x = mod_position, y = y_val,
                                         text = peptide, 
                                         intensity_label = intensity_store_mod/intensity_value,
                                         fill = PTM), 
                 size =2, stroke = NA) 
  return(plot)
}

add_PTM_layer_stacked_inten <- function(plot, PTM_df){
    plot <- plot + 
      geom_point(data = PTM_df, shape = 16, aes(x = mod_position, y = intensity_value, text = peptide, 
                                         intensity_label = intensity_store_mod/intensity_value,
                                         fill = PTM), 
                 size =2, stroke = NA) 
  return(plot)
}



