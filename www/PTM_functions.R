
#function to pull out PTMs from PEAKS files
#returns vector containing all PTMS the same length as number of amino acids
create_PTM_vec_PEAKS <- function(peptide_df, protein){
  PTM_list <- list()
  for(i in 1:nrow(peptide_df)){
    mod_peptide<-peptide_df$Peptide[i]
    peptide <- peptide_df$sequence[i]
    PTM_df <-str_locate_all(mod_peptide, "\\(.?.?.?.?.?.?\\)")[[1]]
    PTM_count <- nrow(PTM_df)
    PTM_value <- vector()
    mod_peptide_strsplit <- strsplit(mod_peptide, "")[[1]]
    for(j in 1: length(mod_peptide_strsplit)){
      if (!grepl("[a-z1-9()+-:.]", mod_peptide_strsplit[j])){
            PTM_value <- c(PTM_value, "F")
            for (k in 1:PTM_count) {
                start2 <- PTM_df[[k,1]]
                end2 <- PTM_df[[k,2]]
                if(j == start-1){
                PTM_value<- PTM_value[-length(PTM_value)]
                PTM_value<- c(PTM_value, substr(mod_peptide, start2, end2))
                }
            }
        PTM_list <- append(PTM_list, list(PTM_value))    
      }
    }}
  peptide_df$PTMs <- PTM_list
  
  return(peptide_df)}
  
create_PTM_vec <- function(pattern, protein, peptide_df) {
  PTM_vector <- rep("", nchar(protein))
  for(i in 1:nrow(peptide_df)){
    peptide <- peptide_df$sequence[i] 
    PTM <- peptide_df$PTMs[i]
    matches_df <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    if(matches_count>0){
      for(j in 1:nrow(matches_df)){
        start <- matches_df[[j,1]]
        end <- matches_df[[j,2]]
        PTM_vector[start:end] <- PTM
        
      }
    }
  }
  
 
  split_AA <- str_split(protein, "")
  AA_df <- data.frame(AA = split_AA[[1]])
  AA_df$AA_index <- 1:nrow(AA_df)
    motifs <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(motifs)
    motif_vector <- rep(F, nrow(AA_df))
    if (matches_count > 0) {
      for (j in 1:matches_count) {
        start <- motifs[[j,1]]
        end <- motifs[[j,2]]
        motif_vector[start:end] <- T
      }}
  return(motif_vector) 
} 


        

#combine PTM vector with AA_df
create_PTM_df <- function(AA_df,
                          PTM_vec) {
  AA_df$PTMs <- PTM_vec
  # print(AA_df)
  return(AA_df)
}


#function to refine PTM df based on user input
#accepts AA_df and returns df with count of specified PTM per AA and column containing 
#regex PTM pattern

refine_PTM_df <- function(AA_df, PTM_pattern, PTM_name){
    for(j in 1:nrow(AA_df)){
      store_vec <- vector()
      count_vec <- vector()
      PTM <- AA_df$PTMs[j]
      matches_df <- str_locate_all(PTM, PTM_pattern)[[1]]
      # print(matches_df)
      matches_count <- nrow(matches_df)
      count_vec <- c(count_vec, matches_count)
      if(matches_count >0){
        store_vec <- c(store_vec, TRUE)
      }else{
        store_vec <- c(store_vec, FALSE)
      }
    }
    AA_df[[paste0(PTM_name, " Count")]] <- count_vec
    AA_df[[PTM_name]] <- store_vec
    PTM_df <- AA_df[AA_df$PTM_names == T,]
    
  return(PTM_df)
}


#add PTM layer to plots 
##add a new layer for eachPTM
add_PTM_layer <- function(plot, PTM_df, color = "red" ){
  plot <- plot +
    geom_vline(xintercept = PTM_df$AA_index, color = color,
               size = 0.5)
  return(plot)
}