#import database
#returns seqinr fasta list
import_db <- function(db_file) {
  fasta_check(db_file)
  db <- read.fasta(file = db_file,
                   seqtype = "AA",
                   as.string = T,
                   whole.header = F)
  return(db)
}


#choose a protein from a database based on an accession number
#returns a protein object
select_prot <- function(db, Accession) {
  protein_names <- names(db)
  protein_index <- which(grepl(Accession, protein_names, ignore.case = T))
  if(length(protein_index) > 1) {
    stop("Multiple proteins found. Include additional characters in Accession number")
  } else if (length(protein_index) < 1) {
    stop("Protein not found. Check that acession number is correct")
  } else {
    protein <- db[[protein_index[1]]]
  }
  return(protein)
}


#choose a protein from a database based on an accession number
#returns a protein object for volcano plot
select_prot_volcano <- function(db, Accession) {
  protein_names <- names(db)
  protein_index <- which(grepl(Accession, protein_names, ignore.case = T))
  if(length(protein_index) > 1) {
    return(NULL)
  } else if (length(protein_index) < 1) {
    return(NULL)
  } else {
    protein <- db[[protein_index[1]]]
    protein_name <- attr(protein, "name")
    protein_accession <- str_split(protein_name, "\\|")[[1]][2]
  }
  return(protein_accession)
}

#generate protein attributes
#takes protein object from database
#returns vector (sequence, name, description)
protein_attributes <- function(protein) {
  protein_name <- attr(protein, "name")
  protein_seq <- protein[1]
  protein_annot <- attr(protein, "Annot")
  protein_accession <- str_split(protein_name, "\\|")[[1]][2]
  protein_attr <- c(protein_seq, protein_name, protein_annot, protein_accession)
  return(protein_attr)
}


#create amino acid sequence dataframe
#returns dataframe with each AA in its own row
#takes raw protein sequence
create_AA_df <- function(protein) {
  split_AA <- str_split(protein, "")
  AA_df <- data.frame(AA = split_AA[[1]])
  AA_df$AA_index <- 1:nrow(AA_df)
  return(AA_df)
}


#import peptide file from peaks
#takes csv file in peaks format and returns dataframe
#options to add a sample name, and a filter which removes any proteins that contain
###that string in the Accession column
read_peptide_csv_PEAKS_bysamp <- function(peptide_file, sample = NA, filter = NA) {
  check_file(peptide_file, "PEAKS")
  peptides <- read.csv(peptide_file)
  filetype(peptides, "Individual", "PEAKS")
  peptides$sequence <- str_remove_all(peptides$Peptide, "[a-z1-9()+-:.]")
  names(peptides)[grepl("Area", names(peptides))] <- "Area"
  names(peptides)[grepl("Intensity", names(peptides))] <- "Intensity"
  names(peptides)[grepl("Spec", names(peptides))] <- "PSM"
  if (!is.na(sample)) {
    peptides$sample <- sample
  }
  if (!is.na(filter)) {
    peptides <- peptides[grepl(filter, peptides$Accession) == F,]
  }
  return_list <- list(peptides, 1)
  return(return_list)
}


#import peptide file from PEAKS (combined, not individual sample)
#takes tsv file in fragpipe format and returns dataframe
#options to add a sample name, and a filter which removes any proteins that contain
###that string in the Accession column
read_peptide_csv_PEAKS_comb <- function(peptide_file, sample_pattern, sample = NA, filter = NA, comb_method = "Sum") {
  check_file(peptide_file, "PEAKS")
  peptide_import <- read.csv(peptide_file)
  filetype(peptide_import, "Combined", "PEAKS")
  names(peptide_import)[names(peptide_import) == "X.Spec"] <- "total_spectra"
  peptides <- peptide_import
  peptides$sequence <- str_remove_all(peptides$Peptide, "[a-z1-9()+-:.]")
  
  if(length(names(peptides)[grepl(sample_pattern, names(peptides))])<=0){
    stop("Sample Pattern not found in file.")
  }
  else {
  PSM_pattern <- paste0("X\\.Spec", ".*", sample_pattern, ".*")
  PSM_df <- peptide_import[,grepl(PSM_pattern, names(peptide_import))]
  PSM_df[is.na(PSM_df)] <- 0
  PSM_vec <- rowSums(as.data.frame(PSM_df), na.rm = T)
  sample_count <- ncol(as.data.frame(PSM_df))
  
  Intensity_pattern <- paste0("Intensity", ".*", sample_pattern, ".*")
  Intensity_df <- peptide_import[,grepl(Intensity_pattern, names(peptide_import))]
  Intensity_df[is.na(Intensity_df)] <- 0
  Intensity_vec <- rowSums(as.data.frame(Intensity_df), na.rm = T)
  
  Area_pattern <- paste0("Area", ".*", sample_pattern, ".*")
  Area_df <- peptide_import[,grepl(Area_pattern, names(peptide_import))]
  Area_df[is.na(Area_df)] <- 0
  Area_vec <- rowSums(as.data.frame(Area_df), na.rm = T)
  
  if (comb_method == "Average") {
    PSM_vec <- PSM_vec / sample_count
    Intensity_vec <- Intensity_vec / sample_count
    Area_vec <- Area_vec / sample_count
  }
  
  peptides$PSM <- PSM_vec
  peptides$Intensity <- Intensity_vec
  peptides$Area <- Area_vec
  if (!is.na(sample)) {
    peptides$sample <- sample
  }
  if (!is.na(filter)) {
    peptides <- peptides[grepl(filter, peptides$Accession) == F,]
  }

  return_list <- list(peptides, sample_count)
  return(return_list)
  }}

#import peptide file from generic
#takes csv file in generic format and returns dataframe
read_peptide_csv_generic_bysamp <- function(peptide_file) {
  #check_file(peptide_file, "PEAKS")
  peptides <- read.csv(peptide_file)
  #filetype(peptides, "Individual", "PEAKS")
  peptides$sequence <- str_remove_all(peptides$Peptide, "[a-z1-9()+-:.]")
  names(peptides)[grepl("Area", names(peptides))] <- "Area"
  names(peptides)[grepl("Intensity", names(peptides))] <- "Intensity"
  names(peptides)[grepl("Spec", names(peptides))] <- "PSM"
  return_list <- list(peptides, 1)
  return(return_list)
}

#import peptide file from generic, multiple samples
#takes csv file in generic format and returns dataframe
read_peptide_csv_generic_comb <- function(peptide_file, sample_pattern, comb_method = "Sum") {
  #check_file(peptide_file, "PEAKS")
  peptide_import <- read.csv(peptide_file)
  #filetype(peptide_import, "Combined", "PEAKS")
  peptides <- peptide_import
  peptides$sequence <- str_remove_all(peptides$Peptide, "[a-z1-9()+-:.]")
  
  if(length(names(peptides)[grepl(sample_pattern, names(peptides))])<=0){
    stop("Sample Pattern not found in file.")
  }
  else {
    if (length(names(peptides)[grepl("PSM", names(peptides))]) > 0) {
      PSM_pattern <- paste0(sample_pattern, ".PSM")
      PSM_df <- peptide_import[,grepl(PSM_pattern, names(peptide_import))]
      PSM_df[is.na(PSM_df)] <- 0
      PSM_vec <- rowSums(as.data.frame(PSM_df))
      sample_count <- ncol(as.data.frame(PSM_df))
      if (comb_method == "Average") {
        PSM_vec <- PSM_vec / sample_count
      }
      peptides$PSM <- PSM_vec
    }
    if (length(names(peptides)[grepl("Intensity", names(peptides))]) > 0) {
      Intensity_pattern <- paste0(sample_pattern, ".Intensity")
      Intensity_df <- peptide_import[,grepl(Intensity_pattern, names(peptide_import))]
      Intensity_df[is.na(Intensity_df)] <- 0
      Intensity_vec <- rowSums(as.data.frame(Intensity_df))
      sample_count <- ncol(as.data.frame(Intensity_df))
      if (comb_method == "Average") {
        Intensity_vec <- Intensity_vec / sample_count
      }
      peptides$Intensity <- Intensity_vec
    }
    if (length(names(peptides)[grepl("Area", names(peptides))]) > 0) {
      Area_pattern <- paste0(sample_pattern, ".Area")
      Area_df <- peptide_import[,grepl(Area_pattern, names(peptide_import))]
      Area_df[is.na(Area_df)] <- 0
      Area_vec <- rowSums(as.data.frame(Area_df))
      sample_count <- ncol(as.data.frame(Area_df))
      if (comb_method == "Average") {
        Area_vec <- Area_vec / sample_count
      }
      peptides$Area <- Area_vec
    }
    return_list <- list(peptides, sample_count)
    return(return_list)
  }}

#import peptide file from MSFragger (not combined, individual sample)
#takes tsv file in fragpipe format and returns dataframe
#options to add a sample name, and a filter which removes any proteins that contain
###that string in the Accession column
read_peptide_tsv_MSFragger_bysamp <- function(peptide_file, sample = NA, filter = NA) {
  check_file(peptide_file, "MSfragger")
  peptides <- read.csv(peptide_file, sep = "\t", header = T)
  filetype(peptides, "Individual", "MSfragger")
  names(peptides)[names(peptides) == "Peptide" | names(peptides) == "Peptide.Sequence"] <- "sequence"
  names(peptides)[grepl("Intensity", names(peptides))] <- "Intensity"
  names(peptides)[grepl("Spec", names(peptides))] <- "PSM"
  peptides <- peptides[peptides$PSM > 0,]
  if (length(names(peptides)[grepl("LFQ", names(peptides))]) >0){
    names(peptides)[grepl("LFQ", names(peptides))] <- "Area"
  }
  if (!is.na(sample)) {
    peptides$sample <- sample
  }
  if (!is.na(filter)) {
    peptides <- peptides[grepl(filter, peptides$Accession) == F,]
  }
  return_list <- list(peptides, 1)
  return(return_list)
}


#import peptide file from MSFragger (combined, not individual sample)
#takes tsv file in fragpipe format and returns dataframe
#options to add a sample name, and a filter which removes any proteins that contain
###that string in the Accession column
read_peptide_tsv_MSFragger_comb <- function(peptide_file, sample_pattern, sample = NA, filter = NA, comb_method = "Sum") {
  check_file(peptide_file, "MSfragger")
  peptide_import <- read.csv(peptide_file, sep = "\t", header = T, )
  filetype(peptide_import, "Combined", "MSfragger")
  peptides <- peptide_import
  names(peptides)[names(peptides)== "Sequence" | names(peptides) == "Peptide.Sequence"] <- "sequence"
  if(length(names(peptides)[grepl(sample_pattern, names(peptides))])<=0){
    stop("Sample Pattern not found in file.")
  }
  else{
    
 
  PSM_pattern <- paste0(".*", sample_pattern, ".*", "Spectral\\.Count")
  PSM_df <- peptide_import[,grepl(PSM_pattern, names(peptide_import))]
  PSM_df[is.na(PSM_df)] <- 0
  PSM_vec <- rowSums(as.data.frame(PSM_df), na.rm = T)
  sample_count <- ncol(as.data.frame(PSM_df))
  
  Intensity_pattern <- paste0(".*", sample_pattern, ".*", "Intensity")
  Intensity_df <- peptide_import[,grepl(Intensity_pattern, names(peptide_import))]
  Intensity_df[is.na(Intensity_df)] <- 0
  Intensity_vec <- rowSums(as.data.frame(Intensity_df), na.rm = T)
  
  if (length(names(peptides)[grepl("LFQ", names(peptides))]) >0){
  Area_pattern <- paste0(".*", sample_pattern, ".*", "MaxLFQ.Intensity")
  Area_df <- peptide_import[,grepl(Area_pattern, names(peptide_import))]
  Area_vec <- rowSums(as.data.frame(Area_df), na.rm = T)
  peptides$Area <- Area_vec
  if (comb_method == "Average") {
    Area_vec <- Area_vec / sample_count
  }
  }
  
  if (comb_method == "Average") {
    PSM_vec <- PSM_vec / sample_count
    Intensity_vec <- Intensity_vec / sample_count
  }
  peptides$PSM <- PSM_vec
  peptides$Intensity <- Intensity_vec
  peptides <- peptides[peptides$PSM > 0,]

  
  if (!is.na(sample)) {
    peptides$sample <- sample
  }
  if (!is.na(filter)) {
    peptides <- peptides[grepl(filter, peptides$Accession) == F,]
  }
  return_list <- list(peptides, sample_count)
  return(return_list)
}}


#import peptide file from MaxQuant (combined, not individual sample)
#takes tsv file in MaxQuant peptides.txt format and returns dataframe
read_peptide_tsv_MaxQuant_comb <- function(peptide_file, sample_pattern, sample = NA, filter = NA, comb_method = "Sum") {
  
  check_file(peptide_file, "MaxQuant")
  peptide_import <- read.csv(peptide_file, sep = "\t", header = T, )
  filetype(peptide_import, "Combined", "MaxQuant")
  names(peptide_import)[names(peptide_import) == "Intensity"] <- "summed_intensity"
  peptides <- peptide_import
  peptides$sequence <- peptides$Sequence
  
  if(length(names(peptides)[grepl(sample_pattern, names(peptides))])<=0){
    stop("Sample Pattern not found in file.")
  }
  else{
    

  PSM_pattern <- paste0("Experiment", ".*", sample_pattern, ".*")
  PSM_df <- peptide_import[,grepl(PSM_pattern, names(peptide_import))]
  PSM_df[is.na(PSM_df)] <- 0
  PSM_vec <- rowSums(as.data.frame(PSM_df), na.rm = T)
  sample_count <- ncol(as.data.frame(PSM_df))
  
  Intensity_pattern <- paste0("Intensity", ".*", sample_pattern, ".*")
  Intensity_df <- peptide_import[,grepl(Intensity_pattern, names(peptide_import))]
  Intensity_df[is.na(Intensity_df)] <- 0
  Intensity_vec <- rowSums(as.data.frame(Intensity_df), na.rm = T)
  
  Area_pattern <- paste0("LFQ.intensity", ".*", sample_pattern, ".*")
  Area_df <- peptide_import[,grepl(Area_pattern, names(peptide_import))]
  Area_df[is.na(Area_df)] <- 0
  Area_vec <- rowSums(as.data.frame(Area_df), na.rm = T)
  
  if (comb_method == "Average") {
    PSM_vec <- PSM_vec / sample_count
    Intensity_vec <- Intensity_vec / sample_count
    Area_vec <- Area_vec / sample_count
  }
  
  peptides$PSM <- PSM_vec
  peptides$Intensity <- Intensity_vec
  peptides$Area <- Area_vec
  peptides <- peptides[peptides$PSM > 0,]
  
  if (!is.na(sample)) {
    peptides$sample <- sample
  }
  if (!is.na(filter)) {
    peptides <- peptides[grepl(filter, peptides$Accession) == F,]
  }
  return_list <- list(peptides, sample_count)
  return(return_list)
}}


#import peptide file from MetaMorpheus (individual Sample)
#takes tsv file in Metamorpheus  format and returns dataframe
read_peptide_tsv_Metamorpheus_bysamp <- function(peptide_file, sample = NA, filter = NA) {
  check_file(peptide_file, "Metamorpheus")
  peptides <- read.csv(peptide_file, sep = "\t", header = T)
  filetype(peptides, "Individual", "Metamorpheus")
  peptides$sequence <- peptides$Base.Sequence
  names(peptides)[grepl("Total.Ion.Current", names(peptides))] <- "Intensity"
  names(peptides)[grepl("PSM.Count.*", names(peptides))] <- "PSM"
  peptides <- peptides[peptides$PSM > 0,]
  if (!is.na(sample)) {
    peptides$sample <- sample
  }
  if (!is.na(filter)) {
    peptides <- peptides[grepl(filter, peptides$Accession) == F,]
  }
  return_list <- list(peptides, 1)
  return(return_list)
}

#import peptide file from MetaMorpheus (combined Sample)
#takes tsv file in Metamorpheus  format and returns dataframe
read_peptide_tsv_Metamorpheus_comb <- function(peptide_file, sample_pattern, sample = NA, filter = NA, comb_method = "Sum") {
  check_file(peptide_file, "Metamorpheus")
  peptides <- read.csv(peptide_file, sep = "\t", header = T)
  filetype(peptides, "Combined", "Metamorpheus")
  if(length(names(peptides)[grepl("Total.Ion.Current", names(peptides))]) >0){
    if(!(any(grepl(sample_pattern, peptides$File.Name)))){
      stop("Sample pattern not found in file.")
    }
    else{
    peptides$sequence <- peptides$Base.Sequence
    names(peptides)[grepl("Total.Ion.Current", names(peptides))] <- "Intensity"
    names(peptides)[grepl("PSM.Count.*", names(peptides))] <- "PSM"
    peptides <- peptides[peptides$PSM > 0,]
    if (!is.na(sample)) {
      peptides$sample <- sample
    }
    if (!is.na(filter)) {
      peptides <- peptides[grepl(filter, peptides$Accession) == F,]
    }
    peptides <- peptides[str_detect(peptides$File.Name, paste0(".*", sample_pattern, ".*")),]
    
    if (comb_method == "Average") {
      peptide_summary <- peptides %>% group_by(sequence) %>% summarise_if(is.numeric, mean, na.rm = TRUE)
    } else {
      peptide_summary <- peptides %>% group_by(sequence) %>% summarise_if(is.numeric, sum, na.rm = TRUE)
    }
    
    sample_count <- length(unique(peptides$File.Name))
    return_list <- list(peptide_summary, sample_count)
    return(return_list) 
  }}
 else{
   if(length(names(peptides)[grepl(sample_pattern, names(peptides))])<=0){
     stop("Sample Pattern not found in file.")
   }
   else{
   names(peptides)[names(peptides) == "Sequence"] <- "Full.Sequence"
   peptides$sequence <- peptides$Base.Sequence
   
   Area_pattern <- paste0("Intensity_", ".*", sample_pattern, ".*")
   Area_df <- peptides[,grepl(Area_pattern, names(peptides))]
   Area_df[is.na(Area_df)] <- 0
   Area_vec <- rowSums(as.data.frame(Area_df), na.rm = T)
   sample_count <- ncol(as.data.frame(Area_df))
   
   if (comb_method == "Average") {
     Area_vec <- Area_vec / sample_count
   }
   
   peptides$Area <- Area_vec
   if (!is.na(sample)) {
     peptides$sample <- sample
   }
   if (!is.na(filter)) {
     peptides <- peptides[grepl(filter, peptides$Accession) == F,]
   }
   return_list <- list(peptides, sample_count)
   return(return_list)
 }}
}


#import peptide file from Proteome Discover (individual Sample)
#takes tsv file in PD  format and returns dataframe
read_peptide_tsv_ProteomeDiscover_bysamp <- function(peptide_file, sample = NA, filter = NA) {
  check_file(peptide_file, "Proteome Discoverer")
  peptides <- read.csv(peptide_file, sep = "\t", header = T)
  peptides$sequence <- peptides$Sequence
  names(peptides)[grepl("X..PSMs", names(peptides))] <- "PSM"
  if (!is.na(sample)) {
    peptides$sample <- sample
  }
  if (!is.na(filter)) {
    peptides <- peptides[grepl(filter, peptides$Accession) == F,]
  }
  return_list <- list(peptides, 1)
  return(return_list)
}



#generates appropriate choices for intensity metric based on uploaded data file
#returns list of choices
intensity_metric_choices <- function(df1){
  choices <- list()
  options <- c("PSM","Area", "Intensity")
  for(i in options){
    if(i %in% names(df1)){
      choices <- append(choices, i)
    }
  }
  return(choices)
}


#create intensity dataframe
#X.Spec = psms
#returns vector with intensities the same length as # of AA
create_intensity_vec <- function(peptide_df,
                                 protein,
                                 intensity = "PSM") {
  vector_list <- list()
  if(intensity %in% names(peptide_df)){
    for (i in 1:nrow(peptide_df)) {
    peptide <- peptide_df$sequence[i]
    intensity_value <- peptide_df[[intensity]][i]
    
    if (is.na(intensity_value)){
      intensity_value <- 0
    }
    matches_df <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    intensity_vector <- rep(0, nchar(protein))
    if (matches_count > 0) {
      for (j in 1:matches_count) {
        start <- matches_df[[j,1]]
        end <- matches_df[[j,2]]
        intensity_vector[start:end] <- intensity_value
      }
    }
    vector_list[[i]] <- intensity_vector
  }
    intensity_df <- as.data.frame(sapply(vector_list, unlist))
    intensity_sum <- rowSums(intensity_df)
    return(intensity_sum)}
  else{
    stop("Intensity metric not found in one or more of your files.")
  }
  
}

#combine intensity vector with amino acid dataframe for plotting
#returns dataframe that can be plotted
combine_AA_intensity <- function(AA_df, intensity_vec) {
  AA_df$intensity <- intensity_vec
  return(AA_df)
}


#create simple plot
#returns plot object
plot_intensity <- function(AA_df, protein, intensity_label = "PSM",
                           font_size = 15,color = "grey18",
                           alpha = 1, plot_type = "line") {
  if (plot_type == "bar") {
    plot <- ggplot(data = AA_df) +
      geom_bar(aes(x = AA_index, y = intensity),
               stat = "identity", position = "dodge",
               fill = color, alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, color = AA),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank()) +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label) +
      scale_x_continuous(breaks = pretty_breaks(10))
  } else {
    plot <- ggplot(data = AA_df) +
      geom_line(aes(x = AA_index, y = intensity),
                color = color, alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, color = AA),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
      scale_x_continuous(breaks = pretty_breaks(10))
  }
  return(plot)
}


#create interactive plot object
create_plotly <- function(plot) {
  plotly <- ggplotly(plot, tooltip = c("AA", "x", "y", "sample", "origin_pep", "repeated",
                                       "fold_change_label"))
  show(plotly)
  return(plotly)
}



#create origin peptide vector
create_peptide_origin_vec <- function(peptide_df,
                                      protein,
                                      intensity = "PSM") {
  vector_list <- list()
  for (i in 1:nrow(peptide_df)) {
    peptide <- peptide_df$sequence[i]
    intensity_value <- peptide_df[[intensity]][i]
    matches_df <- str_locate_all(protein, peptide)[[1]]
    matches_count <- nrow(matches_df)
    origin_vector <- rep("", nchar(protein))
    if (matches_count > 0) {
      for (j in 1:matches_count) {
        start <- matches_df[[j,1]]
        end <- matches_df[[j,2]]
        origin_vector[start:end] <- paste0(peptide, "|")
      }
    }
    vector_list[[i]] <- origin_vector
  }
  origin_df <- as.data.frame(sapply(vector_list, unlist))
  origin_vec <- apply(origin_df, 1, paste, collapse="")
  
  #create unique seq counts
  origin_vector_unique <- rep("", length(origin_vec))
  for (i in 1:length(origin_vec)) {
    if (origin_vec[i] != "") {
      split_str <- str_split(origin_vec[i], "\\|")[[1]]
      split_str <- split_str[split_str != ""]
      unique_peps <- unique(split_str)
      final_string <- paste(unique_peps, collapse = "\n", sep = "")
      origin_vector_unique[i] <- final_string
    }
  }
  return(origin_vector_unique)
}

#combine peptide origin vector with AA_df
add_origin_peptide_vector <- function(AA_df,
                                      origin_vec) {
  AA_df$origin_pep <- origin_vec
  return(AA_df)
}


#plot origin plot
plot_origin <- function(AA_df, protein,
                        intensity_label = "PSM",
                        font_size = 15, 
                        color = "grey18", alpha = 1, plot_type = "line") {
  if (plot_type == "bar") {
    plot <- ggplot(data = AA_df) +
      geom_bar(aes(x = AA_index, y = intensity),
               stat = "identity", position = "dodge",
               fill = color, alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, color = AA, fill = origin_pep),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
      scale_x_continuous(breaks = pretty_breaks(10))
  } else {
    plot <- ggplot(data = AA_df) +
      geom_line(aes(x = AA_index, y = intensity),
                color = color, alpha = alpha) +
      geom_point(aes(x = AA_index, y = intensity, color = AA, fill = origin_pep),
                 size = 0, alpha = 0) +
      theme_bw(base_size = font_size) +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(title = protein, x = "Amino Acid Number", y = intensity_label)+
      scale_x_continuous(breaks = pretty_breaks(10))
  }
  
  show(plot)
  return(plot)
}

#calculate percent coverage based on an intensity dataframe
#returns list that has: nAA - total number of Amino Acids
#                       covered_AA - number of Amino Acids observed
#                       percent_coverage - percent coverage
calculate_percent_cov <- function(intensity_df) {
  nAA <- nrow(intensity_df)
  covered_AA <- sum(intensity_df$intensity > 0)
  percent_coverage <- round((covered_AA / nAA) * 100, 2)
  return(list(nAA = nAA,
              covered_AA = covered_AA,
              percent_coverage = percent_coverage))
  
}


#a function to replace the newline character with ; in origin peptide column
covert_to_download_df <- function(AA_df, intensity_name) {
  AA_df$origin_pep <- gsub("[\r\n]", ";", AA_df$origin_pep)
  names(AA_df)[3] <- intensity_name
  return(AA_df)
}

#to format y-axis for log scale

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(c(1, max(x, na.rm = TRUE))), log = TRUE, n = n)
  }
}