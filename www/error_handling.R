#error handling function to check each db has valid sequence, area, intensity, psm
last_check <- function(peptides) {
  if (!is.numeric(peptides$Area)) {
    stop("Insuffient Area Data")
  }else if (!is.numeric(peptides$Intensity)){
    stop("Insufficient Intensity Data")
  }else if(!is.numeric(peptides$PSM)){
    stop("Insuffiecient PSM Data")
  }
  for (i in peptides$sequence){
    if(grepl("[^A-Z]",i)){
      stop("Insufficient Sequence Data")}
  }
}

#check combined vs individual
filetype <- function(df, type_file, search_engine){
  if (search_engine == "MSfragger"){
    if(!"Protein.Description" %in% names(df)){
      stop("You have not uploaded an MSFragger File.")
    }
  }
  if (search_engine == "MaxQuant"){
    if(!("Last.amino.acid" %in% names(df))){
      stop("You have not uploaded a MaxQuant File.")
    }
  }
  
  if (search_engine == "Metamorpheus"){
    if(!"Base.Sequence" %in% names(df)){
      stop("You have not uploaded a Metamorpheus File.")
    }
    if(length(names(df)[grepl("File.Name", names(df))]) >0){ 
      store <- list()
    for (i in df$File.Name){
      if(!i %in% store){
        store<- append(store, i)}
    }
    if(type_file == "Combined" && length(store)<= 1){
      stop("You have uploaded an individual File.")
    }
    if(type_file == "Individual" && length(store) >1){
      stop("You have uploaded a combined File.")
    }}
    else{
      if(type_file == "Individual"){
        stop("You have uploaded a combined file.")
      }
    }
  }
  store <- df[, grepl("Intensity", names(df))]
  if(type_file == "Individual" && search_engine != "Metamorpheus"){
    if(length(names(store)) >1){
      stop("You have uploaded a combined file.")
    }}
  if(type_file == "Combined" && search_engine != "Metamorpheus"){
    if(length(names(store))<=1){   
      stop("You have uploaded an individual file.")
    }}
    }

#Function to check whether file type is correct
check_file <- function(file_name, search_engine){
  if(search_engine == "MSfragger"){
    if( !"\t" %in% strsplit(readLines(file_name, n=1)[1], split="")[[1]] ){ 
      stop("For the search software MSFragger, the expected file type is .tsv")}
    
  }
  
  if(search_engine == "MaxQuant"){
    if(!"\t" %in% strsplit(readLines(file_name, n=1)[1], split="")[[1]] ){
      stop("For the search software MaxQuant, the expected file type is .tsv")}
  }
  if(search_engine == "Metamorpheus"){
    if( !"\t" %in% strsplit(readLines(file_name, n=1)[1], split="")[[1]] ){ 
      stop("For the search software MetaMorpheus, the expected file type is .tsv")}

  }
  if(search_engine == "PEAKS")  { 
    store <- strsplit(readLines(file_name, n=1)[1], split = "")[[1]]    
  if(!"," %in% store[1:15]){
      stop("For the search software PEAKS, the expected file type is .csv")
    }
  }}
  
#check that fasta file is correct
fasta_check <- function(db_file){
  store <- strsplit(readLines(db_file, n=1)[1], split = " ")[[1]]
  if(!grepl('>', substr(store[1], 1, 1))){
    stop("For the database file input, a fasta file is expected.")
  }
}





  
