
#check combined vs individual
filetype <- function(df, type_file, search_engine){
  store <- df[, grepl("Intensity", names(df))]
  store2 <- df[, grepl("Area", names(df))]
  store3 <- df[, grepl("PSM", names(df))]
  if (search_engine == "Metamorpheus"){
    if(length(names(df)[grepl("File.Name", names(df))]) >0){ 
      if(uniqueN(df, by = "File.Name")>1 && type_file == "Individual"){
        stop('You have uploaded a combined file.')
      }
      if(uniqueN(df, by = "File.Name") <= 1 && type_file == "Combined"){
        stop('You have uploaded an individual file.')}
     }
    else{
      if(type_file == "Individual"){
        stop("You have uploaded a combined file.")
      }
    }
  }
  if(search_engine == "Generic"){
    if(type_file == "Individual"){
      if(length(names(store)) >1 |length(names(store2)) >1 | length(names(store3)) >1){
        stop("You have uploaded a combined file.")
    }
    }
    if(type_file == "Combined"){
      if(length(names(store)) ==1 |length(names(store2)) ==1 | length(names(store3)) ==1){
        stop("You have uploaded an individual file.")
      }
    }
  }
  
  if(type_file == "Combined" && search_engine == "Proteome Discover"){
    stop("You have uploaded an Individual File.")
  }
  
  if(search_engine == "PEAKS"){
    if(type_file == "Individual"){
      if(length(names(store)) >1 |length(names(store2)) >1 | length(names(store3)) >1){
        stop("You have uploaded a combined file.")
      }
    }
    if(type_file == "Combined"){
      if(length(names(store)) ==1 |length(names(store2)) ==1 | length(names(store3)) ==1){
        stop("You have uploaded an individual file.")
      }
    }
  }
  
 
  if(type_file == "Individual" && (search_engine != "Metamorpheus" | search_engine != "Generic")){
    if(length(names(store)) >1){
      if (search_engine != "PEAKS") {
        stop("You have uploaded a combined file.")
      }
    }}
  if(type_file == "Combined" && (search_engine != "Metamorpheus" | search_engine != "Generic")){
    if(length(names(store))<=1){   
      if (search_engine != "PEAKS") {
        stop("You have uploaded an individual file.")
      }
    }}
  }

#Function to check whether file type is correct
check_file <- function(file_name, search_engine){
  store1 <- strsplit(readLines(file_name, n=1)[1], split="")[[1]]
  store2 <- strsplit(readLines(file_name, n=1)[1], split="\t")[[1]]
  store3 <- strsplit(readLines(file_name, n=1)[1], split="\"\t\"")[[1]]
  store4 <- strsplit(readLines(file_name, n=1)[1], split=",")[[1]]
  store5 <- strsplit(readLines(file_name, n =1)[1], split = "\",\"")[[1]]
  
  
  if(search_engine == "MSfragger"){
    if(!any(store1 == "\t")) { 
      stop("For the search software MSFragger, the expected file type is .tsv")}
    if(!any(store2 == "Protein Description")){
      stop("You have not uploaded an MSFragger file.")}
    }
  
  if(search_engine == "MaxQuant"){
    if(!any(store1 == "\t")){
      stop("For the search software MaxQuant, the expected file type is .tsv")}
    if(!any(store2 == "Last amino acid")){
      stop("You have not uploaded a MaxQuant File.")
    }
  }
  if(search_engine == "Metamorpheus"){
    if(!any(store1 == "\t" )){ 
      stop("For the search software MetaMorpheus, the expected file type is .tsv")}
    if(!any(store2 == "Base Sequence")){
      stop("You have not uploaded a Metamorpheus File.")
    }
  }
  if(search_engine == "Proteome Discoverer"){
    if(!any(store1 == "\t" )){ 
      stop("For the search software Proteome Discoverer, the expected file type is .tsv")}
    if(!any(store3 == "Qvality PEP")){
      stop("You have not uploaded a Proteome Discoverer File.")
    }
  }
  if(search_engine == "PEAKS")  { 
  if(!any(store1[1:15] ==",")){
      stop("For the search software PEAKS, the expected file type is .csv")
  }
    if(! (any(store4 == "RT") | any(store5 == "RT") | any(store5 == "-10LgP") | any(store4 == "-10LgP"))){
      stop("You have not uploaded a PEAKS file.")
    }
  }
  if(search_engine == "Generic"){
    if(!any(store1[1:15] ==",")){
      stop("For a generic file, the expected file type is .csv")
    }
    if(any(store4 == "Scan") | any(store5 == "Scan")){
      stop("You have not uploaded a compatible generic file")
    }
  }
  }
  
#check that fasta file is correct
fasta_check <- function(db_file){
  store <- strsplit(readLines(db_file, n=1)[1], split = " ")[[1]]
  if(!grepl('>', substr(store[1], 1, 1))){
    stop("For the database file input, a fasta file is expected.")
  }
}




  
