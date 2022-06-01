#error handling function to check each db has valid sequence, area, intensity, psm
errorfun <- function(peptides) {
  if(!"sequence" %in% names(peptides) | !"PSM" %in% names(peptides) | !"Area"
     %in% names(peptides) | !"Intensity" %in% names(peptides)) {
    stop("Missing Information")
  }else if (!is.numeric(peptides$Area)) {
    stop("Area problem")
  }else if (!is.numeric(peptides$Intensity)){
    stop("Intensity problem")
  }else if(!is.numeric(peptides$PSM)){
    stop("PSM problem")
  }
  for (i in peptides$sequence){
    if(grepl("[^A-Z]",i)){
      stop("sequence problem")}
  }
}

#Function to check whether file type is correct
check_file <- function(file_name, search_engine){
  if( !"\t" %in% strsplit(readLines(file_name, n=1)[1], split="")[[1]] ) {
    if(search_engine == "MSfragger"){
      stop("For the search software MSFragger, the expected file type is .tsv")
    }
  }
  else if(!"\t" %in% strsplit(readLines(file_name, n=1)[1], split="")[[1]] ){
    if(search_engine == "MaxQuant"){
      stop("For the search software MaxQuant, the expected file type is .tsv")
    }
  }
  else if( !"," %in% strsplit(readLines(file_name, n=1)[1], split="")[[1]] ) { 
    if(search_engine == "PEAKS"){
      stop("For the search software PEAKS, the expected file type is .csv")
    }
  }
}



  
