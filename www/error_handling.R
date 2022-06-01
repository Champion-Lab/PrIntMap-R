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