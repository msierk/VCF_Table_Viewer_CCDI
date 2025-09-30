# Overwrite vcfR extract.info() function to account for '=' appearing
extract_info_UP <- function(x, element, as.numeric=FALSE, mask=FALSE){
  
  #if( class(x) == 'chromR' ){
  if( inherits(x, 'chromR') ){
    mask <- x@var.info$mask
    x <- x@vcf
  }
  #if( class(x) != 'vcfR' ){
  if( !inherits(x, 'vcfR') ){
    stop("Expecting an object of class vcfR or chromR.")
  }
  
  #  values <- unlist(
  #    lapply(strsplit(unlist(
  #      lapply(strsplit(x@fix[,'INFO'], split=";"),
  #             function(x){grep(paste("^", element, "=", sep=""), x, value=TRUE)})),
  #      split="="), function(x){x[2]})
  #  )
  
  values <- strsplit(x@fix[,'INFO'], split=";")
  values <- lapply(values, function(x){grep(paste("^", element, "=", sep=""), x, value=TRUE)})
  values <- lapply(values, function(x){ unlist( strsplit(x, split="=") ) })
  values <- lapply(values, function(x){if(length(x) == 2) { x[2] } else { paste0(x[2:length(x)],collapse="=")}})
  values <- lapply(values, function(x){ if(is.null(x)){NA}else{x} })
  values <- unlist(values)
  
  
  if(as.numeric == TRUE){
    values <- as.numeric(values)
  }
  #  if( mask != FALSE & !is.null(mask) ){
  #    values <- values[x@var.info$mask]
  #    values <- values[mask]
  #  }
  values
}

# Overwrite vcfR INFO2df() function to use updated extract_info_UP function
INFO2df_UP <- function(x){
  
  if( inherits(x, "chromR") ){
    x <- x@vcfR
  }  
  
  metaINFO <- metaINFO2df(x)
  
  # Initialize a data.frame for the INFO data
  INFOdf <- data.frame( matrix( nrow=nrow(x@fix),  ncol=nrow(metaINFO) ) )
  names(INFOdf) <- metaINFO[,'ID']
  
  for( i in 1:nrow(metaINFO) ){
    tmp <- extract_info_TJO(x, element = metaINFO[,'ID'][i])
    if( metaINFO[,'Type'][i] == "Integer" & metaINFO[,'Number'][i] == "1" ){
      tmp <- as.integer( tmp )
    }
    if( metaINFO[,'Type'][i] == "Float" & metaINFO[,'Number'][i] == "1" ){
      tmp <- as.numeric( tmp )
    }
    INFOdf[,metaINFO[,'ID'][i]] <- tmp
  }
  return(INFOdf)
}