read.shimadzu = function(filepath) {
  library(dplyr)
  con = file(filepath, "r")
  mcpeaktable = list()
  simsearch = list()
  filelist = list()
  f = 0
  alreadyonheader = FALSE
  while ( TRUE ) {
    if(!alreadyonheader){
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) { break } 
    }
    if(line=="[Header]") {
      alreadyonheader=FALSE
      f = f + 1
      line = readLines(con, n = 1)
      filelist[[f]] <- basename(gsub("\\\\", "/", strsplit(line, "\t")[[1]][2]) )
    } 
    if(line=="[MC Peak Table]") {
      line = readLines(con, n = 1) # # of Peaks
        npeaks <- as.integer(strsplit(line[[1]], "\t", fixed=T)[[1]][[2]])
        line = readLines(con, n = 1) #Mass	TIC
        line = readLines(con, n = npeaks+1)#include the table header
        mcpeaktable[[f]] = read.delim(text=line)
        line = readLines(con, n = 1) #empty line after table
    }
    if(line == "[MS Similarity Search Results for Spectrum Process Table]") { 
      line = readLines(con, n=1) # # of spectra
      line = readLines(con, n=1) # table header
      l = 1
      simsearchlist = list()
      while (length(line) != 0 && line != "[Header]" ){
        simsearchlist[[l]] <- line
        l = l+1 
        line <- readLines(con, n=1)
      }
      simsearchstring <- paste(unlist(simsearchlist), collapse = "\n")
      simsearch[[f]] = read.delim(text=simsearchstring)
      names(simsearch[[f]])[[5]] <- "Synonyms"
      if(length(line) != 0 && line=="[Header]") alreadyonheader= TRUE
      if ( length(line) == 0 ) { break }
    }
  }
  close(con)
  names(mcpeaktable) <- unlist(filelist)
  simsearch1 <- lapply(simsearch, FUN = function(tab) tab[tab$Hit.. == 1,])
  mergedtable <- mapply(FUN = left_join, MoreArgs = list(by=c("Peak."="Spectrum.")), mcpeaktable, simsearch1, SIMPLIFY=FALSE)
  mergedtable <- lapply(mergedtable, mutate_if, is.factor, as.character)
  peaktable <- bind_rows(mergedtable, .id = "Filename") %>% mutate_if(is.character, as.factor)
  return(peaktable)
}


read.shimadzu.quant = function(filepath) {
  library(dplyr)
  con = file(filepath, "r")
  quanttable = list()
  filelist = list()
  f = 0
  alreadyonheader = FALSE
  while ( TRUE ) {
    if(!alreadyonheader){
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) { break } 
    }
    if(line=="[Header]") {
      alreadyonheader=FALSE
      f = f + 1
      line = readLines(con, n = 1)
      filelist[[f]] <- gsub("\\\\", "/", strsplit(line, "\t")[[1]][2]) 
    } 
    if(line=="[MS Quantitative Results]") {
      line = readLines(con, n = 1) # # of IDs
      npeaks <- as.integer(strsplit(line[[1]], "\t", fixed=T)[[1]][[2]])
      line = readLines(con, n = npeaks+1)#include the table header
      quanttable[[f]] = read.delim(text=line, na.strings="Not Identified")
      colnames(quanttable[[f]]) <- c(colnames(quanttable[[f]])[-1], "emptycol") #read.delim screws up row names
      quanttable[[f]]$emptycol <- NULL
      rownames(quanttable[[f]]) <- NULL
      line = readLines(con, n = 1) #empty line after table
    }

  }
  close(con)
  names(quanttable) <- unlist(filelist)
  peaktable <- bind_rows(quanttable, .id = "Filename") %>% 
    mutate(Dirname=dirname(Filename), Filename=basename(Filename)) %>% 
    mutate_if(is.character, as.factor)
  return(peaktable)
}
