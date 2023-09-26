library(xml2)
library(bit64)
library(tools)
library(fuzzyjoin)
library(tidyverse)

#Markes sequence aligner
#Author: John Powers
#Purpose: match up the Markes sequence and the list of Shimadzu chromatogram files to infer skips
#fuzzy joins based on the desorb file of the Markes traps and the creation time of the qgd files

# Markes sequence ---------------------------------------------------------

#Convert from Markes timestamp to UTC
markestime  <- function(x) as.POSIXct(as.POSIXlt(as.integer64(x)/10000000, origin = "0001-01-01",tz="UTC"))

#Fill in gaps in a list of desorb start times by adding a time at the midpoint
fill_starttime <- function(ds) {
  ds[is.na(ds)] <- sapply(which(is.na(ds)), "+", c(-1, 1)) %>% #get matrix of indices of preceding and following entries
    apply(2, function(x) ifelse(sum(x)>2 && as.double(diff(ds[x]), units = "mins")>50, mean(ds[x]), NA)) %>% #if 50 min apart, take midpoint
    unlist() %>% as.POSIXct(origin = lubridate::origin)
  ds
}

#Convert Markes sequence from XML format to a dataframe
#From https://stackoverflow.com/questions/34273132/r-how-to-convert-xml-to-dataframe-in-r-with-the-correct-structure
xml2df <- function(sequences.xml) { #Markes XML structure: Sequence(Sample(Method, Result), File)
  samples.xml <- xml_children(sequences.xml)
  samples <- map_dfr(1:xml_length(sequences.xml), function(sample.index) {
    results <- xml_children(samples.xml[[sample.index]]) %>% 
      map_dfr(~ as.list(xml_attrs(.))) %>% 
      mutate(markes_n = sample.index)
    sample.attr <- as.list(xml_attrs(samples.xml[[sample.index]]))
    bind_cols(sample.attr, results)
  })
  sequence.attr <- as.list(xml_attrs(sequences.xml)) %>% within(rm(name, file))
  bind_cols(sequence.attr, samples)
}

#Copy this file from Markes computer: C:\ProgramData\Markes International\TD\Sequence.prj
#rename with date and add .xml to file name
#this contains a list of all traps ever run in Markes autosampler batches, in XML format
#reading in and converting the large file is slow
filedate <- "230908" #update to YYMMDD date the file was copied off
sequ <- read_xml(paste0("data/Sequence", filedate, ".prj.xml")) %>% 
  xml_find_all("//Sequence") %>% 
  map_df(xml2df) %>% 
  select(id, timestamp, result, type, status, name, value, markes_n) %>% 
  pivot_wider(names_from=name, values_from=value) %>% 
  type.convert(as.is=F) %>% 
  rename_with(make.names) %>%
  rename(Tube = X0.0.0.Markes.TD..SeqVar_Tube,
         RecollectionType = X0.0.0.Markes.TD..SeqVar_RecollectionType) %>% 
  mutate(sequence.start = markestime(timestamp),
         across(c(Desorb.End.Time, Desorb.Start.Time, Trap.Fire.Time), as.POSIXct),
         Desorb.Start.Time = fill_starttime(Desorb.Start.Time)) #fill in missing times with midpoint of surrounding times

sequ$Desorb.Start.Time[3239:3240] <- sequ$Desorb.Start.Time[3238] + c(0.5,1)*60*60 #patch a couple NAs

# Chromatogram data files -------------------------------------------------------------------

#run the following script on GC-MS computer to get a list of the files
#need the file creation time to be accurate so it won't work after copying off the files
#C:/GCMSsolution/Data/Project1/powershell_dir.bat contains the following Windows Powershell script:
#powershell "Get-ChildItem -Recurse c:\GCMSsolution\Data\Project1 | ForEach-Object {$_ | add-member -name "Owner" -membertype noteproperty -value (get-acl $_.fullname).owner -passthru} | Sort-Object fullname | Select FullName,CreationTime,LastWriteTime,Length,Owner | Export-Csv -Force -NoTypeInformation c:\GCMSsolution\Data\Project1\dir.csv"
#Rename the resulting dir.csv with the date to dirYYMMDD.csv

#TODO will this timezone work for everyone?
windowstime <- function(x) as.POSIXct(strptime(x, format = "%m/%d/%Y %I:%M:%S %p"), tz=Sys.timezone())

qgdfiles <-  read_csv(paste0("data/dir",filedate,".csv")) %>% 
  mutate(across(c(CreationTime, LastWriteTime), windowstime),
         FullName = gsub("\\\\", "/", as.character(FullName)),
         ext = file_ext(FullName), 
         FileName = basename(FullName),
         Folder = str_remove(dirname(FullName),"C:/GCMSsolution/Data/Project1/"),
         create_last = LastWriteTime - CreationTime, #length of time the chromatogram was written to
         taildigits = strsplit(FileName, split="_") %>% map_chr(tail, 1),
         GC_n = ifelse(nchar(taildigits) == 6, as.integer(substr(taildigits,1,3)), NA)) %>% #position in the GC run is recorded in the filename
  filter(ext=="qgd") #only want Shimadzu chromatograms

users <- read_tsv("data/GCMS_users.tsv")

# Merge sequence and files ------------------------------------------------

# fuzzy join the two lists by Markes desorb start time and GC-MS file creation time with a certain offset and tolerance
#time difference shift each year, these are the optimum settings:
#up to 2021: +0 min +- 16 min
#2022: +18 min +- 15 min
#2023: +0 min +- 16 min

cd_offset <- 0 # files created at least this many minutes after desorb start time 
cd_tolerance <- 16 # fuzziness before or after (min)

sequ.file <- difference_full_join(sequ, qgdfiles %>% mutate(CreationTime = CreationTime - minutes(cd_offset)),
                                  by=c("Desorb.Start.Time"="CreationTime"), max_dist = cd_tolerance*60, distance_col="create_desorb") %>% 
  add_count(result, name="fuzzy_n") %>%  #result uniquely identifies each row in sequ, so this counts the created duplicates
  #slice_min(create_desorb) %>% # drop all but the best time match
  mutate(markes_GC = markes_n - GC_n, #difference between position in Markes and GC runs - should match unless a skip happened
         CreationTime = CreationTime + minutes(cd_offset), #offset by this many minutes
         CreationDate = date(CreationTime),
         SequenceDate = date(sequence.start),
         eithertime = pmin(CreationTime, Desorb.Start.Time, na.rm=T), #miniumum of the two times, or just one if the other is missing
         either_n = pmin(markes_n, GC_n, na.rm=T), #minimum of the two positions, or just one if the other is missing
         status = factor(paste(Tube.Status,	Unity.Deviation)), #collect Markes error messages from two columns
         project = str_extract(FullName, paste0(users$project, collapse = "|"))) %>% #parts of filename indicate the project
  left_join(users)

sequ.summary <- sequ.file %>% 
  select(c(sequence.start, Desorb.Start.Time, CreationTime, eithertime, status, Tube, markes_n,
           GC_n, either_n, markes_GC, create_desorb,  FileName, user, FullName, id, fuzzy_n)) %>% #columns needed for annotating verdicts
  arrange(sequence.start, markes_n, eithertime) %>%
  mutate(desorb.Start.diff = c(NA, diff(Desorb.Start.Time))) %>% #time differences between desorptions
  write_csv(paste0("output/sequsummary",filedate,".csv"))

save(sequ, qgdfiles, sequ.file, cd_offset, cd_tolerance, 
     file=paste0("output/markes_sequence",filedate,".rda"))
