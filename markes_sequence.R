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
#these settings (+18 min +- 15 min) are optimized for picking out 2022 samples
#but time difference shift each year
#up to 2021 the settings were offset: +0 min tolerance: +- 16 min
cd_offset <- 18 # files created at least this many minutes after desorb start time 
cd_tolerance <- 15 # fuzziness before or after (min)

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

# Plots -----------------------------------------------------------
sequ.file %>% count(fuzzy_n, year(eithertime)) #large count is result=NA from no Markes match, eithertime=NA are skips
#sequ.file %>% filter(fuzzy_n > 1, fuzzy_n < 100) %>% View()
nrow(sequ); nrow(qgdfiles); nrow(sequ.file) #TODO some duplicates created in the merge

plot(y~ x,   data=data.frame(x = tail(sort(qgdfiles$CreationTime - minutes(cd_offset)), 80), y = 1))
points(y~ x, data=data.frame(x = tail(sort(sequ$Desorb.Start.Time), 80), y = 1.02), add=T, col="blue")

tibble(create_desorb = as.numeric(diff(sort(c(qgdfiles$CreationTime, sequ$Desorb.Start.Time))), units="mins")) %>% 
  mutate(n=row_number()) %>% 
  ggplot(aes(x=n, y=create_desorb)) + geom_point(size=0.2) +
  geom_hline(yintercept=c(0, cd_offset), color="blue") +
  geom_hline(yintercept=c(cd_offset-cd_tolerance, cd_offset+cd_tolerance), color="red") +
  scale_y_continuous(limits = c(0,60), n.breaks=30)

sequ.file %>% drop_na(create_desorb) %>% 
  mutate(n=row_number(), create_desorb=as.numeric(create_desorb, units="mins")) %>% 
  ggplot(aes(shape=factor(year(eithertime)), x=n, y=create_desorb, color=fct_lump_n(Folder,10))) +
  geom_point() + labs(y="GC-MS file creation - Markes desorb start time (min)", shape="Year", color="Folder")+
  geom_hline(yintercept=c(0,16), color="blue")
#guessing the second peak in each year is the longer 10deg GC method #2022: not sure if this is the case

sequ.file %>% filter(CreationDate> as.POSIXct("2019-06-19")) %>% 
  select(CreationTime, Desorb.Start.Time) %>% 
  map_dfr(~ as.numeric(diff(.), units="mins")) %>% mutate(n = row_number()) %>% 
  pivot_longer(contains("Time")) %>% 
  ggplot(aes(x=n, y=value, color=name, shape=name)) + 
  geom_point() +  ylim(c(0,90))+ theme_minimal()+ labs(y="Time between samples (min)")

sequ.file %>% arrange(eithertime) %>% mutate(n=row_number()) %>%
  ggplot(aes(x=n, y=markes_GC)) + geom_point() +
  labs(y="Markes batch index - GC batch index")

sequ.file %>% drop_na(eithertime) %>% 
  ggplot(aes(x=SequenceDate, fill=fct_lump_n(Folder,10))) + 
  geom_histogram(binwidth=1) + facet_wrap(vars(year(eithertime)), ncol=1, scales="free_x") + labs(fill="Folder")

sequ.file  %>% drop_na(eithertime) %>% filter(status !="Desorbed NA") %>% 
  ggplot(aes(x=SequenceDate, fill=status)) + 
  geom_histogram(binwidth=1) + facet_wrap(vars(year(eithertime)), ncol=1, scales="free_x")

# Read Shimadzu output ----------------------------------------------------
#adapted from maxfield_gc.R
setwd("../RMBL Batches")
source("../read_shimadzu.R")
#2022 files
datafiles <- list.files(pattern=".txt")
datafiles.22 <- grep("_22", datafiles, value=T)
#safemadzu <- safely(read.shimadzu)
#shimadzu.data.22 <- map(datafiles.22, safemadzu)
#names(shimadzu.data.22) <- datafiles.22
#map(shimadzu.data.22, "error") # deleted partial entry at bottom of "PG_Ihyb_220804_7.txt"
shimadzu.data.22 <- datafiles.22 %>% set_names() %>% map(read.shimadzu) %>% bind_rows(.id="batch")
save(shimadzu.data.22, file = "shimadzu_data_22.rda")

library(reshape2)
all.22 <- dcast(shimadzu.data.22, Filename~Name, sum, value.var="Area")
rownames(all.22) <- all.22[,1]
all.22[,1] <- NULL
all.22.cut <- all.22[,colSums(all.22)>5e8]#arbitrary cutoff

#k-means 
library(vegan)
k <- 40
set.seed(1)
km <- kmeans(decostand(all.22.cut, method="log"), k, nstart=3)

all.km <- tibble(FileName=rownames(all.22)) %>% 
  mutate(rowSum = rowSums(all.22),
         Project = str_extract(FileName, "Blank|[aA]ir") %>% replace_na("sample"),
         Type = fct_collapse(Project, blank="Blank", other_level = "sample"),
         nameBlank = Type=="blank",
         runYear = str_extract(FileName, "2018|2019|2020|2021|2022") %>% replace_na("2018") %>% factor,
         Cluster = km$cluster) %>% # Figure out which k-means clusters are the blanks
  mutate(kBlank = Cluster %in% (count(., nameBlank, Cluster) %>% filter(nameBlank, n>2) %>% pull(Cluster)),
         Mixup = nameBlank != kBlank)

with(all.km, table(kBlank, nameBlank)) #kmeans didn't do well at finding blanks (3:3)

allgc <- sequ.summary %>% #get entire batch if it had a sample that matches
  left_join(all.km %>% select(FileName, nameBlank, Mixup, kBlank, Cluster)) %>% 
  mutate(verdict="", sample="", index=row_number()) %>% 
  left_join(shimadzu.data.22 %>% rename(FileName=Filename) %>% group_by(FileName,batch) %>% tally(name="n_peaks")) %>% #get batch file names from Shimadzu output
  select(c("index", "sequence.start", "batch", "Desorb.Start.Time", "CreationTime", "eithertime", "status", 
           "Tube", "markes_n", "GC_n", "either_n", "markes_GC", "create_desorb", "desorb.Start.diff", 
           "Mixup", "nameBlank", "kBlank", "Cluster", "n_peaks", "verdict", "FileName", "sample", "user", "FullName", "id","fuzzy_n")) %>% 
  filter((!is.na(eithertime) & eithertime > ymd("2022-06-01")) | (!is.na(sequence.start) & sequence.start > ymd("2022-06-01"))) %>% 
  write_csv("../Inventory/2022gc220929.csv")

# Shimadzu batch tables ---------------------------------------------------
setwd("../RMBL Batches")
batchfolder <- "C:\\GCMSsolution\\Data\\Project1\\RMBL Batches\\"
get_batchtxt <- function(vec) {
  batchfolder.regex <- paste0("^", str_replace_all(batchfolder, fixed("\\"),"\\\\"),".*\\.txt")
  str_subset(vec, batchfolder.regex)[1] %>% 
    str_remove(fixed(batchfolder))
}
startline <- "C:\\GCMSsolution\\Data\\Project1\\RMBL Reports" #line before list of sample names
endline <- "Markes test 1.qgm" #line after inaccurate tabbed list of sample names
get_filelines <- function(vec, start, end) {
  vec[(start+3):(end-3)] %>% paste(collapse =";") %>% 
    str_split(fixed(";Unknown;")) %>% unlist()
}
batchtables <- list.files("../RMBL Batches", pattern = ".qgb") %>% 
  file.info() %>% select(size, mtime) %>% rownames_to_column("qgbfile") %>%  
  mutate(strings = map(qgbfile, ~ system(paste0("strings \"", .x, "\""), intern=TRUE)), 
         batch = map_chr(strings, get_batchtxt),
         n_lines = map_int(strings, length),
         n_sections = map_int(strings, ~sum(.x == startline)),
         start_pos = map_int(strings, ~match(startline, .x)),
         end_pos =   map_int(strings, ~match(endline, .x)),
         between = end_pos - start_pos,
         size = size/(2^12)) %>% 
  filter(between > 5) %>% 
  mutate(file_lines = pmap(list(strings, start_pos, end_pos), get_filelines)) %>% select(-strings) %>% 
  unnest(file_lines) %>% group_by(qgbfile) %>% mutate(seq_n = row_number()) %>% ungroup() %>% 
  arrange(mtime, seq_n) %>% 
  write_csv("../Inventory/qgb_batches.csv", quote="all")
