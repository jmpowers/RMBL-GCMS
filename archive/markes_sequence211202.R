library(xml2)
library(bit64)
library(tools)
library(fuzzyjoin)
library(tidyverse)
library(lubridate)
library(googlesheets4)

# Markes sequence ---------------------------------------------------------

markestime  <- function(x) as.POSIXct(as.POSIXlt(as.integer64(x)/10000000, origin = "0001-01-01",tz="UTC"))

fill_starttime <- function(ds) {
  ds[is.na(ds)] <- sapply(which(is.na(ds)), "+", c(-1, 1)) %>% #get matrix of indices of preceding and following entries
    apply(2, function(x) ifelse(sum(x)>2 && as.double(diff(ds[x]), units = "mins")>50, mean(ds[x]), NA)) %>% #if 50 min apart, take midpoint
    unlist() %>% as.POSIXct(origin = lubridate::origin)
  ds
}

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

# Sequence.prj file contains listing of all traps ever run in autosampler batches in XML format
# Get it off Markes laptop with a SD card
setwd("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/Inventory")
sequ <- read_xml("Sequence210826.prj.xml") %>% 
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
         Desorb.Start.Time = fill_starttime(Desorb.Start.Time))

sequ$Desorb.Start.Time[3239:3240] <- sequ$Desorb.Start.Time[3238] + c(0.5,1)*60*60 #patch a couple NAs

# Data Files -------------------------------------------------------------------
#run in Windows Powershell on GC-MS machine (need file CreationTime to be accurate)
#powershell "Get-ChildItem -Recurse c:\GCMSsolution\Data\Project1 | ForEach-Object {$_ | add-member -name "Owner" -membertype noteproperty -value (get-acl $_.fullname).owner -passthru} | Sort-Object fullname | Select FullName,CreationTime,LastWriteTime,Length,Owner | Export-Csv -Force -NoTypeInformation c:\GCMSsolution\Data\Project1\dir190707.csv"

windowstime <- function(x) as.POSIXct(strptime(x, format = "%m/%d/%Y %I:%M:%S %p"), tz=Sys.timezone())

qgdfiles <-  read_csv("dir210827.csv") %>% 
  mutate(across(c(CreationTime, LastWriteTime), windowstime),
         FullName = gsub("\\\\", "/", as.character(FullName)),
         ext = file_ext(FullName), 
         FileName = basename(FullName),
         Folder = str_remove(dirname(FullName),"C:/GCMSsolution/Data/Project1/"),
         create_last = LastWriteTime - CreationTime,
         taildigits = strsplit(FileName, split="_") %>% map_chr(tail, 1),
         GC_n = ifelse(nchar(taildigits) == 6, as.integer(substr(taildigits,1,3)), NA))  %>% 
  filter(ext=="qgd")

# Merge sequence and files ------------------------------------------------
# fuzzy join them by Markes start time and GC-MS file creation time (some minutes later) with tolerance of 16 min
sequ.file <- difference_full_join(sequ, qgdfiles, by=c("Desorb.Start.Time"="CreationTime"), max_dist = 16*60, distance_col="create_desorb") %>% 
  mutate(markes_GC = markes_n - GC_n,
         CreationDate = date(CreationTime),
         SequenceDate = date(sequence.start),
         eithertime = pmin(CreationTime, Desorb.Start.Time, na.rm=T),
         either_n = pmin(markes_n, GC_n, na.rm=T),
         status = factor(paste(Tube.Status,	Unity.Deviation)),
         user = 
          ifelse(grepl("SCH_|GLOB_|KAHO_|STEL_", FileName), "John Powers : Sakai/Weller", 
          ifelse(grepl("Paul", FileName, fixed=T), "Lucas Sharret : Diane Campbell",
          ifelse(grepl("Carvajal", FileName, fixed=T), "Nayeli Carvajal : Kailen Mooney",
          ifelse(grepl("Blank|bake", FileName, ignore.case=T), "Blank : RMBL", 
          ifelse(grepl("STD|Alk|indole", FileName, ignore.case=T), "Standard : RMBL",NA))))))

sequ.summary <- sequ.file %>% select(c(sequence.start, Desorb.Start.Time, CreationTime, eithertime, status, Tube, markes_n, GC_n, either_n, markes_GC, create_desorb,  FileName, user, FullName, id)) %>% 
  arrange(sequence.start, markes_n, eithertime) %>%
  mutate(desorb.Start.diff = c(NA, diff(Desorb.Start.Time))) 
write_csv(sequ.summary, "sequfile.csv")
save.image("markes_sequence.rda")

# Plots -----------------------------------------------------------
nrow(sequ); nrow(qgdfiles); nrow(sequ.file) #TODO some duplicates created in the merge

tibble(create_desorb = as.numeric(diff(sort(c(qgdfiles$CreationTime, sequ$Desorb.Start.Time))), units="mins")) %>% 
  mutate(n=row_number()) %>% 
  ggplot(aes(x=n, y=create_desorb)) + geom_point(size=0.2) + ylim(c(0,60)) +
  geom_hline(yintercept=c(0,16), color="blue")

sequ.file %>% drop_na(create_desorb) %>% 
  mutate(n=row_number(), create_desorb=as.numeric(create_desorb, units="mins")) %>% 
  ggplot(aes(shape=factor(year(eithertime)), x=n, y=create_desorb, color=fct_lump_n(Folder,10))) +
  geom_point() + labs(y="GC-MS file creation - Markes desorb start time (min)", shape="Year", color="Folder")
#guessing the second peak in each year is the longer 10deg GC method

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