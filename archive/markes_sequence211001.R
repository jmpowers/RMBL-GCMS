library(xml2)
library(bit64)
library(tools)
library(fuzzyjoin)
library(tidyverse)
library(lubridate)
library(googlesheets4)

# Markes sequence ---------------------------------------------------------
setwd("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/Inventory")
seq.xml <- read_xml("Sequence210826.prj.xml")

markestime  <- function(x) as.POSIXct(as.POSIXlt(as.integer64(x)/10000000, origin = "0001-01-01",tz="UTC"))

#From https://stackoverflow.com/questions/34273132/r-how-to-convert-xml-to-dataframe-in-r-with-the-correct-structure
sequ <- bind_rows(lapply(xml_find_all(seq.xml, "//Sequence"), function(sequence) {
  parent <- data.frame(as.list(xml_attrs(sequence)), stringsAsFactors=FALSE)
  kids.xml <- xml_children(sequence)
  kids <- bind_rows(lapply(1:xml_length(sequence), function(sample.index) {
    grandkids <- bind_rows(lapply(xml_children(kids.xml[[sample.index]]), function(result) {
      as.list(xml_attrs(result))
    }))
    grandkids$markes_n = rep(sample.index, nrow(grandkids))
    cbind.data.frame(as.list(xml_attrs(kids.xml[[sample.index]])), grandkids, stringsAsFactors=FALSE)
  }))
  cbind.data.frame(parent, kids, stringsAsFactors=FALSE)
})) %>% select(id, timestamp, result, type, status, name=name...11, value, markes_n) %>% 
  spread(name, value) %>% 
  type.convert %>% 
  mutate(sequence.start = markestime(timestamp)) %>%
  mutate_at(vars(`Desorb End Time`, `Desorb Start Time`, `Trap Fire Time`), as.POSIXct) %>%
  rename_all(make.names)  %>%
  rename(Tube = X0.0.0.Markes.TD..SeqVar_Tube) %>%
  rename(RecollectionType = X0.0.0.Markes.TD..SeqVar_RecollectionType)
ds <- sequ$Desorb.Start.Time
ds[is.na(ds)] <- as.POSIXct(unlist(apply(sapply(which(is.na(ds)), "+", c(-1, 1)), 2,  function(x) ifelse(sum(x)>2 && as.double(diff(ds[x]), units = "mins")>50 , mean(ds[x]), NA))), origin = lubridate::origin)
sequ$Desorb.Start.Time <- ds
sequ$Desorb.Start.Time[3239:3240]<-sequ$Desorb.Start.Time[3238]+c(0.5,1)*60*60

# Data Files -------------------------------------------------------------------
#run in Windows Powershell on GC-MS machine (need file CreationTime to be accurate)
#powershell "Get-ChildItem -Recurse c:\GCMSsolution\Data\Project1 | ForEach-Object {$_ | add-member -name "Owner" -membertype noteproperty -value (get-acl $_.fullname).owner -passthru} | Sort-Object fullname | Select FullName,CreationTime,LastWriteTime,Length,Owner | Export-Csv -Force -NoTypeInformation c:\GCMSsolution\Data\Project1\dir190707.csv"

windowstime <- function(x) as.POSIXct(strptime(x, format = "%m/%d/%Y %I:%M:%S %p"),tz=Sys.timezone())

qgdfiles <-  read_csv("dir210827.csv") %>% 
  mutate_at(vars(CreationTime, LastWriteTime), windowstime) %>% 
  mutate(FullName = gsub("\\\\", "/", as.character(FullName)),
         ext = file_ext(FullName), 
         FileName = basename(FullName),
         create_last = LastWriteTime-CreationTime,
         taildigits <- strsplit(FileName, split="_") %>% map_chr(tail, 1),
         GC_n = ifelse(nchar(taildigits) == 6, as.integer(substr(taildigits,1,3)), NA))  %>% 
  filter(ext=="qgd")

# Merge sequence and files ------------------------------------------------
sequ.file <- difference_full_join(sequ, qgdfiles, by=c("Desorb.Start.Time"="CreationTime"), max_dist = 16*60, distance_col="create_desorb") %>% #16 min
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

ggplot(sequ.file %>% drop_na(eithertime), aes(x=create_desorb, fill=factor(year(eithertime)))) + 
  geom_histogram(binwidth=30) + facet_wrap(vars(year(eithertime)), ncol=1)
#guessing the second peak in each year is the longer 10deg GC method

sequ.summary <- sequ.file %>% select(c(sequence.start, Desorb.Start.Time, CreationTime, eithertime, status, Tube, markes_n, GC_n, either_n, markes_GC, create_desorb,  FileName, user, FullName, id)) %>% 
  arrange(sequence.start, markes_n, eithertime) %>%
  mutate(desorb.Start.diff = c(NA, diff(Desorb.Start.Time))) 
write_csv(sequ.summary, "sequfile.csv")
  
# Schiedea runs -----------------------------------------------------------
Powers <- grepl("Powers RMBL Data", sequ.summary$FullName, fixed=T)
schiedea.ids <- sequ.summary %>% filter(Powers) %>% select(id) %>% unique %>% na.omit
schgc <- sequ.summary %>% filter(id %in% schiedea.ids$id | Powers)
schgc$vial <- as.integer(str_match(schgc$FileName, "_V(\\d{1,3})")[,2])
schgc$dupe <- !is.na(schgc$vial) & schgc$vial %in% schgc$vial[duplicated(schgc$vial)]
#Run  kmeans first
schgc <- cbind(schgc, sch.km[,-1])
rownames(schgc) <- 1:nrow(schgc)
setwd("../Schiedea")
notebook <- read.delim("leak_notebook.csv")
schgc <- cbind(schgc, notebook[match(schgc$vial, as.integer(as.character(notebook$Vial)), incomparables=NA),-4])
schgc$TubeMatch <- schgc$Tube == schgc$SteelTube
write.csv(schgc, "schiedea_all190815.csv")
#read back the verdicts
rmbl <- gs_title("Schiedea RMBL GC-MS Inventory")
schgc.verdict <- gs_read(rmbl)
schgc$verdict <- schgc.verdict$verdict
schgc$FileName2 <- ifelse(is.na(schgc.verdict$sample), schgc$FileName, schgc.verdict$sample)
schgc$vial2 <- as.integer(str_match(schgc$FileName2, "_V(\\d{1,3})")[,2])
schgc$dupe2 <- !is.na(schgc$vial2) & schgc$vial2 %in% schgc$vial2[duplicated(schgc$vial2)]
#schgc.use <- schgc[!(schgc$verdict %in% c("alreadyrun", "empty","leak-blank", "notmine", "notrun", "skip-blank", "skip-notrun")),]
#match up with metadata
sampl <-gs_title("Schiedea Volatiles Sampling")
s <- gs_read(sampl, col_types = cols(.default = col_character(),
                                     Flow = col_double(),
                                     Mass = col_double(),
                                     RunDate = col_date(format = ""),
                                     SampleDate = col_date(format = ""),
                                     Start = col_time(format = ""),
                                     Stop = col_time(format = ""),
                                     Equi = col_double(),
                                     Duration = col_double(),
                                     Total = col_double()
))
s.nr <- s[s$GC=="NR",]
s.nr$Vial <- as.integer(s.nr$Vial)
sort(union(s.nr$Vial, schgc$vial2))
sort(setdiff(s.nr$Vial, schgc$vial2)) #not run yet
sort(setdiff(schgc$vial2, s.nr$Vial)) #run but not entered  
broken <- c(52,53,56,86,117,124,125,129,174) #broken traps - not run
empty <- c(7,16,18,28,40,313,314,315,316) # presumed clean - run
notrun <- c(5,8,21,30,33,50,51,55,57,64,66,74,75,80,81,82,83,84,88,89,91,95,100,102,103,104,107,108,109,11,12,22,23,24,25,26,58,59,60,61,62,63,58,92,94,96,97,98,9, 90, 79, 87, 72, 99)
setdiff(setdiff(setdiff(s.nr$Vial, schgc$vial2), broken), notrun)
s.nr$DN <- factor(ifelse(s.nr$Start > 16*60*60, "Night", "Day"))
#s.nr$FileName <- schgc$FileName2[match(s.nr$Vial, schgc$vial2)]
sv.all <- merge(s.nr, schgc, by.x = "Vial", by.y="vial2", all.x = T, all.y = T)
write.csv(sv.all, "s.nr.schgc190815.csv")

# Maxfield runs -----------------------------------------------------------
setwd("../Maxfield/")
Maxfield <- grepl("Maxfield|07262018|June28_07292018|Ipomopsis 2020|Corydalis|8112021|8122021|8132021", sequ.summary$FullName)
#Try to get all 2021 data to annotate
#811-8132021 are for picking out 2021 Lupinus argenteus samples
#There are some 2021 GNA, VF runs in Test instrument runs too, but adding that folder would break pre-2021 verdict order
maxfield.ids <- sequ.summary %>% filter(Maxfield) %>% select(id) %>% unique %>% na.omit
maxfgc <- sequ.summary %>% filter(id %in% maxfield.ids$id | Maxfield)

#now get kmeans results from maxfield_gc.R
maxfgc <- maxfgc %>% cbind(maxf.km[,-1]) %>% 
  mutate(verdict="", sample="", index=row_number()) %>% 
  left_join(maxf.data %>% rename(FileName=Filename) %>% group_by(FileName,batch) %>% tally(name="n_peaks")) %>% #get batch file names from Shimadzu output
  select(c("index", "sequence.start", "batch", "Desorb.Start.Time", "CreationTime", "eithertime", "status", 
           "Tube", "markes_n", "GC_n", "either_n", "markes_GC", "create_desorb", "desorb.Start.diff", 
           "Mixup", "nameBlank", "kBlank", "Cluster", "n_peaks", "verdict", "FileName", "sample", "user", "FullName", "id"))
write_csv(maxfgc, "maxfield_all210827.csv")

#Output only Maxfield samples to split filenames into parts
data_inventory <- "1X8oo7qZlo1p6MVl_CBeBe6CUTHEAcd-FWQzfHud3Qws"
read_sheet(data_inventory, sheet="maxfield_all210827annot", na="NA") %>% 
  filter(user=="D Campbell : Maxfield") %>% 
  mutate(sample = na_if(sample,"") %>% coalesce(FileName) %>% str_remove(".qgd")) %>% 
  select(index, batch, Desorb.Start.Time, verdict, sample) %>% 
  separate(sample, into=paste0("file",1:9), remove=F) #%>% 
  #range_write(ss=data_inventory, sheet="maxfield_meta",range="A:N")

#get hand-split filename metadata
maxfmeta <- read_sheet(data_inventory, sheet="maxfield_meta", guess_max=2000, col_types="c") %>% 
                                   select(index, type:vial) %>% as.data.frame %>% 
  drop_na(type) %>% #ditch the gunk for now
  distinct(index, .keep_all = T) %>% # TODO investigate these dupes with two Shimadzu batches
  left_join(maxfgc %>% distinct(index, .keep_all = T) %>% mutate(index=as.character(index))) %>%
  mutate(plotid=ifelse(type=="floral", str_sub(plantid, 1,2),NA),
                  plant= ifelse(type=="floral", str_sub(plantid, 3),NA),
                  sampledate=ymd(ifelse(year(sequence.start)==2018, paste0("2018",sampledate), sampledate)))

# Plots and tables-----------------------------------------------------------
t(with(sequ.file, table(user, CreationDate)))
t(with(sequ.file, table(status, SequenceDate)))

plot(diff((sequ.file%>% filter(CreationDate> as.POSIXct("2019-06-19")))$CreationTime)/60,ylim=c(0,90), col="blue", pch="|")
points(diff((sequ.file%>% filter(CreationDate> as.POSIXct("2019-06-19")))$Desorb.Start.Time)/60, col="red", pch=".", cex=3)

ggplot(sequ.file, aes(SequenceDate, status)) + geom_point(alpha=0.8, shape="|") + theme_minimal() #misleading because depends on use

plot(c(tail(sort(qgdfiles$CreationTime), 60),tail(sequ$Desorb.Start.Time, 60)))
plot(diff(sort(c(qgdfiles$CreationTime, sequ$Desorb.Start.Time)))/60, ylim=c(0,60), col=rep(2:3), pch=".")
abline(16,0)
abline(0,0)

hist(as.numeric(sequ.file$create_desorb)/60)
table(sequ.file$markes_GC)


