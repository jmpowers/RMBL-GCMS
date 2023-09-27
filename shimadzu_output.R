library(tidyverse)
library(reshape2)
library(vegan)

# Read Shimadzu output ----------------------------------------------------
#adapted from maxfield_gc.R
source("read_shimadzu.R")
#TODO had to chop off VM_Linaria_230903.txt, some files were not NIST searched - see "2023/olderversion" subfolder
thisyear <- 2023 
setwd(paste0("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/RMBL Batches/", thisyear))
datafiles <- list.files(pattern=".txt")
safemadzu <- safely(read.shimadzu)
#shimadzu.data <- map(datafiles, safemadzu)
#names(shimadzu.data) <- datafiles
#map(shimadzu.data, "error") 
shimadzu.data <- datafiles %>% set_names() %>% map(read.shimadzu) %>% bind_rows(.id="batch")
setwd("~/MyDocs/MEGA/UCI/RMBL 2020/RMBL-GCMS/")
save(shimadzu.data, file = paste0("data/shimadzu_data_",thisyear,".rda"))

all <- dcast(shimadzu.data, Filename~Name, sum, value.var="Area")
rownames(all) <- all[,1]
all[,1] <- NULL
all.cut <- all[,colSums(all)>5e8]#arbitrary cutoff

#k-means 
k <- 40
set.seed(1)
km <- kmeans(decostand(all.cut, method="log"), k, nstart=3)

all.km <- tibble(FileName=rownames(all)) %>% 
  mutate(rowSum = rowSums(all),
         Project = str_extract(FileName, "Blank|[aA]ir") %>% replace_na("sample"),
         Type = fct_collapse(Project, blank="Blank", other_level = "sample"),
         nameBlank = Type=="blank",
         runYear = str_extract(FileName, paste0(2018:thisyear, collapse="|")) %>% replace_na("2018") %>% factor(),
         Cluster = km$cluster) %>% # Figure out which k-means clusters are the blanks
  mutate(kBlank = Cluster %in% (count(., nameBlank, Cluster) %>% filter(nameBlank, n>2) %>% pull(Cluster)),
         Mixup = nameBlank != kBlank)

with(all.km, table(kBlank, nameBlank)) #lots of false positives

filedate <- "230908"
load(paste0("output/markes_sequence",filedate,".rda"))
allgc <- sequ.file %>% 
  arrange(sequence.start, markes_n, eithertime) %>%
  mutate(desorb.Start.diff = c(NA, diff(Desorb.Start.Time))) %>% 
  left_join(all.km %>% select(FileName, nameBlank, Mixup, kBlank, Cluster)) %>% 
  mutate(verdict="", sample="", index=row_number()) %>% 
  left_join(shimadzu.data %>% rename(FileName=Filename) %>% group_by(FileName,batch) %>% tally(name="n_peaks")) %>% #get batch file names from Shimadzu output
  select(c("index", "sequence.start", "batch", "Desorb.Start.Time", "CreationTime", "eithertime", "status", 
           "Tube", "markes_n", "GC_n", "either_n", "markes_GC", "create_desorb", "desorb.Start.diff", 
           "Mixup", "nameBlank", "kBlank", "Cluster", "n_peaks", "verdict", "FileName", "sample", "user", "FullName", "id","fuzzy_n")) %>% 
  filter((!is.na(eithertime) & eithertime > ymd(paste0(thisyear,"-01-01"))) | (!is.na(sequence.start) & sequence.start > ymd(paste0(thisyear,"-01-01")))) %>% 
  write_csv(paste0("output/gc",filedate,".csv"))

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

