library(tidyverse)
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
