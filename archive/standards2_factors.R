library(reshape2)
library(tidyverse)
library(lubridate)
library(vegan)
library(viridis)

setwd("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/Standards")
source("../read_shimadzu.R")

datafiles <- list.files("./data", full.names=TRUE) #grab all the files in the data folder - add new Shimadzu output here and give it 2 digit prefix number
stds.data <- do.call(rbind, lapply(datafiles, read.shimadzu)) 
stds.data <- stds.data[stds.data$Filename!="Blank_07192019.qgd",]

stds.all <- dcast(stds.data, Filename~Name, sum, value.var="Area") #long to wide
rownames(stds.all) <- stds.all[,1]
stds.all[,1] <- NULL

mix <- c("3-Hexen-1-ol, (E)-", "3-Hexen-1-ol, (Z)-", "3-Hexen-1-ol", "4-Hexen-1-ol, (E)-", "4-Hexen-1-ol, (Z)-",  "4-Hexen-1-ol",
         "(1R)-2,6,6-Trimethylbicyclo[3.1.1]hept-2-ene", "(1S)-2,6,6-Trimethylbicyclo[3.1.1]hept-2-ene", ".alpha.-Pinene",
         "Linalool", "Methyl salicylate", "Caryophyllene", ".alpha.-Farnesene", "Benzaldehyde", "Indole","2,6,10-Dodecatrien-1-ol, 3,7,11-trimethyl-")
stds <- select(stds.all, mix) %>% 
  mutate(Hexenol = `3-Hexen-1-ol, (E)-`+`3-Hexen-1-ol, (Z)-`+`3-Hexen-1-ol`+`4-Hexen-1-ol, (E)-`+`4-Hexen-1-ol, (Z)-`+`4-Hexen-1-ol`,
         `alpha-Pinene`=`(1R)-2,6,6-Trimethylbicyclo[3.1.1]hept-2-ene` + `(1S)-2,6,6-Trimethylbicyclo[3.1.1]hept-2-ene` + `.alpha.-Pinene`,
         .keep = "unused") %>%   relocate(`alpha-Pinene`) %>% 
  rename(`alpha-Farnesene`=`.alpha.-Farnesene`,`Farnesol`=`2,6,10-Dodecatrien-1-ol, 3,7,11-trimethyl-`)

rownames(stds) <- rownames(stds.all)
mergedmix <- colnames(stds)
write.csv(stds, "./metadata/filelist.csv") #write out all the filenames and areas of the standards
write.table(rownames(stds.all), "./metadata/file_list.csv", quote=FALSE, row.names=FALSE, col.names="Filename") #just the filenames

###then put the new filelist entries into this google doc:
# https://docs.google.com/spreadsheets/d/1eJGRwjfDWT5znfkJADR-wS1tVUWV2bHC/edit#gid=2061496758
# and download as shimadzu_standards.csv

stds.meta <- read.csv("./metadata/shimadzu_standards.csv")[,1:8] %>% #1:8 to get jsut the metadata columns
  mutate(DateBatch = factor(paste(RunDate, Batch)), 
         RunDateP = as.POSIXct(RunDate), 
         Year = factor(paste0(year(RunDateP), ifelse(year(RunDateP)==2018 & RunDateP < as.POSIXct("2018-07-12", origin=lubridate::origin), " leak"," "), ifelse(RunDate=="2021-07-23"," fil2",""))))
stds <- stds[match(stds.meta$FileName, rownames(stds)),] #reorder the metdata to match order of the peak areas. 
meta.mix <- cbind(stds.meta, stds) #joins the metadata and the peak areas
#write.csv(meta.mix, "./metadata/shimadzu_standards.csv", row.names=F)

#stds.rel <- stds / stds["STD3uL0.1rep_07242017.qgd",][col(stds)] / meta.mix$Conc /10
stds.rel <- stds / (stds.meta$Conc * 10) #area relative to the 0.1 mg/mL standards
meta.mix.rel <- cbind(stds.meta, stds.rel)

#skip all the quant stuff for now, go to the Correlations section

######Quantitative integrations########
quantmix <- c("alpha-Pinene"="(1R)-2,6,6-Trimethylbicyclo[3.1.1]hept-2-ene", "Linalool"="Linalool", "Caryophyllene"="Caryophyllene", 
               "alpha-Farnesene"=".alpha.-Farnesene", "Indole"="Indole", "Methyl salicylate"="Methyl salicylate", "Hexenol"="4-Hexen-1-ol, (Z)-", "Benzaldehyde"="Benzaldehyde")#"Farnesol"="1,6,10-Dodecatrien-3-ol, 3,7,11-trimethyl-"
quantmix.short <- set_names(names(quantmix), c("a-pinene","linalool","b-caryophyllene","a-farnesene","indole","methyl 2-hydroxybenzoate","(E)-hex-3-en-1-ol", "benzaldehyde"))#"farnesol"
quantpath <- "~/MyDocs/MEGA/UCI/RMBL 2020/ipomopsis-temp/data/volatiles/"
quant_slopes <- read_tsv(paste0(quantpath, "quant_slopes.tsv")) %>% left_join(enframe(quantmix, name="std", value="Name")) %>% 
  drop_na(std) %>% select(std, value) %>% deframe()
  
quant.long <- read_tsv(paste0(quantpath, "quant.tsv")) %>% 
  filter(Filename %in% stds.meta$FileName | str_detect(Filename, "[Bb]lank"), Name %in% quantmix) %>% 
  rename(FileName=Filename)
quant.stds <- quant.long %>% pivot_wider(names_from = "Name", values_from="Area", values_fn=mean) %>% #a couple duplicate blanks 
   rename(!!!quantmix)
meta.quant <- left_join(stds.meta %>% filter(Type=="Standards", !is.na(Conc), !Conc %in% c(0.0166, 0.0333)),
                        quant.stds) %>% drop_na(batch) #need to integrate a few more

blank.quant <- quant.stds %>% filter(str_detect(FileName, "^Blank")) %>% 
  mutate(Conc=0,  Year=str_match(FileName,"(20[1-2][0-9])[._]")[,2], 
         daymonth= str_match(str_remove_all(FileName, "-"), paste0("([0-9]{2,4})", Year))[,2] %>%
           str_remove("^0"), 
         RunDate = ymd(paste(Year, str_sub(daymonth,1,1), str_sub(daymonth, 2, -1)))) %>% 
  drop_na(Year) %>% mutate(Year = if_else(Year == "2021" & RunDate >= ymd("2021-07-22"), "2021  fil2", Year))

batch_slopes <- meta.quant %>% pivot_longer(all_of(names(quantmix))) %>% 
  group_by(name, DateBatch) %>% nest() %>% 
  mutate(n=map_int(data, nrow), model = map(data, ~lm(sqrt(value) ~ sqrt(Conc), data=.x)), slope = map_dbl(model, ~coef(.x)[2])) %>% select(-data, -model) %>% drop_na(slope) %>% 
  group_by(name) %>% mutate(slope = slope/max(slope))
  #drop_na(slope) %>% filter(slope < 5000) %>% pull(DateBatch)

#heatmap of slopes
ggplot(batch_slopes, aes(x=name, y=DateBatch, fill=slope)) + geom_tile() + scale_fill_viridis(option="inferno")

#heatmap of amounts
meta.quant %>% 
  pivot_longer(all_of(names(quantmix))) %>% 
  group_by(Conc, name) %>% 
  mutate(value = value/max(value)) %>% 
  ggplot(aes(x=name, y=paste(DateBatch, Sort), fill=value)) + geom_tile() + scale_fill_viridis(option="inferno") +
  facet_wrap(vars(Conc), scales="free_y")

#look at heatmaps to exclude these dud batches:
batches_sans_stds <- c("2018 leak", "newxfer_leak",  "Diane", "2019fresh","2018old", "Warmup")
batches_sans_indole <- c(batches_sans_stds, "freshstds_noindole")
datebatches_sans_stds <- c("2018-09-03 2018 no leak")
datebatches_sans_indole <- c(datebatches_sans_stds, "2019-06-21 2018fresh", "2017-06-30 2017 good")

####Compare quant and qual####
stds.quantqual <- left_join(meta.quant %>% pivot_longer(all_of(names(quantmix)), values_to="quant"),
                            meta.mix %>% pivot_longer(all_of(mergedmix), values_to="qual")) %>% 
  mutate(combined = if_else(qual==0, quant, qual*quant_slopes[name])) #leave combined unscaled

stds.quantqual %>% group_by(name, Year, Conc) %>% summarize(across(c("quant", "qual", "combined"), ~mean(.x>0))) %>% 
  pivot_longer(c("quant", "qual", "combined"), names_to="measure") %>% 
  ggplot(aes(x=factor(Conc), y=value, color=name, shape=measure, linetype=measure)) + facet_wrap(vars(Year)) + 
  geom_point() + geom_path(aes(group=paste(name, measure)))+ theme_minimal()

ggplot(stds.quantqual, aes(x=qual, y=quant)) + facet_wrap(vars(name)) + geom_point(aes(color=Year))+
  scale_x_sqrt()+ scale_y_sqrt() + scale_color_brewer(palette = "Set1") + theme_minimal()

#replace standard's unscaled quant values with combined unscaled quant + qual
meta.quant2 <- stds.quantqual %>% select(-quant, -qual) %>% pivot_wider(values_from="combined")

#TODO the quant regressions aren't always working in 2019
# for the high values of
# indole: peak is shifted late (13.25 vs. 12.96)
# caryophyllene: double peak - first is small at 14.85 (target 14.8), second is huge at 15.10
# linalool: overloaded peak, arrives late

####ambient/floral metadata####
#grab from quant from ipotemp.Rmd
#combine quant values with combined quant and qual
floral.quant.wide <- bind_rows(quant =   quant[as.character(metadata$sample),names(quantmix.short)] %>%
                                 rownames_to_column("Filename"),
                               qual =  ipo.all[as.character(metadata$sample),names(quantmix.short)] %>% 
                                 rownames_to_column("Filename"),
                               .id="measure") %>% 
  rename_with(~quantmix.short[.x], .cols=all_of(names(quantmix.short))) %>% 
  pivot_longer(all_of(names(quantmix))) %>% 
  pivot_wider(names_from="measure") %>% 
  mutate(combined = if_else(qual==0, quant, qual*quant_slopes[name])) %>% #leave combined unscaled
  select(-quant, -qual) %>% pivot_wider(values_from = "combined")

floral.quant <- left_join(metadata %>% select(Year, Folder, type, Expt, PlantType, DN, RunDate, Filename=sample),
                          floral.quant.wide) %>%
  filter(type %in% c("floral","ambient")) %>% mutate(Year=as.character(Year))

#combine standards, blanks, and floral/ambient samples
mega.quant <- bind_rows(meta.quant2 %>% mutate(RunDate=ymd(RunDate)), blank.quant, floral.quant) %>% 
  mutate(ConcType = coalesce(as.character(Conc), as.character(type)) %>% factor() %>% 
           fct_relevel(c("0","ambient","floral","1e-04"))) %>% 
  drop_na(ConcType, Year) %>% mutate(Year=trimws(Year)) %>% 
  filter(!Year %in% c("2020","2018 leak"),
         !(Batch %in% batches_sans_indole | DateBatch %in% datebatches_sans_indole))

with(mega.quant, table(Year, ConcType))

#### Quant regressions ####
lapply(names(quantmix), function(chem) {
  ggplot(mega.quant %>% rename(all_of(c(Chem=chem))), aes(x=Conc, group=DateBatch, color=Year, y=Chem/quant_slopes[chem])) +
    scale_x_log10("Concentration (mg/mL)") + scale_y_log10("Area") +
    stat_smooth(method="lm", se=F, geom="line", alpha=0.6, size=1) + 
    geom_point(pch=24, size=2) + 
    scale_color_brewer("Year", palette="Set1", type="qual") + 
    theme_minimal() + ggtitle(chem)
  ggsave(filename=paste0("./regressions-quant/regression",chem,".svg"), width=8, height=6)
})


#### Quant boxplots ####

#limit of blank, see Clin Biochem Rev. 2008 Aug; 29(Suppl 1): S49–S52
#one-sided test at 5% level = 1.645 Function: ~mean(.x) + qnorm(.95)*sd(.x)
#but data is not normally distributed, use quantile instead
LoB <- blank.quant %>% group_by(Year) %>% summarize(across(all_of(names(quantmix)), ~quantile(.x, 0.95)))
ggplot(blank.quant, aes(x=Indole, fill=Year)) + geom_density()+ scale_x_sqrt()+facet_wrap(vars(Year), ncol=1) +
  geom_point(aes(y=0, color=Year))+
  geom_vline(data=LoB, aes(xintercept=Indole)) +
  guides(fill="none", color="none")
#limit of ambient
LoA <- floral.quant %>% filter(type=="ambient") %>%
  group_by(Year) %>% summarize(across(all_of(names(quantmix)), ~quantile(.x,0.95)))
 
lapply(names(quantmix), function(chem) {
    ggplot(mega.quant %>% rename(all_of(c(Chem=chem))),
           aes(x=ConcType, color=factor(ConcType), y=Chem/quant_slopes[chem])) + 
    facet_wrap(vars(Year), nrow=1)+ geom_boxplot(size=1)+ 
    geom_hline(data = LoB %>% rename(all_of(c(Chem=chem))), aes(yintercept=Chem/quant_slopes[chem], color="0"),  linetype=2) + 
    geom_hline(data = LoA %>% rename(all_of(c(Chem=chem))), aes(yintercept=Chem/quant_slopes[chem], color="ambient"),  linetype=2) +
    geom_text(data = mega.quant %>% group_by(ConcType, Year)  %>% rename(all_of(c(Chem=chem))) %>% 
                summarize(nonzero=mean(Chem>0)), 
              aes(y=250*nonzero, label=round(nonzero*100)), show.legend=F)+
    scale_y_log10(n.breaks=6, minor_breaks=NULL) + scale_color_brewer(palette="Set2")+
    labs(color="Concentration (mg/mL)", y="Quantitative integration", title=paste(chem, "in blanks and standards"), 
         subtitle="Numbers show percent detection, nondetects excluded from boxplots") +
    theme_bw() + theme(axis.text.x=element_blank(), panel.grid.major.x = element_blank(), 
                       axis.ticks=element_blank(), axis.title.x=element_blank(), 
                       legend.position = "top") 
  ggsave(filename=paste0("./boxplots-quant/boxplot",chem,".svg"), width=12, height=7)
})


####Correlations######
plot(stds) #generally points should not be on the edges (=nondetect)
barplot(t(stds[order(rowSums(stds)),]), col=rainbow(ncol(stds)), las=2, cex.names=0.3, names.arg=unite(stds.meta[order(rowSums(stds)),],TypeConc, Type, Conc, RunDate)$TypeConc) #eye candy

####NMDS######
nmds.stds <- metaMDS(stds[rowSums(stds)>0,], dist="bray")
ordiplot(nmds.stds, display="sites", type = "t") #not very helpful


####Regression######
lapply(mergedmix, function(chem) {
  ggplot(meta.mix, aes_string(x="Conc", group="DateBatch", color="Year", y=as.name(chem))) +
    scale_x_log10("Concentration (mg/mL)") + 
    scale_y_log10("Area") +
    stat_smooth(method="lm", se=F, geom="line", alpha=0.6, size=1) + 
    geom_point(pch=24, size=2) + 
    scale_color_brewer("Year", palette="Set1", type="qual") + 
    theme_minimal() + ggtitle(chem)
  ggsave(filename=paste0("./regressions/regression",chem,".png"), width=8, height=6)
})

# Regressions for 2019 before transfer line replacement
meta.mix %>% filter(Year=="2019 ", RunDateP < ymd("2019-07-08"), 
                    !DateBatch %in% c("2019-06-21 Indole tests1","2019-06-19 ")) %>% 
  pivot_longer(`alpha-Pinene`:Hexenol, names_to="standard", values_to="Area") %>% filter(Area>0) %>% 
  write_csv("./metadata/regressions_2019_prexfer.csv")
  
meta.mix.early2019.long <- read_csv("./metadata/regressions_2019_prexfer_filtered.csv") %>% 
  mutate(exclude = factor(replace_na(exclude,0)), standard=str_replace_all(standard, "[- ]","_"))

meta.mix.early2019 <- meta.mix.early2019.long %>% filter(exclude=="0") %>% 
  pivot_wider(names_from="standard", values_from="Area") %>% 
  mutate(mass_ng = Conc * 3 * 1000) # (conc. ug/uL) * (3 uL) * (1000 ng / ug)

slopes.early2019 <- tibble(standard =meta.mix.early2019.long %>% 
  filter(exclude=="0") %>% pull(standard) %>% unique %>% set_names) %>% 
  mutate(model = map(standard, ~ lm(formula(paste(.x, "~ mass_ng + 0")), data=meta.mix.early2019)),
         area_per_ng = map_dbl(model, coef) %>% round(),
         r2 = map_dbl(model, ~ summary(.x)$r.squared) %>% round(3))
slopes.early2019 %>% select(-model) %>% 
  write_csv("./metadata/regressions_2019_prexfer_filtered_slopes.csv")
  
ggplot(meta.mix.early2019.long, aes(x=Conc, color=DateBatch, y=Area, shape=factor(exclude))) + 
  facet_wrap(vars(standard), scales="free_y")+
  scale_x_continuous("Concentration (mg/mL)") + scale_y_continuous("Area") + 
  scale_shape_manual("Exclude", values=c(19,1)) +
  geom_point(size=2) + geom_path(aes(group=paste(DateBatch, Replicate)))+
  #geom_abline(data = slopes.early2019, aes(slope=slope, intercept=1e-6)) +
  theme_minimal() 

# Regressions from files picked by DRC, although she added manual integrations:
# 2018 post leak fix, pre column change, and 2019 post transfer line and fix
meta.mix.DC.long <- bind_rows(
  meta.mix %>% filter(Year=="2018 ", RunDate %in% c("2018-07-15","2018-08-06","2018-08-07") &
                        !FileName %in% c("5STD3uL0.1_07152018.qgd",
                                         "STD3uL0.0001_08062018.qgd",
                                         "STD3uL0.0001try2_08062018.qgd",
                                         "STD3uL0.01rep2_08062018.qgd")), 
  meta.mix %>% filter(Year=="2019 ", 
                      RunDate %in% c("2019-07-24","2019-07-25"), 
                      Type %in% c("Standards","PineneOctanone"))) %>% 
  pivot_longer(`alpha-Pinene`:Hexenol, names_to="standard", values_to="Area") %>% 
  left_join(read_csv("./metadata/standards corrected - DC for export.csv") %>% 
               pivot_longer(`alpha-Pinene`:Hexenol, names_to="standard", values_to="Area_DC")) %>% 
  arrange(RunDateP, Conc) %>% filter(!standard %in% c("alpha-Farnesene", "Benzaldehyde")) %>% 
  mutate(diff_DC = Area_DC - Area) %>% 
  write_csv("./metadata/regressions_20182019_DC.csv")
meta.mix.DC <- meta.mix.DC.long %>% pivot_wider(names_from="standard", values_from="Area_DC") %>% 
  mutate(mass_ng = Conc * 3 * 1000)

library(RColorBrewer)
ggplot(meta.mix.DC.long %>% filter(Area_DC>0) %>% mutate(diffs=ifelse(diff_DC==0,"same","different")), 
       aes(x=Conc, color=DateBatch, y=Area_DC, shape=diffs)) + 
  facet_wrap(vars(standard))+
  scale_x_log10("Concentration (mg/mL)") + scale_y_log10("Area") + 
  scale_shape_manual("Match", values=c(1,19,2)) + 
  scale_color_manual(values=brewer.pal(6, "Dark2")[c(2,3,4,1,5,6)])+
  geom_point(data=meta.mix.DC.long %>% filter(Area>0), aes(y=Area, shape="skip"), size=2) + 
  geom_point(size=2) + geom_path(aes(group=paste(DateBatch, Replicate)))+
  geom_segment(aes(xend=Conc, yend=Area), size=1, color="black") + theme_minimal()
library(plotly)
ggplotly()

slopes.DC.18 <- tibble(standard =meta.mix.DC.long %>% pull(standard) %>% unique) %>% 
  mutate(model = map(standard, ~ lm(formula(paste0("\`",.x,"\`", "~ mass_ng + 0")), data=meta.mix.DC %>% filter(Year=="2018 "))),
         area_per_ng = map_dbl(model, coef) %>% round(),
         r2 = map_dbl(model, ~ summary(.x)$r.squared) %>% round(3))
slopes.DC.18 %>% select(-model) %>% 
  write_csv("./metadata/regressions_2018_DC_slopes.csv")

slopes.DC.19 <- tibble(standard =meta.mix.DC.long %>% pull(standard) %>% unique) %>% 
  mutate(model = map(standard, ~ lm(formula(paste0("\`",.x,"\`", "~ mass_ng + 0")), data=meta.mix.DC %>% filter(Year=="2019 "))),
         area_per_ng = map_dbl(model, coef) %>% round(),
         r2 = map_dbl(model, ~ summary(.x)$r.squared) %>% round(3))
slopes.DC.19 %>% select(-model) %>% 
  write_csv("./metadata/regressions_2019_DC_slopes.csv")

# Regressions for 2021 before filament replacement
meta.mix.2021.long <- meta.mix %>% filter(Year=="2021 ", Type =="Standards", !FileName %in% c("STD3uL0.0001rep2_20210628_3.qgd")) %>% 
  rename(alphaPinene = `alpha-Pinene`, Methylsalicylate = `Methyl salicylate`,alphaFarnesene=`alpha-Farnesene`) %>% 
  pivot_longer(alphaPinene:Hexenol, names_to="standard", values_to="Area") %>% # %>% filter(Area>0) 
  filter(!(Batch=="2021fresh" & standard == "Methylsalicylate"))

meta.mix.2021 <- meta.mix.2021.long %>% pivot_wider(names_from="standard", values_from="Area") %>% 
  mutate(mass_ng = Conc * 3 * 1000)

slopes.2021 <- tibble(standard =meta.mix.2021.long %>% pull(standard) %>% unique %>% set_names) %>% 
  mutate(model = map(standard, ~ lm(formula(paste(.x, "~ mass_ng + 0")), data=meta.mix.2021)),
         area_per_ng = map_dbl(model, coef) %>% round(),
         ng_per_area = (1/map_dbl(model, coef)) %>% round(10),
         r2 = map_dbl(model, ~ summary(.x)$r.squared) %>% round(3)) %>% 
  filter(!standard %in% c("Benzaldehyde","Farnesol"))
slopes.2021 %>% select(-model) %>% 
  write_csv("./metadata/regressions_2021_filtered_slopes.csv")

# Regressions for 2021 after filament replacement
meta.mix.2021fil2.long <- meta.mix %>% filter(Year=="2021  fil2", Type =="Standards") %>% 
  rename(alphaPinene = `alpha-Pinene`, Methylsalicylate = `Methyl salicylate`,alphaFarnesene=`alpha-Farnesene`) %>% 
  pivot_longer(alphaPinene:Hexenol, names_to="standard", values_to="Area") # %>% filter(Area>0) 

meta.mix.2021fil2 <- meta.mix.2021fil2.long %>% pivot_wider(names_from="standard", values_from="Area") %>% 
  mutate(mass_ng = Conc * 3 * 1000)

slopes.2021fil2 <- tibble(standard =meta.mix.2021fil2.long %>% pull(standard) %>% unique %>% set_names) %>% 
  mutate(model = map(standard, ~ lm(formula(paste(.x, "~ mass_ng + 0")), data=meta.mix.2021fil2)),
         area_per_ng = map_dbl(model, coef) %>% round(),
         ng_per_area = (1/map_dbl(model, coef)) %>% round(10),
         r2 = map_dbl(model, ~ summary(.x)$r.squared) %>% round(3)) %>% 
  filter(!standard %in% c("Benzaldehyde","Farnesol"))
slopes.2021fil2 %>% select(-model) %>% 
  write_csv("./metadata/regressions_2021fil2_filtered_slopes.csv")

# Combine years
convert.2021 <- c("alpha-Pinene"="alphaPinene", "Methyl salicylate"= "Methylsalicylate", "alpha-Farnesene"="alphaFarnesene")
convert.all <- c("alpha_Pinene" = "alpha-Pinene", "Methyl_salicylate" = "Methyl salicylate", "Farnesol" = "alpha-Farnesene")
bind_rows(`2018`=slopes.DC.18, `2019`=slopes.DC.19, 
          `2021`=     slopes.2021     %>% mutate(standard=fct_recode(standard, !!!convert.2021)),
          `2021fil2`= slopes.2021fil2 %>% mutate(standard=fct_recode(standard, !!!convert.2021)),
          .id="year") %>% 
  mutate(standard=fct_recode(standard, !!!convert.all)) %>% 
  select(-model, -ng_per_area) %>% write_csv("./metadata/regressions_181921_filtered_slopes.csv")

# Regressions for 2023 
meta.mix.2023.long <- meta.mix %>% filter(Year=="2023 ", Type =="Standards") %>% 
  rename(alphaPinene = `alpha-Pinene`, Methylsalicylate = `Methyl salicylate`,alphaFarnesene=`alpha-Farnesene`) %>% 
  pivot_longer(alphaPinene:Hexenol, names_to="standard", values_to="Area") # %>% filter(Area>0) 

meta.mix.2023 <- meta.mix.2023.long %>% pivot_wider(names_from="standard", values_from="Area") %>% 
  mutate(mass_ng = Conc * 3 * 1000)

slopes.2023 <- tibble(standard =meta.mix.2023.long %>% pull(standard) %>% unique %>% set_names) %>% 
  mutate(model = map(standard, ~ lm(formula(paste(.x, "~ mass_ng + 0")), data=meta.mix.2023)),
         area_per_ng = map_dbl(model, coef) %>% round(),
         ng_per_area = (1/map_dbl(model, coef)) %>% round(10),
         r2 = map_dbl(model, ~ summary(.x)$r.squared) %>% round(3)) %>% 
  filter(!standard %in% c("Benzaldehyde"))
slopes.2023 %>% select(-model) %>% 
  write_csv("./metadata/regressions_2023_slopes.csv")

####Area/Height######
stds.A.H <- dcast(stds.data, Filename~Name, median, value.var="A.H") #Area/Height
lapply(mix, function(chem) {
ggplot(cbind(stds.meta, stds.A.H) %>% filter(Type %in% ifelse(chem=="Indole",  c("Standards", "Indole"), "Standards")), aes_string(x="RunDate", y=as.name(chem), color="log10(Conc)"))  + geom_point() + scale_color_viridis() + theme_dark() + theme(axis.text.x = element_text(angle=90, hjust=-1)) + geom_line(aes(group=factor(Conc))) + ggtitle(paste(chem,"Area/Height"))
  ggsave(filename=paste0("./height/height_",chem,".png"), width=8, height=6)
})

####Area######
lapply(mergedmix, function(chem) {
  ggplot(meta.mix %>% filter(Type %in% switch(chem, "Indole"=c("Indole","Standards"), "alpha-Pinene"=c("PineneOctanone", "Standards"), "Standards")), aes_string(x="RunDate", y=as.name(chem), color="log10(Conc)"))  + geom_point() + scale_color_viridis() + theme_dark() + theme(axis.text.x = element_text(angle=90, hjust=-1)) + scale_y_log10() + geom_line(aes(group=factor(Conc))) + ggtitle(paste(chem,"Area"))
  ggsave(filename=paste0("./area/area_",chem,".png"), width=8, height=6)
})


######Alkanes#####
alk.data <- read_csv("metadata/alkanes.csv") %>% left_join(stds.data) %>% 
  filter(!grepl("Adams", Filename),grepl("Alk", Filename, ignore.case=T), n < 20) %>% 
  select(-Synonyms) %>% write_csv("metadata/alkanes_raw.csv")

alk.data <- read_csv("metadata/alkanes_fixed.csv")
alk.data.median <- alk.data %>% drop_na(n) %>% group_by(n) %>% summarize(Ret.Time=median(Ret.Time)) 
alk.data.median %>% mutate(Ret.Time=round(Ret.Time,2)) %>% write_csv("metadata/alkanes_ladder.csv")
get_kovats <- with(alk.data.median, approxfun(x = Ret.Time, y = n * 100)) #Linear interpolation
save(alk.data.median, get_kovats, file="metadata/get_kovats.rda")
ggplot(alk.data, aes(x=n, y=Ret.Time, color=Filename))  + scale_x_continuous(breaks=1:26) + scale_y_continuous(breaks=1:20) +
  geom_point() + geom_line(size=0.3) +   geom_line(data=alk.data.median, color="black")

alk.data %>% left_join(alk.data.median %>% rename(Ret.Time.Median = Ret.Time)) %>% mutate(offset=Ret.Time-Ret.Time.Median) %>%
                         ggplot(aes(x=n, y=offset, color=Filename)) + geom_line() + ylim(c(-0.25,0.25)) + geom_point()

ggplot(alk.data, aes(x=Ret.Time, y=Filename))  + scale_x_continuous(breaks=1:26) + 
 geom_rug(sides="b") + geom_text(aes(label=n, color=factor(n)), fontface=2) +theme_minimal()


alk.data.all <- stds.data %>% filter(!grepl("Adams", Filename),grepl("Alk", Filename, ignore.case=T), 
                                     Ret.Time < 20, Ret.Time > 2.5, Height > 800000)
ggplot(alk.data.all, aes(x=Ret.Time, y=Height, color=Filename, 
                         shape=ifelse(Name %in% read_csv("metadata/alkanes.csv")$Name, "Alkane", 
                                      ifelse(Name=="", "No ID","Not alkane")))) +
  geom_point() + geom_line() + geom_rug(sides="b") + geom_text(aes(label=Name)) + 
  scale_x_continuous(breaks=1:26) +  scale_shape_manual("Alkane", values=c(19,4,1)) + scale_y_sqrt()

