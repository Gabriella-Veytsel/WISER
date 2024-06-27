library(tidyverse)
library(ggtree)
library(ggplot2)
library(readxl)
library(phylotools)

source("C:/Users/Gabriella Veytsel/OneDrive/Documents/Clean/functions.R")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")

#########
#Clinical
#########

#UGA Clinical: 557 samples
uga_clinical_samples <- read_excel("D:/wastewater/Clinical_sample_metadata.xlsx") #includes metadata without associated samples
uga_clinical_samples <- uga_clinical_samples %>%
  mutate(seq.name = paste(Sample, lineage, date, `S#`, sep = "|"))

uga_clinical_samples_fasta <- read.fasta("D:/wastewater/analysis2/clinical/all_clincal_seqs_N50_aligned_formatted.fasta")
uga_clinical_samples_join <- uga_clinical_samples_fasta %>% inner_join(uga_clinical_samples)
uga_clinical_samples_join <- uga_clinical_samples_join %>%
  mutate(zip = as.double(zip))

uga_clinical_samples_join_who <- who_name(uga_clinical_samples_join)
uga_clinical_samples_join_who$WHO_name <- ifelse(is.na(uga_clinical_samples_join_who$WHO_name), "Other", uga_clinical_samples_join_who$WHO_name) #No WHO name
uga_clinical_samples_join_who$WHO_name <- as.factor(uga_clinical_samples_join_who$WHO_name)

range(uga_clinical_samples_join$date) #dates range from 2021-03-10 to 2022-02-25
ggplot(uga_clinical_samples_join, aes(x=date)) + geom_histogram()
ggsave("D:/wastewater/figures/UGA Sequences over Time.pdf")

ggplot(uga_clinical_samples_join_who, aes(x=WHO_name)) + geom_histogram(stat="count")
ggsave("D:/wastewater/figures/UGA Variants.pdf", width=10, height=7)

#GISAID
metadata_district10 <- read_tsv("D:/georgia/metadata_district10.tsv")
range(metadata_district10$date) #dates range from 2020-12-26 to 2022-10-19

metadata_district10 %>% group_by(county) %>% summarise(n=n()) #680 samples in Clarke County
metadata_clarke <- metadata_district10 %>% filter(county == "Clarke County")

#If I match the dates of UGA clinical samples
metadata_clarke <- metadata_clarke %>%
  filter(date <= as.Date("2022-02-25")) %>%
  filter(date >= as.Date("2021-03-10")) #680 samples in Clarke County -> 530 samples in Clarke County during same time period
metadata_clarke$WHO_name <- as.factor(metadata_clarke$WHO_name)
#metadata_clarke %>% group_by(major_city) %>% summarise(n=n()) #490 samples in Athens
metadata_clarke %>% group_by(coverage) %>% summarise(n=n()) 

ggplot(metadata_clarke, aes(x=date)) + geom_histogram()
ggsave("D:/wastewater/figures/GISAID Sequences over Time.pdf")

ggplot(metadata_clarke, aes(x=WHO_name)) + geom_histogram(stat="count")
ggsave("D:/wastewater/figures/GISAID Variants.pdf", width = 10, height=7)

#Compare
metadata_clarke$dataset = "GISAID"
uga_clinical_samples_join_who$dataset = "UGA"
combined_data <- bind_rows(metadata_clarke, uga_clinical_samples_join_who)

ggplot(combined_data, aes(x = date, fill = dataset)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Overlayed Histograms", x = "date", y = "Frequency") +
  theme_minimal()
ggsave("D:/wastewater/figures/UGA vs. GISAID Sequences over Time.pdf")

ggplot(combined_data, aes(x = WHO_name, fill = dataset)) + 
  geom_histogram(stat="count", bins = 30) +
  labs(title = "Stacked Bar Chart", x = "Variant", y = "Frequency") +
  theme_minimal()
ggsave("D:/wastewater/figures/UGA vs. GISAID Variants.pdf", height=7, width=12)

#Get the FASTA for GISAID Clarke County same time period
write_tsv(metadata_clarke, "D:/wastewater/metadata_clarke_sametimeperiod.tsv")

complete <- read.fasta("D:/georgia old/gisaid/complete delta genomes/gisaid_sequences.combined.fasta")
clarke <- metadata_clarke %>% left_join(complete, by = c("strain" = "seq.name"))
clarke$seq.name <- paste(clarke$strain, clarke$date, clarke$pangolin_lineage, clarke$WHO_name, sep = "|")
dat2fasta(clarke, "D:/wastewater/analysis2/clinical gisaid full length/gisaid_clarke.fasta")

###########
#Analysis 3
###########
ww <- read.fasta("D:/wastewater/ww full length/ww_full_removeNs.fasta")
ww_delta <- ww %>%
  filter(grepl("AY|B.1.617.2", seq.name)) #15
dat2fasta(ww_delta, "D:/wastewater/analysis3/ww_full_removeNs_delta.fasta")

gisaid <- read.fasta("D:/wastewater/analysis2/clinical gisaid full length/gisaid_clarke.fasta")
gisaid_delta <- gisaid %>%
  filter(grepl("AY|B.1.617.2", seq.name)) #263
dat2fasta(gisaid_delta, "D:/wastewater/analysis3/gisaid_delta.fasta")

#I aligned these in Nextalign, so it'll be different, just pull the IDs and align in MAFFT
cluster_IDs <- read.fasta("D:/georgia/clusters/clusters_fasta.fasta") %>% pull(seq.name)
contextual <- read.fasta("D:/georgia/nextstrain/contextual_sequences.fasta")
focal_weighted <- read.fasta("D:/georgia/weighted_sub.fasta")
gisaid_input <- bind_rows(contextual, focal_weighted) 

#Trying to match
gisaid_input$seq.name <- gsub("hCoV-19/", "", gisaid_input$seq.name)
cluster_IDs <- gsub("USA_", "USA/", cluster_IDs)
cluster_IDs <- gsub("_2021/", "/2021/", cluster_IDs)
cluster_IDs <- gsub("_2022/", "/2022/", cluster_IDs)
matching <- gisaid_input %>%
  filter(seq.name %in% cluster_IDs)

cluster_IDs <- cluster_IDs %>% as.data.frame()
notmatching <- cluster_IDs %>% filter(grepl("OOS", cluster_IDs$.))
notmatching_USA <- notmatching %>% filter(grepl("USA/", notmatching$.))
notmatching_Other <- notmatching %>% filter(!grepl("USA/", notmatching$.))
notmatching_Other$. <- sub("_", "/", notmatching_Other$.) #Sub() matches only first occurence of a pattern
#After troubleshooting, found that Czech_Republic is causing trouble
notmatching_Other$. <- ifelse(notmatching_Other$. == "Czech/Republic_25_9_2021_89/2021/2021-09-20/38/OOS",
                                                      "Czech_Republic/25_9_2021_89/2021/2021-09-20/38/OOS",
                                                      notmatching_Other$.)
notmatching_Other$. <- ifelse(notmatching_Other$. == "Czech/Republic_BC_21102021_11/2021/2021-10-07/40/OOS",
                                                      "Czech_Republic/BC_21102021_11/2021/2021-10-07/40/OOS",
                                                       notmatching_Other$.)

notmatching_USA <- notmatching_USA %>% 
  separate_wider_delim(".", delim = "/", names = c("country", "ID", "year", "date", "week", "loc"))
notmatching_Other <- notmatching_Other %>% 
  separate_wider_delim(".", delim = "/", names = c("country", "ID", "year", "date", "week", "loc"))

notmatching_USA <- notmatching_USA %>%
  mutate(seq.name = paste(country, ID, year, sep = "/")) 
notmatching_Other <- notmatching_Other %>%
  mutate(seq.name = paste(country, ID, year, sep = "/")) 

notmatching <- bind_rows(notmatching_USA, notmatching_Other) 
notmatching_name <- notmatching %>% pull(seq.name)
matching2 <- gisaid_input %>%
  filter(seq.name %in% notmatching$seq.name)
gisaid_match <- rbind(matching, matching2)
gisaid_match <- gisaid_match %>% left_join(notmatching)
gisaid_match <- gisaid_match %>%
  mutate(seq.name.new = paste(seq.name, date, week, loc, sep = "/")) %>%
  mutate(seq.name.new = gsub("/NA/NA/NA", "", seq.name.new)) %>%
  select(-seq.name) %>% rename(seq.name = seq.name.new)

gisaid_clarke_delta <- read.fasta("D:/wastewater/analysis3/gisaid_delta.fasta")
gisaid_clarke_delta_ID <- gisaid_clarke_delta %>% 
  separate_wider_delim("seq.name", delim = "/", names = c("hcov", "country", "ID", "rest")) %>% pull("ID")

gisaid_match <- gisaid_match %>%
  mutate(ID = str_extract(seq.name, "(?<=/)[^/]+")) #extract everything between first and second slash

gisaid_match_clarke <- gisaid_match %>%
  filter(gisaid_match$ID %in% gisaid_clarke_delta_ID) #There's only 40 sequences from the Northeast district > 16 from Clarke County

gisaid_match <- gisaid_match %>%
  filter(gisaid_match$ID %!in% gisaid_clarke_delta_ID) #There's only 40 sequences from the Northeast district > 16 from Clarke County
dat2fasta(gisaid_match, "D:/wastewater/analysis3/clusters.fasta")

##
library(treeio)
tree2 <- read.tree("D:/wastewater/analysis2/ww full length/iqtree/ww_full_removeNs_aligned_matched4_rmgaps_trim_date.fasta.treefile")
tree2 <- root(tree2, "MN908947.3")
write.tree(tree2, "D:/wastewater/analysis2/ww full length/iqtree/ww_full_removeNs_aligned_matched4_rmgaps_trim_date_root.fasta.treefile")

tree1 <- read.tree("D:/wastewater/analysis1/all full length/a1_output_mis4_trim.fasta.timetree.nwk")
ggtree(tree1) + geom_tiplab(aes(color=grep("Alpha", label)))
tree1$tip.label
