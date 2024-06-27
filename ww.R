library(phylotools)
library(tidyverse)
library(stringdist)
library(treeio)

####################################################################################
ww <- read.fasta("D:/wastewater/ww full length/ww_full.fasta")

#Remove Ns and gaps
ww <- ww %>% mutate(seq.text = gsub("N", "", seq.text)) 
str_detect(ww$seq.text, "-") #Looks like gaps come after alignment

dat2fasta(ww, "D:/wastewater/ww full length/ww_full_removeNs.fasta")

####################################################################################

#Analysis 1
#Read alignment
output_aligned1 <- read.fasta("D:/wastewater/analysis1/all full length/mafft/output.fasta")
ref_seq1 <- output_aligned1[1, ]
ww1 <- output_aligned1 %>% filter(grepl("^W", output_aligned1$seq.name)) #ww only
ww1 <- rbind(ref_seq1, ww1) #ww + ref
output_aligned1 <- output_aligned1 %>% filter(!grepl("^W", output_aligned1$seq.name)) #clinical only

#Analysis 2
ww2 <- read.fasta("D:/wastewater/analysis2/ww full length/output.fasta")
ref_seq2 <- ww2[1, ]

#Analysis 3
#Read alignment
output_aligned3 <- read.fasta("D:/wastewater/analysis3/output.fasta")
ref_seq3 <- output_aligned3[1, ]
ww3 <- output_aligned3 %>% filter(grepl("^W", output_aligned3$seq.name)) #ww only
ww3 <- rbind(ref_seq3, ww3) #ww + ref
output_aligned3 <- output_aligned3 %>% filter(!grepl("^W", output_aligned3$seq.name)) #clinical only

#Read alignment - subset
output_aligned3 <- read.fasta("D:/wastewater/analysis3/output_subset.fasta")
ref_seq3 <- output_aligned3[1, ]
ww3 <- output_aligned3 %>% filter(grepl("^W", output_aligned3$seq.name)) #ww only
ww3 <- rbind(ref_seq3, ww3) #ww + ref
output_aligned3 <- output_aligned3 %>% filter(!grepl("^W", output_aligned3$seq.name)) #clinical only

#Frankenstein detection
frankenseq <- function(data, window_size) {
  # Extract the reference sequence from the first row of the dataset
  refseq <- data[1, 2]
  # Determine the number of rows in the dataset
  n <- nrow(data)
  # Create an empty list to store the compared sequences
  compared_seqs <- vector("list", n - 1)
  
  # Iterate through each row of the dataset (excluding the first row, which contains the reference sequence)
  for (i in 2:n) {
    # Extract the test sequence from the current row
    test.seq <- data[i, 2]
    # Calculate the total length of the reference sequence
    total_length <- nchar(refseq)
    # Calculate the number of iterations needed for the sliding window approach
    stop_length <- total_length - window_size + 1
    # Initialize a new sequence variable with the test sequence
    new.seq <- test.seq
    
    # Iterate through each position in the test sequence using a sliding window approach
    for (j in 1:stop_length) {
      # Extract the reference and test windows for comparison
      ref.window <- substr(refseq, j, j + window_size - 1)
      test.window <- substr(test.seq, j, j + window_size - 1)
      # Calculate the Levenshtein distance between the reference and test windows
      distance <- stringdist::stringdist(ref.window, test.window, method = "lv")
      
      # If the distance is greater than or equal to 4, replace the corresponding portion of the test sequence with gaps
      if (distance >= 4) {
        new.seq <- paste0(substr(new.seq, 1, j - 1), strrep("-", window_size), substr(new.seq, j + window_size, nchar(new.seq)))
      }
    }
    
    # Store the modified test sequence in the list of compared sequences
    compared_seqs[[i - 1]] <- new.seq
  }
  
  # Remove the first row (reference sequence) from the original dataset
  new_data <- data[-1, ] 
  # Add a new column to the dataset containing the compared sequences
  new_data$compared_sequences <- compared_seqs
  # Return the modified dataset
  return(new_data)
}

#Analysis 1
data_with_compared_seqs1 <- frankenseq(ww1, 10)
data1 <- data_with_compared_seqs %>% select(seq.name, compared_sequences) %>% rename("seq.text" = "compared_sequences")
data1 <- rbind(data1, output_aligned1)
dat2fasta(data1, "D:/wastewater/analysis1/all full length/mafft/output_mismatch4.fasta")

#Analysis 2
data_with_compared_seqs2 <- frankenseq(ww2, 10)
data2 <- data_with_compared_seqs2 %>% select(seq.name, compared_sequences) %>% rename("seq.text" = "compared_sequences")
data2 <- rbind(data2, ref_seq2)
dat2fasta(data2, "D:/wastewater/analysis2/ww full length/output_mismatch4.fasta")

#Analysis 3
data_with_compared_seqs3 <- frankenseq(ww3, 10)
data3 <- data_with_compared_seqs3 %>% select(seq.name, compared_sequences) %>% rename("seq.text" = "compared_sequences")
data3 <- rbind(data3, output_aligned3)
#dat2fasta(data3, "D:/wastewater/analysis3/output_mismatch4.fasta")
dat2fasta(data3, "D:/wastewater/analysis3/output_mismatch4_subset.fasta")

#Analysis 1
tree1 <- read.tree("D:/wastewater/analysis1/all full length/iqtree/output_mismatch4_dates.fasta.treefile")
tree1 <- root(tree1, "MN908947.3")
write.tree(tree1, "D:/wastewater/analysis1/all full length/iqtree/output_mismatch4_dates_rooted.fasta.treefile")

#Analysis 3
tree3 <- read.tree("D:/wastewater/analysis3/iqtree/output_mismatch4_dates.fasta.treefile")
tree3 <- root(tree3, "MN908947.3")
write.tree(tree3, "D:/wastewater/analysis3/iqtree/output_mismatch4_dates_rooted.fasta.treefile")

#TMRCA
time_length(ymd("2022-09-05")-ymd("2020-08-28"), "years")
time_length(ymd("2022-09-05")-ymd("2020-10-19"), "years")
time_length(ymd("2022-09-05")-ymd("2021-10-09"), "years")
time_length(ymd("2022-09-05")-ymd("2019-12-26"), "years")
time_length(ymd("2022-09-05")-ymd("2019-12-26"), "years")

time_length(ymd("2022-02-24")-ymd("2020-08-28"), "years")
time_length(ymd("2022-02-24")-ymd("2020-10-19"), "years")
time_length(ymd("2022-02-24")-ymd("2021-10-09"), "years")
time_length(ymd("2022-02-24")-ymd("2019-12-26"), "years")
time_length(ymd("2022-02-24")-ymd("2019-12-26"), "years")
