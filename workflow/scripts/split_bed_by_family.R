#
# 1_split_bed.r
#
# This script reads a master BED file, filters it, and splits it
# into multiple individual BED files based on the repeat 'family'.
# It is designed to be called from Snakemake and accepts command-line arguments.
#
# Example:
# Rscript 1_split_bed.r --input_bed /app/Data/input.bed --output_dir /app/Data/bed_type
#

library(argparse)
library(data.table)
library(tidyverse)
library(stringr) # For str_split

# --- 1. Define Command-Line Arguments ---
parser <- ArgumentParser(description = "Split a master BED file into per-family BED files.")
parser$add_argument("--input_bed", type = "character", required = TRUE,
                    help = "Path to the input master BED file (e.g., 250512_CHM13_RPMSK.bed)")
parser$add_argument("--output_dir", type = "character", required = TRUE,
                    help = "Path to the directory to save the split BED files")
parser$add_argument("--intermediate", type = "character", required = TRUE,
                    help = "Path to the directory to save the intermediate summery files")

args <- parser$parse_args()

# --- 2. Create Output Directory ---
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}
if (!dir.exists(args$intermediate)) {
  dir.create(args$intermediate, recursive = TRUE)
}

# --- 3. Load and Process Data ---
print(paste("Loading BED file:", args$input_bed))
dat <- fread(args$input_bed)

# Define column names for the input BED
bed_col_names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "ConsensusStart", "ConsensusEnd", "filler")
colnames(dat) <- bed_col_names[1:ncol(dat)]

# Select only the columns needed for the BED split (for bedtools getfasta)
data <- dat %>% select(chrom, chromStart, chromEnd, name)

# Extract family and subfamily info from the 'name' column
# MER4A1#LTR/ERV1
data$subfam<-sapply(str_split(data$name,"#"),"[",2)#as [2]
data$fam<-sapply(str_split(data$subfam,"/"),"[",1)

# --- 4. Filter out unwanted families ---

print("Filtering out 'Unknown', 'Simple_repeat', 'Low_complexity'...")
initial_rows <- nrow(data)
print(paste("Curret: ", initial_rows, "rows."))
data_filtered <- data %>%
  filter(!fam %in% c("Unknown", "Simple_repeat", "Low_complexity"))
print(paste("Removed", initial_rows - nrow(data_filtered), "rows."))

stats<-data
stats$width<-stats$chromEnd - stats$chromStart
stats<-stats %>% group_by(fam) %>% summarize(total=sum(width))
stats$total_mb<-stats$total/1000000
stats<-stats %>% select(fam,total_mb)
inte_path <- file.path(args$intermediate, paste0("1_Repeat_length_MBp.tsv"))
fwrite(stats, inte_path, sep = "\t", row.names = FALSE, quote = FALSE)
#This leaves 57 subfams from 13 families, comprising 1287 types
#937 of these appear on the y chromosome and 1179 of them appear on the x chromosome

d1<-data %>% group_by(name,`chrom`)%>%summarize(c=n())#count of same repeat per chr
d2<-data %>% group_by(name)%>%summarize(t=n())#count of same repeat in genome

d1<-inner_join(d1,d2,by="name")
d1$p<-d1$c/d1$t#proportion of each repeat name of chromosome to genome
inte_path <-  file.path(args$intermediate, paste0("1_Repeat_appearance.tsv"))
fwrite(d1, inte_path, sep = "\t", row.names = FALSE, quote = FALSE)
inst<-d1 %>% filter(`chrom`=="chrX" & p>.1) ### 106 appear on the x chr more than 10% of the time
inte_path <-  file.path(args$intermediate, paste0("1_chrX_enr_repeats.tsv"))
fwrite(inst, inte_path, sep = "\t", row.names = FALSE, quote = FALSE)
inte_path <-  file.path(args$intermediate, paste0("1_chrY_enr_repeats.tsv"))
inst<-d1 %>% filter(`chrom`=="chrY" & p>.05) ### 77 appear on the y chr more than 10% of the time
fwrite(inst, inte_path, sep = "\t", row.names = FALSE, quote = FALSE)

names<-unique(data_filtered$name)
print(paste("Found", length(names), "unique types to split into files."))
for (i in 1:length(names)) {
  d<-data_filtered %>% filter(name==names[i])
  nm<-gsub("/","_",names[i])
  nm <- gsub("[\\?*<>|:]", "", nm)
  output_filename <- file.path(args$output_dir, paste0(nm, ".bed"))
  fwrite(d,output_filename,quote=F,row.names = F,col.names = F,sep="\t")
}

print(paste("Successfully split BED file into", length(names), "files in", args$output_dir))
