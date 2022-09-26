list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)

# list all paths to metrics_summary.csv
# it looks into all the folders in the current working directory to find metrics_summary.csv files
metrics_files <- list.files(pattern = "metrics_summary.csv$", recursive = TRUE)

# read the first file in the list
file1 <- fread(metrics_files[1], header = TRUE)

# concatenate the other files to the first file
for (file in metrics_files[2:length(metrics_files)]){
    other_file <- fread(file, header=TRUE)
    file1 <- rbind(file1, other_file)
}

# trim file name to be only sample name
samples <- lapply(metrics_files, function(x) strsplit(x, "/")[[1]][1])

# add sample name as first column
combined_table <- cbind(samples, file1)

# save file
fwrite(x = combined_table, file =  "cellranger_count_metrics_allsamples.tsv", sep = "\t")
