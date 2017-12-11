#!/usr/bin/RScript
basedir <- "/home/htseq_counts/"
setwd(basedir)

counts_input_dir <- paste(basedir, sep="/")
pat <- ".txt"
star_alignments.all <- list.files(path = counts_input_dir,
                         pattern = pat,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

myfiles <- as.vector(star_alignments.all)
data_array <- list()

# read each file as array element of data_array and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
  infile = paste(counts_input_dir, myfiles[i], sep = "/")
  data_array[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*)_all_counts.txt", "\\1", myfiles[i])
  colnames(data_array[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
data <- data_array[[myfiles[1]]]

# Join all the other tables by ID
for (i in 2:length(myfiles)) {
  y <- data_array[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# The ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

## add total counts per sample
data <- rbind(data, tot.counts=colSums(data))

####################################
# take summary rows to a new table
# ( not starting with ENS with invert=TRUE )

# transpose table for readability
data.all.summary <- data[grep("^ENS", rownames(data), perl=TRUE, invert=TRUE), ]

# transpose table
t(data.all.summary)

# write summary to file
write.csv(data.all.summary, file = "htseq_counts_all-summary.csv")

####################################
# take all data rows to a new table

data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]

# write data to file
write.csv(data.all, file = "htseq_counts_all.csv")

# cleanup intermediate objects
rm(y, z, i, data_array)

