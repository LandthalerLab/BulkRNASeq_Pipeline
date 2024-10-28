# library
library(dplyr)

# set arguements
args = commandArgs(trailingOnly=TRUE)

# list of files
list_of_table <- strsplit(args[1], split = " ")[[1]]

# merge tables
alldata <- NULL
for(i in 1:length(list_of_table)){
  print(list_of_table[i])
  cur_data <- read.table(list_of_table[i], header = T, sep = "\t")
  if(i == 1){
    alldata <- cur_data
  } else {
    cur_data <- cur_data[,!colnames(cur_data) %in% "Total.Reads"]
    alldata <- alldata %>% full_join(cur_data, by = c("Sample" = "Sample"))
  }
}

# write table
alldata <- alldata[,!apply(alldata,2,function(x) all(is.na(x)))]
write.table(alldata, args[2], quote = FALSE, row.names = FALSE, sep = "\t")
