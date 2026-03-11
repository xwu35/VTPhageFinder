#!/usr/bin/env Rscript

# this script is used to combine the assembly statistics from quast

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))

option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="path to quast results directory", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output combined table [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE);
opt = parse_args(opt_parser);

# check input
if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("path to quast results directory must be supplied (input)", call.=FALSE)
}

# list files
quast_files <- list.files(path = opt$dir, 
                          pattern = "transposed_report.tsv", 
                          recursive =  TRUE, 
                          full.names = TRUE)

# read in table if it is not empty
df_lists <- lapply(quast_files, function(x){
  if (file.info(x)$size > 0) {
    df <- fread(x, header = TRUE)
    return(df)
  } else {
    message(paste("Skippig empty file:", x))
    return(NULL)
  }
})

# combine all tables
results <- do.call(rbind, df_lists)

# save results
write.table(results, opt$out,
            sep="\t", row.names = F, col.names=T, quote = F)