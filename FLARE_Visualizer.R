
#### ==== CHECK and LOAD LIBRARIES ==== ####

libs <- c("data.table","tidyverse","optparse","randomcoloR","foreach","doParallel")

installed_packages <- installed.packages()[,1]

libs_installed <- libs[(libs %in% installed_packages)]
libs_not_installed <- libs[!(libs %in% installed_packages)]

#install and load packages not installed
for (lib in libs_not_installed) {
  install.packages(lib)
  if (!require(lib,character.only = TRUE)) {
    stop("Error: could not install new R package:", lib)
  }
}
#load installed packages
for (lib in libs_installed) {
  if (!require(lib,character.only = TRUE)) {
    stop("Error: could not load existing R package:", lib)
  }
}
#### ==== PARSE OPTIONS ==== ####

option_list = list(
  make_option(c("-i", "--input-file"), type="character", default=NULL, 
              help="path of the ancestry file produced by FLARE", metavar="character"),
  
  make_option(c("-d", "--do-parallel"), type="character", default=FALSE,
              action="store_true",
              help="", metavar="character"),
  
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="name of the output file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$`input-file`)){
  print_help(opt_parser)
  stop("Error: --input-file argument should be given", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("Error: --output argument should be given", call.=FALSE)
}