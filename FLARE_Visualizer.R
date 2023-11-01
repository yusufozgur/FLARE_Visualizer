#TODO: consider not including parallel and foreach packages to package load if this option FALSE
#TODO: instead of read_line arbitrary 50 lines and getting populations and
#      tested individuals from that text, get metadata and headers, which are lines
#      starting with ## and # respectively and extract populations and individuals
#      that way
#TODO: add / to opt$`output-dir` automatically if not added by user
#### ==== CHECK and LOAD LIBRARIES ==== ####

libs <- c("data.table", "tidyverse", "optparse", "randomcoloR", "foreach", "doParallel", "cli")

installed_packages <- installed.packages()[, 1]

libs_installed <- libs[(libs %in% installed_packages)]
libs_not_installed <- libs[!(libs %in% installed_packages)]

#install and load packages not installed
for (lib in libs_not_installed) {
	install.packages(lib)
	if (!require(lib, character.only = TRUE)) {
		stop("Error: could not install new R package:", lib)
	}
}
#load installed packages
for (lib in libs_installed) {
	if (!require(lib, character.only = TRUE)) {
		stop("Error: could not load existing R package:", lib)
	}
}

rm(installed_packages, libs_installed, libs_not_installed)#not libs, used in foreach
#### ==== PARSE ARGUMENTS ==== ####

option_list = list(
	make_option(
		c("-i", "--input-file"),
		type = "character",
		default = NULL,
		help = "path of the ancestry file produced by FLARE. it can either be a vcf file or vcf.gz file.",
		metavar = "PATH"
	),

	make_option(
		c("-d", "--do-parallel"),
		type = "character",
		default = FALSE,
		action = "store_true",
		help = "Perform parallel computing for faster processing. Turned off by default"
	),

	make_option(
		c("-o", "--output-dir"),
		type = "character",
		default = NULL,
		help = "path of the output directory",
		metavar = "PATH"
	)
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$`input-file`)) {
	print_help(opt_parser)
	stop("Error: --input-file argument should be given", call. = FALSE)
}
if (is.null(opt$`output-dir`)) {
	print_help(opt_parser)
	stop("Error: --output-dir argument should be given", call. = FALSE)
}
if (!dir.exists(opt$`output-dir`)) {
	cli_inform(paste("Output directory does not exist, creating it:", opt$`output-dir`))
	dir.create(opt$`output-dir`)
}
#TODO: sep should use platform appropriate path separator
admixture_proportions_table_path <- paste(opt$`output-dir`, "admixture_proportions.tsv", sep = "/")
if (file.exists(admixture_proportions_table_path)) {
	cli_inform("admixture_proportions.tsv already exists in the output dir, deleting it")
	unlink(admixture_proportions_table_path)
}
file.create(admixture_proportions_table_path)#create it now to be able to append later
#### ==== tagore constants ==== ####
tagore_feature = 3
tagore_size = 1

#### ==== Extract Ancestries and assign distinct colors ==== ####

#extract ancestries with their respective ids from vcf metadata
ancestries_and_colors <-
	read_lines(file = opt$`input-file`, n_max = 50) |>
	str_extract("##ANCESTRY.*>") |>
	na.omit() |>
	tail(n = 1) |>
	str_extract("<(.*)>") |>
	(\(x) gsub("[<>]", "", x))() |>
	strsplit(",") |>
	(\(x) x[[1]])() |>
	(\(x) separate_wider_delim(
		as.data.frame(x),
		x,
		delim = "=",
		names = c("name", "id")
	))() |>
	as.data.frame() |>
	mutate(id = as.numeric(id))

#assign unique colors to each ancestry
ancestries_and_colors$color <- ancestries_and_colors |> nrow() |> distinctColorPalette()

#### ==== Extract tested individuals ==== ####
inds <- read_lines(file = opt$`input-file`, n_max = 50) |>
	str_extract("#CHROM.*") |>
	na.omit() |>
	strsplit(split = "\t") |>
	(\(x) x[[1]][10:length(x) - 10])() |> # remove first 9 columns, which are the fixed columns
	sort()

#### ==== functions: visualization for a single person ==== ####
# loading the whole FLARE output, which can be in magnitude of gigabytes, and
# processing that can require too much ram, instead, creating a function to only
# load data for the designated individual and visualize chromosomes of that
# individual only, and running this function in a loop proved to be better for
# performance.

pop_id_to_pop_name <- function(id, ref_table) {
	lookup <- match(id, ref_table$id)
	return(ref_table$name[lookup])
}

process_base_data <- function(base_data) {
	base_data <- base_data |>
		rename(CHROM = `#CHROM`) |>
		select(-ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT)

	ind_cols <- colnames(base_data)[3:ncol(base_data)]

	base_data <- base_data |>
		pivot_longer(
			cols = all_of(ind_cols),
			names_to = "person",
			values_to = "inferred_ancestry"
		)

	colon_times <- str_count(base_data$inferred_ancestry[1], ":")

	if (colon_times == 2) {
		cat("Colons count is ", colon_times, ". probs=False option detected.\n")
		base_data$inferred_ancestry <-
			base_data$inferred_ancestry |>
			lapply(\(x) str_extract(as.character(x), "\\d+:\\d+$"))
	} else if (colon_times == 4) {
		cat("Colons count is ", colon_times, ". probs=True option detected.\n")
		base_data$inferred_ancestry <-
			base_data$inferred_ancestry |>
			lapply(\(x) str_extract(as.character(x), "\\d+\\|\\d+:\\d+:\\d+")) |>
			lapply(\(x) str_extract(as.character(x), "\\d+:\\d+$"))
	} else {
		cat("Colons count is ", colon_times, ". Data seems to be in an incorrect format.\n")
		stop()
	}

	return(as.data.frame(base_data))
}

create_admixture_proportions_table <- function(processed_data) {
	processed_data <- processed_data |>
		separate_longer_delim(cols = inferred_ancestry, delim = ":") |>
		count(person, inferred_ancestry, name = "counts") |>
		mutate(inferred_ancestry = pop_id_to_pop_name(inferred_ancestry, ancestries_and_colors))

	total_count_each_person = select(processed_data, person, counts) |>
		group_by(person) |>
		summarise(total_counts = sum(counts))

	processed_data <- processed_data |>
		mutate(proportions = counts / total_count_each_person$total_counts[match(person, total_count_each_person$person)])

	processed_data |>
		pivot_wider(
			id_cols = c(person),
			names_from = c(inferred_ancestry),
			values_from = proportions
		) |>
		as.data.frame() |>
		mutate_all(.funs = \(x) replace_na(x, 0))
}

prepare_data_for_tagore_inputfiles <- function(processed_data) {
	processed_data |>
		separate_wider_delim(
			cols = inferred_ancestry,
			delim = ":",
			names = c("h1", "h2")
		) |>
		mutate(
			CHROM = paste("chr", CHROM, sep = ""),
			POS2 = POS,
			feature = tagore_feature,
			size = tagore_size,
			h1 = as.numeric(h1),
			h2 = as.numeric(h2)
		) |>
		pivot_longer(
			cols = all_of(c("h1", "h2")),
			names_to = "chr_copy",
			values_to = "color"
		) |>
		mutate(chr_copy = str_sub(chr_copy, start = -1, end = -1)) |>
		mutate(color = ancestries_and_colors$color[match(color, ancestries_and_colors$id)]) |>
		relocate(CHROM, POS, POS2, feature, size, color, chr_copy) |>
		select(-person) #if single person, which is always true for this program
}

single_person_visualizer <- function(person_name) {
	cli_inform("Reading file")
	base_data <- fread(
		file = opt$`input-file`,
		skip = "#CHROM",
		select = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", person_name)
	)

	cli_inform("Processing file")
	base_data <- process_base_data(base_data)

	cli_inform("Creating admixture proportions table")
	admixture_proportions_table <- create_admixture_proportions_table(base_data)
	#For using in tagore
	nonzero_ancestries <- colnames(admixture_proportions_table)[-1]
	#fill in the 0 proportion ancestries as otherwise append will append wrong

	for (ancestry in ancestries_and_colors$name) {
		if (!(ancestry %in% colnames(admixture_proportions_table))) {
			admixture_proportions_table[ancestry] = 0
		}
	}
	#reorder rows so append will be right
	ordered_columns <-
		c("person", (colnames(admixture_proportions_table)[-1] |> sort()))
	admixture_proportions_table <-
		admixture_proportions_table[, ordered_columns]


	cli_inform("Writing admixture proportions table")
	fwrite(
		admixture_proportions_table,
		file = admixture_proportions_table_path,
		sep = "\t",
		append = T
	)

	cli_inform("Writing tagore file")

	tagore_inputfile <- prepare_data_for_tagore_inputfiles(base_data)

	tagore_inputfile_path <- paste(
		opt$`output-dir`,
		"/",
		person_name,
		".tagorefile",
		sep = ""
	)

	fwrite(
		tagore_inputfile,
		file = tagore_inputfile_path,
		sep = "\t",
		col.names = F,
		row.names = F
	)

	print("Executing tagore")

	tagore_legend <- ancestries_and_colors |>
		filter(name %in% nonzero_ancestries) |>
		mutate(merge = paste(name, color, sep = ":")) |>
		(\(x) x$merge)() |>
		paste(collapse = ",")

	paste(
		"tagore --input ",
		tagore_inputfile_path,
		" --prefix tagore_images/",
		person_name,
		" -vf -b hg38 ",
		"-l ",
		"'",
		tagore_legend,
		"'",
		sep = ""
	) |> print()
}

#### ==== Main LOOP ==== ####
#either run the single person visualization in parallel for a performance boost, or run it in a single process
if (opt$`do-parallel`) {
	registerDoParallel(cores = detectCores())
	foreach (name = inds) %dopar% {
		#load libs again
		for (lib in libs) {
			require(lib, character.only = T)
		}

		cat(name, " : ", match(name, inds), "/", length(inds), "\n")
		single_person_visualizer(name)
	}
} else {
	for (name in inds) {
		cat(name, " : ", match(name, inds), "/", length(inds), "\n")
		single_person_visualizer(name)
	}
}
