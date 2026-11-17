

# load packages
packages <- c("doParallel", "dplyr", "purrr",
              "data.table", "taxonomizr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages],
                   repos = "http://cran.us.r-project.org")
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Set up SQL database for taxonomizer function:

prepareDatabase('accessionTaxa.sql', getAccessions = FALSE)

registerDoParallel(cores=24)

# reads output of daa-meganizer (.csv) as data.table
# params: fp - path to file
# assigns taxonomy to taxids
# returns: a data.table with counts, taxids, and taxonomy
taxa_to_dt <- function(fp){
  if (file.info(fp)$size == 0){
    counts_dt <- data.table(V1 = 1, V2 = 0)
  } else {
  counts_dt <- fread(fp) 
  }
  setnames(counts_dt, old = c("V1", "V2"),
           new = c("taxid", basename(fp)))
  print("Assigning Taxonomy.")
  print(basename(fp))
  mat_taxaID <- getTaxonomy(counts_dt$taxid, 'accessionTaxa.sql')  # matrix
  taxaID <- as.data.table(mat_taxaID)
  taxaID[, taxid := counts_dt$taxid]
  setkey(taxaID, taxid)
  setkey(counts_dt, taxid)
  returns = list("taxonomy" = taxaID, 
                 "otus" = counts_dt)
  return(returns) # returns a data.table
}

# where to get files from:
#file_paths <- unlist(snakemake@input)
#out_path <- snakemake@output[[1]]
file_paths <- list.files(path = "~/viridi_workflow/results/tax_files", 
                         full.names = TRUE, pattern = ".+_me00001_sup1_tp01")
out_path <- "~/viridi_workflow/results/ft/me00001_sup1_tp01.csv"

my_otu_table <- data.table(taxid=numeric(),
                           key='taxid')
my_tax_table <- data.table(taxid=numeric(),
                           superkingdom=character(),
                           phylum=character(),
                           class=character(),
                           order=character(),
                           family=character(),
                           genus=character(),
                           species=character(),
                           key='taxid')



results <- foreach (i = seq_along(file_paths)) %dopar% {
  taxa_to_dt(file_paths[i])
} # output is a list of counts tables

otu_dt_list <- results %>% map("otus") # make list of all otu dts
keyed_dts <- lapply(otu_dt_list, setkey, 'taxid')
my_otu_table <- Reduce(function(...) merge(..., all=TRUE), keyed_dts)

tax_dt_list <- results %>% map("taxonomy") #list of all tax dts
my_tax_table <- rbindlist(tax_dt_list)
my_tax_table <- unique(my_tax_table)

setkey(my_tax_table, taxid)
otu_w_tax <- my_tax_table[my_otu_table]

fwrite(otu_w_tax, out_path)
