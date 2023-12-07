options(warn=-1)
suppressPackageStartupMessages(library('readr'))
suppressPackageStartupMessages(library('dplyr'))

# Capture command line arguments
args = commandArgs(trailingOnly = TRUE)

# Assuming BASE is the first argument
BASE = args[1]

# Construct var1 and var2 using BASE
rda_path = file.path(BASE, "data", "raw", "tftargets", "tftargets.rda")
tmp_csv_path = file.path(BASE, "data", "tmp", "tftargets", "tftargets.csv")

# Check if the file already exists
if (!file.exists(rda_path)) {
    # File doesn't exist, proceed with download
    download.file('https://github.com/slowkow/tftargets/raw/master/data/tftargets.rda', rda_path)
    message("tftargets downloaded successfully.")
} else {
    # File already exists, skip download
    message("tftargets already exists. Skipping download.")
}

load(rda_path)

#  ENCODE, ITFP, Marbach2016, Neph2012, TRED, TRRUST
all_tf_df = data.frame(
  Source = c(),
  TF = c(),
  Gene = c()
)

tf_dataset_single = c(
  "ENCODE", # entrez
  "ITFP", #
  "Marbach2016",
  "TRED", # entrez
  "TRRUST")
# tf_dataset_double = c("Neph2012")

for (dataset in tf_dataset_single){
  temp_dataset = get(dataset)
  for (i in 1:length(temp_dataset)){
    # row_count = row_count + 1
    TF_dataset = dataset
    TF_name = names(temp_dataset)[i]
    TF_gene = paste0(temp_dataset[i][[1]],collapse = ",")
    temp_df = data.frame(
      Source = TF_dataset,
      TF = TF_name,
      Gene = TF_gene
    )
    all_tf_df = rbind(all_tf_df,temp_df)
    #all_tf_df[nrow(all_tf_df) + 1,] = c(TF_dataset,TF_name,TF_gene)
  }
}
all_tf_df_merged = all_tf_df %>% 
  group_by(TF, Gene) %>% 
  summarize(
    Source = paste(Source, collapse = ", ")
  )

write_csv(all_tf_df_merged, tmp_csv_path, col_names=TRUE)
