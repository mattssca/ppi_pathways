# Function to download and extract TSV files
download_hpa_data <- function(url, description) {
  cat(paste0("Downloading ", description, "...\n"))
  temp_zip <- tempfile(fileext = ".zip")
  temp_dir <- tempdir()

  download.file(url, temp_zip, mode = "wb")
  unzip(temp_zip, exdir = temp_dir)

  # Get the extracted file name
  tsv_file <- list.files(temp_dir, pattern = "\\.tsv$", full.names = TRUE)
  tsv_file <- tsv_file[grep(basename(url), basename(tsv_file), ignore.case = TRUE)]

  if (length(tsv_file) == 0) {
    tsv_file <- list.files(temp_dir, pattern = "\\.tsv$", full.names = TRUE)[1]
  }

  cat(paste0("Reading data from: ", basename(tsv_file), "\n"))
  data <- read_tsv(tsv_file, show_col_types = FALSE)

  return(data)
}
