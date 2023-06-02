is_ncbi_nr_installed <- function(folder_path) {
  if (!fs::dir_exists(folder_path))
    stop("The folder '", folder_path, " does not seem to exist!", call. = FALSE)

  message("Checking if a DIAMOND2 formatted ncbi nr  database with name 'nr_diamond' is available:")
  if (!fs::file_exists(file.path(folder_path, "nr_diamond"))) {
    stop("", call. = FALSE)
  } else {
    message("[................] Check passed successfully!")
  }
  message("\n")
  message("\n")

  message("Checking if a taxonmap file with name 'prot.accession2taxid' is available:")
  if (!fs::file_exists(file.path(folder_path, "prot.accession2taxid"))) {
    stop("", call. = FALSE)
  } else {
    message("[................] Check passed successfully!")
  }
  message("\n")
  message("\n")

  message("Checking if a taxonnodes file with name 'nodes.dmp' is available:")
  if (!fs::file_exists(file.path(folder_path, "nodes.dmp"))) {
    stop("", call. = FALSE)
  } else {
    message("[................] Check passed successfully!")
  }
  message("\n")
  message("\n")


}
