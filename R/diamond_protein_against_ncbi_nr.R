#' @title Perform Protein to Protein DIAMOND2 Searches (BLASTP) against the NCBI non-redundant (NR) database 
#' @description Run protein to protein DIAMOND2 of reference sequences
#' against a blast-able database or fasta file.
#' @param query path to input file in fasta format.
#' @param ncbi_nr_folder_path_path path to the directory that either stores the raw NCBI NR database
#' with assumed name \code{nr.gz} (requires \code{make_diamond_db = TRUE}) or to the already formatted \code{nr_diamond.dmnd} or \code{nr} database (assumes default: \code{make_diamond_db = FALSE}).
#' @param make_diamond_db logical specifying whether or not the NCBI NR database at \code{ncbi_nr_folder_path_path} should be formatted with \code{diamond makedb} (\code{make_diamond_db = FALSE}; default).
#' @param task protein search task option. Options are:
#' \itemize{
#' \item \code{task = "blastp"} : Standard protein-protein comparisons (default).
#' }
#' @param sensitivity_mode specify the level of alignment sensitivity. The higher the sensitivity level, the more deep homologs can be found, but at the cost of reduced computational speed.
#' \itemize{
#'   \item \code{sensitivity_mode = "faster"} : fastest alignment mode, but least sensitive (default). Designed for finding hits of >70% identity and short read alignment such as NGS reads.
#'   \item \code{sensitivity_mode = "default"} : Default mode. Designed for finding hits of >70% identity and short read alignment such as NGS reads.
#'   \item \code{sensitivity_mode = "fast"} : fast alignment mode, but least sensitive (default). Designed for finding hits of >70% identity and short read alignment such as NGS reads.
#'   \item \code{sensitivity_mode = "mid-sensitive"} : fast alignments between the \code{fast} mode and the sensitive mode in sensitivity.
#'   \item \code{sensitivity_mode = "sensitive"} : fast alignments, but full sensitivity for hits >40% identity.
#'   \item \code{sensitivity_mode = "more-sensitive"} : more sensitive than the \code{sensitive} mode.
#'   \item \code{sensitivity_mode = "very-sensitive"} : sensitive alignment mode.
#'   \item \code{sensitivity_mode = "ultra-sensitive"} : most sensitive alignment mode (sensitivity as high as BLASTP).
#' }
#' @param use_arrow_duckdb_connection shall DIAMOND2 hit output table be transformed to an in-process (big data disk-processing) arrow connection to DuckDB? This is useful when the DIAMOND2 output table to too large to fit into memory. Default is \code{use_arrow_duckdb_connection = FALSE}.
#' Please consult the Installation Vignette for details.
#' @param evalue Expectation value (E) threshold for saving hits (default: \code{evalue = 0.001}).
#' @param out_format a character string specifying the format of the file in which the DIAMOND results shall be stored.
#' Available options are:
#'  \itemize{
#'  \item \code{out_format = "pair"} : Pairwise
#'  \item \code{out_format = "xml"} : XML
#'  \item \code{out_format = "csv"} : Comma-separated file
#'  }
#' @param cores number of cores for parallel DIAMOND searches.
#' @param max_target_seqs maximum number of aligned sequences that shall be retained. Please be aware that \code{max_target_seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param hard_mask shall low complexity regions be hard masked with TANTAN? Default is \code{db_hard_mask = TRUE}.
#' @param diamond_exec_path a path to the DIAMOND executable or \code{conda/miniconda} folder.
#' @param add_makedb_options a character string specifying additional makedb options that shall be passed on to the diamond makedb command line call, e.g. \code{add_make_options = "--taxonnames"} (Default is \code{add_diamond_options = NULL}).
#' @param add_diamond_options a character string specifying additional diamond options that shall be passed on to the diamond command line call, e.g. \code{add_diamond_options = "--block-size 4.0 --compress 1 --no-self-hits"} (Default is \code{add_diamond_options = NULL}).
#' @author Hajk-Georg Drost
#' @examples
#' \dontrun{
#' # run diamond assuming that the diamond executable is available
#' # via the system path ('diamond_exec_path = NULL') and using
#' # sensitivity_mode = "ultra-sensitive"
#' diamond_example <- diamond_protein_against_ncbi_nr(
#'               query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
#'               ncbi_nr_folder_path_path = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
#'               sensitivity_mode = "ultra-sensitive",
#'               use_arrow_duckdb_connection  = FALSE)
#'
#' # look at DIAMOND results
#' diamond_example
#'
#' # run diamond assuming that the diamond executable is available
#' # via the miniconda path ('diamond_exec_path = "/opt/miniconda3/bin/"')
#' # and using 2 cores as well as sensitivity_mode = "ultra-sensitive"
#' diamond_example_conda <- diamond_protein_against_ncbi_nr(
#' query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
#' ncbi_nr_folder_path_path = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
#' sensitivity_mode = "ultra-sensitive", diamond_exec_path = "/opt/miniconda3/bin/",
#' use_arrow_duckdb_connection  = FALSE, cores = 2)
#'
#' # look at DIAMOND results
#' diamond_example_conda
#'
#' # run diamond assuming that the diamond executable is available
#' # via the system path ('diamond_exec_path = NULL') and using
#' # sensitivity_mode = "ultra-sensitive" and adding command line options:
#' # "--block-size 4.0 --compress 1 --no-self-hits"
#' diamond_example_ultra_sensitive_add_diamond_options <- diamond_protein_against_ncbi_nr(
#' query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
#' ncbi_nr_folder_path_path = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
#' sensitivity_mode = "ultra-sensitive",
#' max_target_seqs = 500,
#' use_arrow_duckdb_connection  = FALSE,
#' add_diamond_options = "--block-size 4.0 --compress 1 --no-self-hits",
#' cores = 1
#' )
#'
#' # look at DIAMOND results
#' diamond_example_ultra_sensitive_add_diamond_options
#'
#' # run diamond assuming that the diamond executable is available
#' # via the system path ('diamond_exec_path = NULL') and using
#' # sensitivity_mode = "ultra-sensitive" and adding makedb command line options:
#' # "--taxonnames"
#' diamond_example_ultra_sensitive_add_makedb_options <- diamond_protein_against_ncbi_nr(
#' query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
#' ncbi_nr_folder_path_path = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
#' sensitivity_mode = "ultra-sensitive",
#' max_target_seqs = 500,
#' use_arrow_duckdb_connection  = FALSE,
#' add_makedb_options = "--taxonnames",
#' cores = 1
#' )
#'
#' # look at DIAMOND results
#' diamond_example_ultra_sensitive_add_makedb_options
#' }
#'
#' @export

diamond_protein_against_ncbi_nr <- function(query,
                                            ncbi_nr_folder_path_path,
                                            make_diamond_db = FALSE,
                                            store_hit_table_in_tmp = TRUE,
                                            store_hit_table_at_path = tempdir(),
                                            task = "blastp",
                                            sensitivity_mode = "fast",
                                            use_arrow_duckdb_connection = TRUE,
                                            evalue   = 1E-20,
                                            out_format = "csv",
                                            cores = 1,
                                            max_target_seqs = "unlimited",
                                            hard_mask = TRUE,
                                            diamond_exec_path = NULL,
                                            add_makedb_options = NULL,
                                            add_diamond_options = NULL) {
  # check if diamond is properly installed on local machine
  if (!is_diamond_installed(diamond_exec_path))
    stop("Please install a valid version of DIAMOND2.", call. = FALSE)

  # check if ncbi nr database files are properly downloaded and stored
  # on local machine
  # if (!is_ncbi_nr_installed(ncbi_nr_folder_path))
  #   stop(
  #     "The files associated with the NCBI NR database could not be found. Please check your path to ncbi_nr_folder_path: ",
  #     ncbi_nr_folder_path,
  #     ".",
  #     call. = FALSE
  #   )

  if (!is.element(
    sensitivity_mode,
    c(
      "faster",
      "default",
      "fast",
      "mid-sensitive",
      "sensitive",
      "more-sensitive",
      "very-sensitive",
      "ultra-sensitive"
    )
  ))
    stop("Please specify a sensitivity_mode that is supported by this function.",
         call. = FALSE)

  if (store_hit_table_in_tmp &
      (store_hit_table_at_path != tempdir()))
    stop(
      "Please be aware that you specified 'store_hit_table_in_tmp = TRUE', meaning that the argument 'store_hit_table_at_path' has to specify the 'tempdir()' path on your system.",
      call. = FALSE
    )

  if (!store_hit_table_in_tmp &
      (store_hit_table_at_path == tempdir()))
    warning(
      "You specified 'store_hit_table_in_tmp = TRUE' but then provided them tempdir() folder path in 'store_hit_table_at_path: ",
      tempdir(),
      " as location to store the DIAMOND2 hit putput table.",
      "\n",
      "Please change the path in 'store_hit_table_at_path' if you wish to store your DIAMOND2 hit output table somewhere else, e.g. at location 'ncbi_nr_folder_path': ",
      ncbi_nr_folder_path,
      call. = FALSE
    )

  if (sensitivity_mode == "fast") {
    sensitivity_mode <- ""
  } else {
    sensitivity_mode <- paste0(" --", sensitivity_mode)
  }

  # make sure that input files contain the correct sequence type
  file_contains_aa(query, "query")
  #file_contains_aa(subject, "subject")

  # determine the number of cores on a multicore machine
  multi_cores <- parallel::detectCores()

  # in case one tries to use more cores than are available
  if (cores > multi_cores)
    stop("You chose more cores than are available on your machine.", call. = FALSE)

  # test if query file exists
  if (!file.exists(query))
    stop("Unfortunately, no query file has been found at ", query, call. = FALSE)

  # test if subject file or database exists
  if (!fs::dir_exists(ncbi_nr_folder_path_path) & !make_diamond_db)
    stop("Unfortunately, the folder has been found at ",
         ncbi_nr_folder_path_path,
         call. = FALSE)

  if (!is.element(task, c("blastp")))
    stop(
      "Please choose a protein-protein comparison task that is supported by DIAMOND2: task = 'blastp'",
      call. = FALSE
    )


  if (hard_mask) {
    hard_mask <- "1"
  } else {
    hard_mask <- "0"
  }

  if (!is.null(diamond_exec_path)) {
    message("\n")
    message(
      "Running diamond2 with '",
      diamond_exec_path,
      "/diamond blastp ",
      " with  query: ",
      query,
      " and subject: ",
      ifelse(
        fs::file_exists(file.path(
          ncbi_nr_folder_path_path, "nr_diamond.dmnd"
        )),
        file.path(ncbi_nr_folder_path_path, "nr_diamond.dmnd"),
        file.path(ncbi_nr_folder_path_path, "nr")
      ),
      " using ",
      cores,
      " core(s) ..."
    )
    message("\n")
    message(
      "Alignment Sensitivity Mode: ",
      ifelse(sensitivity_mode == "", "--fast", sensitivity_mode),
      " and max-target-seqs: ",
      max_target_seqs
    )
    message("Masking of low-complexity regions: ",
            ifelse(hard_mask == "1", "TANTAN", "none"))
    message("\n")
    message("\n")
    message(
      "The parameter 'make_diamond_db = ",
      make_diamond_db,
      "' was set. It is assumed that the subject NCBI NR database is already formatted either as *.dmnd or BLAST database file."
    )
  }

  if (is.null(diamond_exec_path)) {
    message("\n")
    message(
      "Starting 'diamond2 blastp ",
      " with  query: ",
      query,
      " and subject: ",
      ifelse(
        fs::file_exists(file.path(
          ncbi_nr_folder_path_path, "nr_diamond.dmnd"
        )),
        file.path(ncbi_nr_folder_path_path, "nr_diamond.dmnd"),
        file.path(ncbi_nr_folder_path_path, "nr")
      ),
      " using ",
      cores,
      " core(s) ..."
    )
    message("\n")
    message(
      "Alignment Sensitivity Mode: ",
      ifelse(sensitivity_mode == "", "--fast", sensitivity_mode),
      " and max-target-seqs: ",
      max_target_seqs
    )
    message("\n")
    message("\n")
    message("Masking of low-complexity regions: ",
            ifelse(hard_mask == "1", "TANTAN", "none"))
    message("\n")
  }


  if (!is.null(diamond_exec_path)) {
    diamond_call <-
      paste0(
        file.path(diamond_exec_path, "diamond"),
        " blastp",
        sensitivity_mode,
        " --query ",
        ws_wrap(query),
        " -o ",
        ws_wrap(file.path(
          ifelse(
            store_hit_table_in_tmp == tempdir(),
            tempdir(),
            store_hit_table_at_path
          ),
          paste0(basename(query),
                 "_diamond_hit_output.tsv")
        )),
        " --db ",
        ifelse(
          fs::file_exists(file.path(
            ncbi_nr_folder_path_path, "nr_diamond.dmnd"
          )),
          file.path(ncbi_nr_folder_path_path, "nr_diamond.dmnd"),
          file.path(ncbi_nr_folder_path_path, "nr")
        )
      )
  }

  if (is.null(diamond_exec_path)) {
    diamond_call <-
      paste0(
        "diamond blastp",
        sensitivity_mode,
        " --query ",
        ws_wrap(query),
        " -o ",
        ws_wrap(file.path(
          ifelse(
            store_hit_table_in_tmp == tempdir(),
            tempdir(),
            store_hit_table_at_path
          ),
          paste0(basename(query),
                 "_diamond_hit_output.tsv")
        )),
        " --db ",
        ifelse(
          fs::file_exists(file.path(
            ncbi_nr_folder_path_path, "nr_diamond.dmnd"
          )),
          file.path(ncbi_nr_folder_path_path, "nr_diamond.dmnd"),
          file.path(ncbi_nr_folder_path_path, "nr")
        )
      )
  }

  output_diamond <-
    file.path(
      ifelse(
        store_hit_table_in_tmp == tempdir(),
        tempdir(),
        store_hit_table_at_path
      ),
      paste0(
        stringr::str_replace(basename(query), "[.]fa", ""),
        "_diamond_hit_output.tsv"
      )
    )

  # output_blast without white space wrap - ws_wrap()
  output_read_diamond <- ws_wrap(file.path(
    ifelse(
      store_hit_table_in_tmp == tempdir(),
      tempdir(),
      store_hit_table_at_path
    ),
    paste0(
      stringr::str_replace(basename(query), "[.]fa", ""),
      "_diamond_hit_output.tsv"
    )
  ))

  tryCatch({
    # format subject into database if make_diamond_db == TRUE
    if (make_diamond_db) {
      message("\n")
      message(
        "The raw NCBI NR database stored at '",
        file.path(ncbi_nr_folder_path_path, "nr.gz"),
        "' will be formatted with 'diamond makedb'."
      )
      if (!is.null(diamond_exec_path)) {
        system(paste0(
          file.path(diamond_exec_path, "diamond"),
          " makedb --in ",
          file.path(ncbi_nr_folder_path_path, "nr.gz"),
          " --db ",
          file.path(ncbi_nr_folder_path_path, "nr_diamond.dmnd"),
          ifelse(
            !is.null(add_makedb_options),
            paste0(" ", add_makedb_options),
            ""
          )
        ))
      }

      if (is.null(diamond_exec_path)) {
        system(paste0(
          "diamond makedb --in ",
          ile.path(ncbi_nr_folder_path_path, "nr.gz"),
          " --db ",
          file.path(ncbi_nr_folder_path_path, "nr_diamond.dmnd"),
          ifelse(
            !is.null(add_makedb_options),
            paste0(" ", add_makedb_options),
            ""
          )
        ))
      }
    }
  }, error = function(e) {
    stop(
      "Something went wring when trying to run 'diamond makedb', please check your input data and whether the diamond executable is installed properly.",
      call. = FALSE
    )
  })

  tryCatch({
    system(
      paste0(
        # ifelse(
        #   is.null(diamond_exec_path),
        #   diamond_call,
        #   paste0("export PATH=$PATH:", diamond_exec_path)
        # ),
        diamond_call,
        " --evalue ",
        evalue,
        ifelse(max_target_seqs == "unlimited", " -k0"  , paste0(" -k ", max_target_seqs)),
        " -o ",
        output_diamond ,
        " --threads ",
        cores,
        " --masking ",
        hard_mask,
        ifelse(
          !is.null(add_diamond_options),
          paste0(" ", add_diamond_options),
          ""
        ),
        paste0(
          " --outfmt 6 qseqid sseqid pident nident length mismatch gapopen gaps positive ppos qstart qend qlen qcovhsp sstart send slen evalue bitscore score"
        )
      )
    )
  }, error = function(e) {
    stop(
      "Something went wring when trying to run 'diamond blastp', please check your input data and whether the diamond executable is installed properly.",
      call. = FALSE
    )
  })

  if (use_arrow_duckdb_connection) {
    diamond_tbl <- read_diamond(
      file = output_read_diamond,
      out_format = out_format,
      use_arrow_duckdb_connection = use_arrow_duckdb_connection
    )
    message("\n")
    message("DIAMOND2 search finished successfully!")
    message("\n")
    message("A Arrow->DuckDB database connection to the DIAMOND2 output file has been generated.")
    return(diamond_tbl)
  } else {
    diamond_tbl <- read_diamond(
      file = output_read_diamond,
      out_format = out_format,
      use_arrow_duckdb_connection = use_arrow_duckdb_connection
    )

    message("\n")
    message("DIAMOND2 search finished successfully! ")
    message("\n")
    message(
      "The DIAMOND2 output file was imported into the running R session. The DIAMOND2 output file has been stored at: ",
      output_diamond
    )
    return(diamond_tbl)
  }
}
