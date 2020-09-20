#' @title Perform Protein to Protein DIAMOND Searches (BLASTP)
#' @description Run protein to protein DIAMOND of reference sequences
#' against a blast-able database or fasta file.
#' @param query path to input file in fasta format.
#' @param subject path to subject file in fasta format or blast-able database.
#' @param output_path path to folder at which DIAMOND output table shall be stored.
#' Default is \code{output_path = NULL} (hence \code{getwd()} is used).
#' @param is_subject_db logical specifying whether or not the \code{subject} file is a file in fasta format (\code{is_subject_db = FALSE}; default)
#' or a \code{fasta} file that was previously converted into a blast-able database using \code{diamond makedb} (\code{is_subject_db = TRUE}).
#' @param task protein search task option. Options are:
#' \itemize{
#' \item \code{task = "blastp"} : Standard protein-protein comparisons (default).
#' }
#' @param sensitivity_mode specify the level of alignment sensitivity. The higher the sensitivity level, the more deep homologs can be found, but at the cost of reduced computational speed.
#' \itemize{
#'   \item \code{sensitivity_mode = "fast"} : fastest alignment mode, but least sensitive (default).
#'   \item \code{sensitivity_mode = "more-sensitive"} : second fastest alignment mode, but second least sensitive.
#'   \item \code{sensitivity_mode = "very-sensitive"} : sensitive alignment mode.
#'   \item \code{sensitivity_mode = "ultra-sensitive"} : most sensitive alignment mode (sensitivity as high as BLASTP).
#' }
#' @param db_import shall the DIAMOND hit output be stored in a PostgresSQL database and shall a connection be established to this database? Default is \code{db_import = FALSE}.
#' In case users wish to to only generate a DIAMOND output file without importing it to the current R session they can specify \code{db_import = NULL}.
#' @param postgres_user when \code{db_import = TRUE} and \code{out_format = "postgres"} is selected, the DIAMOND output is imported and stored in a
#' PostgresSQL database. In that case, users need to have PostgresSQL installed and initialized on their system.
#' Please consult the Installation Vignette for details.
#' @param evalue Expectation value (E) threshold for saving hits (default: \code{evalue = 0.001}).
#' @param out_format a character string specifying the format of the file in which the DIAMOND results shall be stored.
#' Available options are:
#'  \itemize{
#'  \item \code{out_format = "pair"} : Pairwise
#'  \item \code{out_format = "xml"} : XML
#'  \item \code{out_format = "csv"} : Comma-separated file
#'  }
#' @param cores number of cores for parallel BLAST searches.
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
#' diamond_example <- diamond_protein_to_protein(
#'               query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
#'               subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
#'               sensitivity_mode = "ultra-sensitive",
#'               output_path = tempdir(),
#'               db_import  = FALSE)
#'
#' # look at DIAMOND results
#' diamond_example
#'
#' # run diamond assuming that the diamond executable is available
#' # via the miniconda path ('diamond_exec_path = "/opt/miniconda3/bin/"')
#' # and using 2 cores as well as sensitivity_mode = "ultra-sensitive"
#' diamond_example_conda <- diamond_protein_to_protein(
#' query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
#' subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
#' sensitivity_mode = "ultra-sensitive", diamond_exec_path = "/opt/miniconda3/bin/",
#' output_path = tempdir(),
#' db_import  = FALSE, cores = 2)
#'
#' # look at DIAMOND results
#' diamond_example_conda
#'
#' # run diamond assuming that the diamond executable is available
#' # via the system path ('diamond_exec_path = NULL') and using
#' # sensitivity_mode = "ultra-sensitive" and adding command line options:
#' # "--block-size 4.0 --compress 1 --no-self-hits"
#' diamond_example_ultra_sensitive_add_diamond_options <- diamond_protein_to_protein(
#' query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
#' subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
#' sensitivity_mode = "ultra-sensitive",
#' max_target_seqs = 500,
#' output_path = tempdir(),
#' db_import  = FALSE,
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
#' diamond_example_ultra_sensitive_add_makedb_options <- diamond_protein_to_protein(
#' query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
#' subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
#' sensitivity_mode = "ultra-sensitive",
#' max_target_seqs = 500,
#' output_path = tempdir(),
#' db_import  = FALSE,
#' add_makedb_options = "--taxonnames",
#' cores = 1
#' )
#'
#' # look at DIAMOND results
#' diamond_example_ultra_sensitive_add_makedb_options
#' }
#'
#' @export

diamond_protein_to_protein <- function(query,
                                       subject,
                                       output_path = NULL,
                                       is_subject_db = FALSE,
                                       task = "blastp",
                                       sensitivity_mode = "ultra-sensitive",
                                       db_import = FALSE,
                                       postgres_user = NULL,
                                       evalue   = 1E-3,
                                       out_format = "csv",
                                       cores = 1,
                                       max_target_seqs = 500,
                                       hard_mask = TRUE,
                                       diamond_exec_path = NULL,
                                       add_makedb_options = NULL,
                                       add_diamond_options = NULL) {

  if (!is_diamond_installed(diamond_exec_path))
    stop("Please install a valid version of DIAMOND.", call. = FALSE)

  if (!is.null(db_import)) {
    if (db_import) {
      if (!is.element(out_format, c("xml", "tab", "csv")))
        stop("Only output formats: 'xml', 'tab', or 'csv' can be imported so far.",
             call. = FALSE)
    }
  }

  if (!is.element(
    sensitivity_mode,
    c(
      "fast",
      "more-sensitive",
      "very-sensitive",
      "ultra-sensitive"
    )
  ))
    stop("Please specify a sensitivity_mode that is supported by this function.",
         call. = FALSE)

  if (sensitivity_mode == "fast") {
    sensitivity_mode <- ""
  } else {
    sensitivity_mode <- paste0(" --", sensitivity_mode)
  }

  # make sure that input files contain the correct sequence type
  file_contains_aa(query, "query")
  file_contains_aa(subject, "subject")

  # determine the number of cores on a multicore machine
  multi_cores <- parallel::detectCores()

  # in case one tries to use more cores than are available
  if (cores > multi_cores)
    stop("You chose more cores than are available on your machine.", call. = FALSE)

  # test if query file exists
  if (!file.exists(query))
    stop("Unfortunately, no query file has been found at ", query, call. = FALSE)

  # test if subject file or database exists
  if (!file.exists(subject) & !is_subject_db)
    stop("Unfortunately, no subject file has been found at ",
         subject,
         call. = FALSE)

  if (!is.element(task, c("blastp")))
    stop(
      "Please choose a protein-protein comparison task that is supported by DIAMOND: task = 'blastp'",
      call. = FALSE
    )

  if (!is.null(output_path)){
    if (!file.exists(output_path))
      stop("Please specify a valid output path. Your specified path does not seem to exist: ", output_path, call. = FALSE)
  }

  if (hard_mask) {
    hard_mask <- "1"
  } else {
    hard_mask <- "0"
  }

  if (!is.null(diamond_exec_path)) {
    message("\n")
    message(
      "Running diamond with '",
      diamond_exec_path,
      "/diamond blastp ",
      " with  query: ",
      query,
      " and subject: ",
      subject,
      " using ",
      cores,
      " core(s) ..."
    )
    message("\n")
    message("Alignment Sensitivity Mode: ", ifelse(sensitivity_mode == "", "--fast", sensitivity_mode), " and max-target-seqs: ", max_target_seqs)
    message("\n")
    message("Masking of low-complexity regions: ", ifelse(hard_mask == "1", "TANTAN", "none"))
    message("\n")
  }

  if (is.null(diamond_exec_path)) {
    message("\n")
    message(
      "Starting 'diamond blastp ",
      " with  query: ",
      query,
      " and subject: ",
      subject,
      " using ",
      cores,
      " core(s) ..."
    )
    message("\n")
    message("Alignment Sensitivity Mode: ", ifelse(sensitivity_mode == "", "--fast", sensitivity_mode), " and max-target-seqs: ", max_target_seqs)
    message("\n")
    message("\n")
    message("Masking of low-complexity regions: ", ifelse(hard_mask == "1", "TANTAN", "none"))
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
        " --db ",
        ws_wrap(file.path(tempdir(), basename(subject)))
      )
  }

  if (is.null(diamond_exec_path)) {
    diamond_call <-
      paste0(
        "diamond blastp",
        sensitivity_mode,
        " --query ",
        ws_wrap(query),
        " --db ",
        ws_wrap(file.path(tempdir(), basename(subject)))
      )
  }

  output_diamond <-
    file.path(
      ifelse(
        is.null(output_path),
        ws_wrap(getwd()),
        ws_wrap(output_path)
      ),
      paste0(
        unlist(stringr::str_split(basename(query), "[.]"))[1],
        "_",
        unlist(stringr::str_split(basename(subject), "[.]"))[1],
        "_",
        task,
        "_eval_",
        evalue,
        ".blast_tbl"
      )
    )

  # output_blast without white space wrap - ws_wrap()
  output_read_diamond <-
    file.path(
      ifelse(is.null(output_path), getwd(), output_path),
      paste0(
        unlist(stringr::str_split(basename(query), "[.]"))[1],
        "_",
        unlist(stringr::str_split(basename(subject), "[.]"))[1],
        "_",
        task,
        "_eval_",
        evalue,
        ".blast_tbl"
      )
    )


tryCatch({

  # format subject into database if is_subject_db == FALSE
  if (!is_subject_db) {
    if (!is.null(diamond_exec_path)) {
      system(
        paste0(
          file.path(diamond_exec_path, "diamond"), " makedb --in ",
          subject,
          " --db ",
          file.path(tempdir(), basename(subject)),
          ifelse(!is.null(add_makedb_options), paste0(" ", add_makedb_options), "")
        )
      )
    }

    if (is.null(diamond_exec_path)) {
      system(
        paste0(
          "diamond makedb --in ",
          subject,
          " --db ",
          file.path(tempdir(), basename(subject)),
          ifelse(!is.null(add_makedb_options), paste0(" ", add_makedb_options), "")
        )
      )
    }
  }
}, error = function(e) {stop("Something went wring when trying to run 'diamond makedb', please check your input data and whether the diamond executable is installed properly.",
                          call. = FALSE)})

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
      " -k ",
      max_target_seqs,
      " -o ",
      output_diamond ,
      " --threads ",
      cores,
      " --masking ",
      hard_mask,
      ifelse(!is.null(add_diamond_options), paste0(" ", add_diamond_options), ""),
      paste0(
        " --outfmt 6 qseqid sseqid pident nident length mismatch gapopen gaps positive ppos qstart qend qlen qcovhsp sstart send slen evalue bitscore score"
      )
    )
  )}, error = function(e) {stop("Something went wring when trying to run 'diamond blastp', please check your input data and whether the diamond executable is installed properly.",
                      call. = FALSE)})

  if (!is.null(db_import)) {
    if (db_import) {
      diamond_tbl <- read_diamond(file = output_read_diamond,
                              out_format = "postgres",
                              postgres_user = postgres_user)
      message("\n")
      message(
        "DIAMOND search finished successfully!"
      )
      message("\n")
      message("A Postgres database connection to the DIAMOND output file has been generated. The DIAMOND output file can be found at: ",
              output_blast)
      return(diamond_tbl)
    } else {
      diamond_tbl <- read_diamond(file = output_read_diamond,
                              out_format = out_format,
                              postgres_user = NULL)

      message("\n")
      message(
        "DIAMOND search finished successfully! "
      )
      message("\n")
      message("The DIAMOND output file was imported into the running R session. The DIAMOND output file has been stored at: ",
              output_diamond)
      return(diamond_tbl)
    }
  } else {
    message("\n")
    message("DIAMOND search finished! DIAMOND output file has been stored at: ",
            output_diamond)
  }
}
