#' @title Import DIAMOND2 output file
#' @description When performing BLAST searches with the \code{blast_*()} functions,
#' the corresponding BLAST output file can be imported into the current R session using this function.
#' All output formats given by BLAST are supported (see e.g. description for details).
#' @param file path to BLAST output file.
#' @param out_format a character string specifying the output format of the BLAST output that shall be imported.
#' Available options are:
#'  \itemize{
#'  \item \code{out_format = "xml"} : XML
#'  \item \code{out_format = "csv"} : Comma-separated values
#'  }
#' @param use_arrow_duckdb_connection shall DIAMOND2 hit output table be transformed to an in-process (big data disk-processing) arrow connection to DuckDB? This is useful when the DIAMOND2 output table to too large to fit into memory. Default is \code{use_arrow_duckdb_connection = FALSE}. 
#' @author Hajk-Georg Drost
#' @seealso \code{\link{diamond_protein_to_protein}}
#' @export

read_diamond <- function(file, out_format, use_arrow_duckdb_connection = FALSE) {

  if (!file.exists(file))
    stop("The DIAMOND2 output file '", file, "' does not exist! Please check what might have went wrong with the DIAMOND2 call.", call. = FALSE)

  if (!is.element(
    out_format,
    c(
      "xml",
      "tsv",
      "csv",
      "json.seq.aln",
      "json.blast.multi",
      "xml2.blast.multi",
      "json.blast.single",
      "xml2.blast.single"
    )
  )
  )
  stop("Sorry, but '", out_format,"' is not an available import type. Please choose an 'out_format' that is supported by this function.", call. = FALSE)

   if ((out_format == "csv" ) & use_arrow_duckdb_connection) {

       # # install local version of Spark if not available yet
       # if (nrow(sparklyr::spark_installed_versions()) == 0) {
       #     sparklyr::spark_install(version = spark_version)
       # }
       #
       # # open local connection to spark
       # sparkconnect <- sparklyr::spark_connect(master = "local")
       #
       # import_blast_tbl <- sparklyr::copy_to(sparkconnect, iris,
       #                                       "spark_blast_tbl",
       #                                        overwrite = TRUE)
     
       #require(RPostgreSQL)

       # postgres_filename <- paste0(unlist(stringr::str_split(basename(file),"[.]"))[1],"_postgres")
       # 
       # connect_db <-
       #     DBI::dbConnect(
       #         DBI::dbDriver("PostgreSQL"),
       #         user = postgres_user,
       #         password = "",
       #         host = "localhost",
       #         port = 5432,
       #         dbname = postgres_user)
       # 
       # DBI::dbWriteTable(
       #     connect_db,
       #     name      = postgres_filename,
       #     value     = file,
       #     row.names = FALSE,
       #     header    = FALSE,
       #     sep       = "\t",
       #     overwrite = TRUE
       # )
       # 
       # blast_sql_db <-
       #     dplyr::src_postgres(
       #         dbname = postgres_user,
       #         host = "localhost",
       #         port = 5432,
       #         user = postgres_user,
       #         password = ""
       #     )
       # 
       # blast_postgres <-
       #     dplyr::tbl(blast_sql_db, postgres_filename)
       # 
       # on.exit({
       #     #sparklyr::spark_disconnect(sparkconnect)
       #     DBI::dbDisconnect(connect_db)
       # 
       # })
     # query_id -> subject_id -> perc_identity -> num_ident_matches -> alig_length -> mismatches -> gap_openings -> n_gaps -> pos_match -> ppos -> q_start -> q_end -> q_len ->
     #   qcovhsp -> s_start -> s_end -> s_len -> evalue -> bit_score -> score_raw -> NULL
     
     diamond_duckdb <-
       arrow::open_tsv_dataset(
         file,
         col_names = FALSE,
         col_types = arrow::schema(
           arrow::field("X1", arrow::string()),
           arrow::field("X2", arrow::string()),
           arrow::field("X3", arrow::float64()),
           arrow::field("X4", arrow::int32()),
           arrow::field("X5", arrow::int32()),
           arrow::field("X6", arrow::int32()),
           arrow::field("X7", arrow::int32()),
           arrow::field("X8", arrow::int32()),
           arrow::field("X9", arrow::int32()),
           arrow::field("X10", arrow::float64()),
           arrow::field("X11", arrow::int32()),
           arrow::field("X12", arrow::int32()),
           arrow::field("X13", arrow::int32()),
           arrow::field("X14", arrow::float64()),
           arrow::field("X15", arrow::int32()),
           arrow::field("X16", arrow::int32()),
           arrow::field("X17", arrow::int32()),
           arrow::field("X18", arrow::float64()),
           arrow::field("X19", arrow::int32()),
           arrow::field("X20", arrow::float64())
         )
       ) |>  arrow::to_duckdb() |> dplyr::rename(
         query_id  = f0,
         subject_id  = f1,
         perc_identity  = f2,
         num_ident_matches  = f3,
         alig_length  = f4,
         mismatches  = f5,
         gap_openings  = f6,
         n_gaps  = f7,
         pos_match  = f8,
         ppos  = f9,
         q_start  = f10,
         q_end  = f11,
         q_len  = f12,
         qcovhsp  = f13,
         s_start  = f14,
         s_end  = f15,
         s_len  = f16,
         evalue  = f17,
         bit_score  = f18,
         score_raw  = f19
       )
     
       return(diamond_duckdb)
   }


    if (out_format == "csv") {
        diamond_csv <- readr::read_delim(file = file, delim = "\t",
                                       col_names = FALSE,
                                       col_types = readr::cols(
                                           "X1" = readr::col_character(),
                                           "X2" = readr::col_character(),
                                           "X3" = readr::col_double(),
                                           "X4" = readr::col_integer(),
                                           "X5" = readr::col_integer(),
                                           "X6" = readr::col_integer(),
                                           "X7" = readr::col_integer(),
                                           "X8" = readr::col_integer(),
                                           "X9" = readr::col_integer(),
                                           "X10" = readr::col_double(),
                                           "X11" = readr::col_integer(),
                                           "X12" = readr::col_integer(),
                                           "X13" = readr::col_integer(),
                                           "X14" = readr::col_double(),
                                           "X15" = readr::col_integer(),
                                           "X16" = readr::col_integer(),
                                           "X17" = readr::col_integer(),
                                           "X18" = readr::col_double(),
                                           "X19" = readr::col_number(),
                                           "X20" = readr::col_double()
                                       ))
        if (nrow(diamond_csv) > 0) {
          if (ncol(diamond_csv) != length(blast_outfmt_colnames()))
            stop("Tne number of diamond output columns and the number of column names does not match! Please check what might have gone wrong.", call. = FALSE)
          colnames(diamond_csv) <- blast_outfmt_colnames()
          return(diamond_csv)
        } else {
          message("Unfortunately, no DIAMOND2 hit was found for '", file, "'.")
          return(FALSE)
        }

    }

}
