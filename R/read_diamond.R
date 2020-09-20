#' @title Import DIAMOND output file
#' @description When performing BLAST searches with the \code{blast_*()} functions,
#' the corresponding BLAST output file can be imported into the current R session using this function.
#' All output formats given by BLAST are supported (see e.g. description for details).
#' @param file path to BLAST output file.
#' @param out_format a character string specifying the output format of the BLAST output that shall be imported.
#' Available options are:
#'  \itemize{
#'  \item \code{out_format = "postgres"} : store BLAST output as Postgres database and generate postgres database connection.
#'  \item \code{out_format = "xml"} : XML
#'  \item \code{out_format = "csv"} : Comma-separated values
#'  }
#' @param postgres_user specify username for RPostgreSQL connection.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{diamond_protein_to_protein}}
#' @export

read_diamond <- function(file, out_format, postgres_user = NULL) {

  if (!file.exists(file))
    stop("The DIAMOND output file '", file, "' does not exist! Please check what might have went wrong with the DIAMOND call.", call. = FALSE)

  if (!is.element(
    out_format,
    c(
      "postgres",
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

   if (out_format == "postgres") {

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


       if (is.null(postgres_user))
           stop("Please specify a 'postgres_user' to import DIAMOND output into PostgresSQL database.", call. = FALSE)

       #require(RPostgreSQL)

       postgres_filename <- paste0(unlist(stringr::str_split(basename(file),"[.]"))[1],"_postgres")

       connect_db <-
           DBI::dbConnect(
               DBI::dbDriver("PostgreSQL"),
               user = postgres_user,
               password = "",
               host = "localhost",
               port = 5432,
               dbname = postgres_user)

       DBI::dbWriteTable(
           connect_db,
           name      = postgres_filename,
           value     = file,
           row.names = FALSE,
           header    = FALSE,
           sep       = "\t",
           overwrite = TRUE
       )

       blast_sql_db <-
           dplyr::src_postgres(
               dbname = postgres_user,
               host = "localhost",
               port = 5432,
               user = postgres_user,
               password = ""
           )

       blast_postgres <-
           dplyr::tbl(blast_sql_db, postgres_filename)

       on.exit({
           #sparklyr::spark_disconnect(sparkconnect)
           DBI::dbDisconnect(connect_db)

       })

       return(blast_postgres)
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
          message("Unfortunately, no DIAMOND hit was found for '", file, "'.")
          return(FALSE)
        }

    }

}
