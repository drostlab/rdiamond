#' @title Perform a DIAMOND Protein to Protein best hit search
#' @description This function performs a DIAMOND search (best hit) of a given set of protein sequences against a given database.
#' @param query a character string specifying the path to the protein sequence file of interest (query organism).
#' @param subject a character string specifying the path to the protein sequence file of interest (subject organism).
#' @param is_subject_db logical specifying whether or not the \code{subject} file is a file in \code{fasta} format (\code{is_subject_db = FALSE}; default)
#' or a \code{fasta} file that was previously converted into a blast-able database using \code{diamond makedb} (\code{is_subject_db = TRUE}).
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format} = \code{"fasta"}.
#' Default is \code{format} = \code{"fasta"}.
#' @param sensitivity_mode specify the level of alignment sensitivity. The higher the sensitivity level, the more deep homologs can be found, but at the cost of reduced computational speed.
#' \itemize{
#'   \item \code{sensitivity_mode = "faster"} : fastest alignment mode, but least sensitive (default). Designed for finding hits of >70% identity and short read alignment such as NGS reads.
#'   \item \code{sensitivity_mode = "default"} : Default mode. Designed for finding hits of >70% identity and short read alignment such as NGS reads.
#'   \item \code{sensitivity_mode = "fast"} : fastest alignment mode, but least sensitive (default). Designed for finding hits of >70% identity and short read alignment such as NGS reads.
#'   \item \code{sensitivity_mode = "mid-sensitive"} : fast alignments between the \code{fast} mode and the sensitive mode in sensitivity.
#'   \item \code{sensitivity_mode = "sensitive"} : fast alignments, but full sensitivity for hits >40% identity.
#'   \item \code{sensitivity_mode = "more-sensitive"} : more sensitive than the \code{sensitive} mode.
#'   \item \code{sensitivity_mode = "very-sensitive"} : sensitive alignment mode.
#'   \item \code{sensitivity_mode = "ultra-sensitive"} : most sensitive alignment mode (sensitivity as high as BLASTP).
#' }
#' @param out_format a character string specifying the format of the file in which the DIAMOND results shall be stored.
#' Available options are:
#'  \itemize{
#'  \item \code{out_format = "pair"} : Pairwise
#'  \item \code{out_format = "xml"} : XML
#'  \item \code{out_format = "csv"} : Comma-separated file
#'  }
#' @param evalue Expectation value (E) threshold for saving hits (default: \code{evalue = 0.001}).
#' @param max_target_seqs maximum number of aligned sequences that shall be retained. Please be aware that \code{max_target_seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param cores number of cores for parallel DIAMOND searches.
#' @param hard_mask shall low complexity regions be hard masked with TANTAN? Default is \code{db_hard_mask = TRUE}.
#' @param diamond_exec_path a path to the DIAMOND executable or \code{conda/miniconda} folder.
#' @param add_makedb_options a character string specifying additional makedb options that shall be passed on to the diamond makedb command line call, e.g. \code{add_make_options = "--taxonnames"} (Default is \code{add_diamond_options = NULL}).
#' @param add_diamond_options a character string specifying additional diamond options that shall be passed on to the diamond command line call, e.g. \code{add_diamond_options = "--block-size 4.0 --compress 1 --no-self-hits"} (Default is \code{add_diamond_options = NULL}).
#' @param output_path a path to the location were the DIAMOND best hit output shall be stored. E.g. \code{output_path} = \code{getwd()}
#' to store it in the current working directory, or \code{output_path} = \code{file.path("put", "your", "path", "here")}.
#' @author Hajk-Georg Drost
#' @details Given a set of protein sequences (query sequences), a best hit diamond search (DBH) is being performed.
#' @examples \dontrun{
#' # performing homology inference using the diamond best hit (DBH) method using protein sequences
#' best_hits <- diamond_protein_to_protein_best_hits(
#' query   = system.file('seqs/ortho_thal_aa.fasta', package = 'rdiamond'),
#' subject = system.file('seqs/ortho_lyra_aa.fasta', package = 'rdiamond'))
#' # look at results
#' best_hits
#'
#' # store the DIAMOND output file to the current working directory
#' best_hits <- diamond_protein_to_protein_best_hits(
#' query  = system.file('seqs/ortho_thal_aa.fasta', package = 'rdiamond'),
#' subject = system.file('seqs/ortho_lyra_aa.fasta', package = 'rdiamond'),
#' output_path  = getwd())
#' # look at results
#' best_hits
#'
#' # run diamond_best_hits() with multiple cores
#' best_hits <- diamond_protein_to_protein_best_hits(
#' query   = system.file('seqs/ortho_thal_aa.fasta', package = 'rdiamond'),
#' subject = system.file('seqs/ortho_lyra_aa.fasta', package = 'rdiamond'),
#' cores   = 2)
#' # look at results
#' best_hits
#'
#' # performing homology inference using the diamond best hit (DBH) method and
#' # specifying the path to the DIAMOND executable (here miniconda path)
#' best_hits <- diamond_protein_to_protein_best_hits(
#' query  = system.file('seqs/ortho_thal_aa.fasta', package = 'rdiamond'),
#' subject = system.file('seqs/ortho_lyra_aa.fasta', package = 'rdiamond'),
#' diamond_exec_path = "/opt/miniconda3/bin/")
#' # look at results
#' best_hits
#' }
#'
#' @return A tibble as returned by the \code{diamond_best_hits} function, storing the \code{query_ids}
#' in the first column and the \code{subject_ids} (best hit homologs) in the second column.
#' @seealso \code{\link{diamond_protein_to_protein_best_reciprocal_hits}}, \code{\link{diamond_protein_to_protein}}
#' @export

diamond_protein_to_protein_best_hits <- function(
    query,
    subject,
    is_subject_db = FALSE,
    format          = "fasta",
    sensitivity_mode = "ultra-sensitive",
    out_format       = "csv",
    evalue            = "1E-5",
    max_target_seqs   = 5000,
    cores           = 1,
    hard_mask       = TRUE,
    diamond_exec_path = NULL,
    add_makedb_options = NULL,
    add_diamond_options = NULL,
    output_path     = NULL){
  
  
  if (!is.element(format, c("fasta")))
    stop("Please specify a query and subject format that is supported by this function: format = 'fasta'.", call. = FALSE)
  
  # performing a DIAMOND search from query against subject: diamond_protein_to_protein(query,subject)
  hit_tbl <- rdiamond::diamond_protein_to_protein(
    query      = query,
    subject    = subject,
    is_subject_db = is_subject_db,
    task = "blastp",
    sensitivity_mode = sensitivity_mode,
    out_format = out_format,
    evalue  = evalue,
    max_target_seqs = max_target_seqs,
    hard_mask = hard_mask,
    diamond_exec_path = diamond_exec_path,
    add_makedb_options = add_makedb_options,
    add_diamond_options = add_diamond_options,
    output_path     = output_path)
  
  query_id <- '.' <- NULL
  hit_tbl <- dplyr::do(dplyr::group_by(hit_tbl, query_id), filter_best_hits(.))
  return(hit_tbl)
}