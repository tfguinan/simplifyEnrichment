

#' Global parameters
#'
#' @param ... Arguments for the parameters, see "details" section.
#' @param RESET Whether to reset to default values.
#' @param READ.ONLY Please ignore.
#' @param LOCAL Please ignore.
#' @param ADD Please ignore.
#' 
#' @details
#' There are the following global options:
#'
#' - `verobse`: Whether to print messages.
#' 
#' @return A `GlobalOptionsFun` object.
#' @export
#' @import GlobalOptions
se_opt = setGlobalOptions(
	verbose = TRUE
)
