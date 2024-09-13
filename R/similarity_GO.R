
env = new.env()

#' Calculate Gene Ontology (GO) semantic similarity matrix
#'
#' @param go_id A vector of GO IDs.
#' @param ont Sub-ontology of GO. Value should be one of "BP", "CC" or "MF". If it is not specified,
#'      the function automatically identifies it by random sampling 10 IDs from `go_id` (see `guess_ont()`).
#' @param db Annotation database. It should be an OrgDb package name from \url{https://bioconductor.org/packages/release/BiocViews.html#___OrgDb}. The value
#'    can also directly be an `OrgDb` object.
#' @param measure Semantic measure for the GO similarity, pass to [`simona::term_sim()`]. All valid values are in [`simona::all_term_sim_methods()`].
#'
#' @details
#' The default similarity method is "Sim_XGraSM_2013". Since the semantic similarities are calculated based on gene annotations to GO terms,
#' I suggest users also try the following methods:
#' 
#' - `"Sim_Lin_1998"`
#' - `"Sim_Resnik_1999"`
#' - `"Sim_Relevance_2006"`
#' - `"Sim_SimIC_2010"`
#' - `"Sim_XGraSM_2013"`
#' - `"Sim_EISI_2015"`
#' - `"Sim_AIC_2014"`
#' - `"Sim_Wang_2007"`
#' - `"Sim_GOGO_2018"`
#'
#' @return
#' `GO_similarity()` returns a symmetric matrix.
#' @export
#' @import simona
#' @import GetoptLong
#' @examples
#' \donttest{
#' go_id = random_GO(100)
#' mat = GO_similarity(go_id)
#' }
GO_similarity = function(go_id, ont = NULL, db = "org.Hs.eg.db", measure = "Sim_XGraSM_2013") {

	if(is.null(ont)) {
		ont = guess_ont(go_id, db)
		if(is.null(ont)) {
			stop_wrap("Cannot determine which GO ontology (BP/CC/MF) you are using. Please manualy set `ont` argument.")
		}
		message(qq("You haven't provided value for `ont`, guess it as `@{ont}`."))
	}

	hash = digest(list(ont = ont, db = db))
	if(is.null(env$go[[hash]])) {
		dag = create_ontology_DAG_from_GO_db(namespace = ont, org_db = db, relations = c("part_of", "regulates"))

		ic = term_IC(dag, method = "IC_annotation")
		all_go_id = names(ic[!is.na(ic)])

		env$go[[hash]] = list(dag = dag, all_go_id = all_go_id)
	} else {
		dag = env$go[[hash]]$dag
		all_go_id = env$go[[hash]]$all_go_id
	}

	go_sim = term_sim(dag, go_id, method = measure)

	attr(go_sim, "measure") = measure
	attr(go_sim, "ontology") = paste0("GO:", ont)
	return(go_sim)
}


split_by_block = function(n, size) {
	size = min(c(n, size))
	REST = n %% size
    LARGE = n - REST
    NBLOCKS = n %/% size
    GROUP = rep(1:NBLOCKS, each = size)
    if (REST > 0) GROUP = c(GROUP, rep(NBLOCKS + 1, REST))
    split(1:n, GROUP)
}


#' @rdname GO_similarity
#'
#' @details
#' In `guess_ont()`, only 10 random GO IDs are checked.
#'
#' @return
#' `guess_ont()` returns a single character scalar of "BP", "CC" or "MF". 
#' If there are more than one ontologies detected. It returns `NULL`.
#' 
#' @export
#' @import AnnotationDbi
#' @examples
#' \donttest{
#' go_id = random_GO(100)
#' guess_ont(go_id)
#' }
guess_ont = function(go_id, db = 'org.Hs.eg.db') {

	if(is.character(db)) {
		db = get(db, asNamespace(db))
	}
	test_go_id = sample(go_id, min(c(length(go_id), 10)))
	suppressMessages(df <- select(db, keys = test_go_id, columns = "ONTOLOGY", keytype = "GO"))
	guess_ont = unique(df$ONTOLOGY)
	guess_ont = guess_ont[!is.na(guess_ont)]
	if(length(guess_ont) != 1) {
		return(NULL)
	} else {
		return(guess_ont)
	}
}

#' @rdname GO_similarity
#'
#' @param n Number of GO IDs.
#' 
#' @details
#' In `random_GO()`, only GO terms with gene annotations are sampled.
#'
#' @return
#' `random_GO()` returns a vector of GO IDs.
#' @export
random_GO = function(n, ont = c("BP", "CC", "MF"), db = "org.Hs.eg.db") {
	ont = match.arg(ont)
	hash = digest(list(ont = ont, db = db))
	if(is.null(env$go[[hash]])) {
		dag = create_ontology_DAG_from_GO_db(namespace = ont, org_db = db, relations = c("part_of", "regulates"))

		ic = term_IC(dag, method = "IC_annotation")
		all_go_id = names(ic[!is.na(ic)])

		env$go[[hash]] = list(dag = dag, all_go_id = all_go_id)
	} else {
		dag = env$go[[hash]]$dag
		all_go_id = env$go[[hash]]$all_go_id
	}

	sample(all_go_id, min(n, length(all_go_id)))
}
