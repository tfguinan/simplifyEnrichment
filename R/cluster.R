
#' Cluster terms based on their similarity matrix
#'
#' @param mat A similarity matrix.
#' @param method The clustering methods. Value should be in [`all_clustering_methods()`].
#' @param control A list of parameters passed to the corresponding clustering function.
#' @param verbose Whether to print messages.
#'
#' @details
#' New clustering methods can be registered by [`register_clustering_methods()`].
#'
#' Please note it is better to directly use `cluster_terms()` for clustering while not the individual ``cluster_by_*`` functions
#' because `cluster_terms()` does additional cluster label adjustment.
#' 
#' By default, there are the following clustering methods and corresponding clustering functions:
#' 
#' - `kmeans` see [`cluster_by_kmeans()`].
#' - `dynamicTreeCut` see [`cluster_by_dynamicTreeCut()`].
#' - `mclust` see [`cluster_by_mclust()`].
#' - `apcluster` see [`cluster_by_apcluster()`].
#' - `hdbscan` see [`cluster_by_hdbscan()`].
#' - `fast_greedy` see [`cluster_by_fast_greedy()`].
#' - `louvain` see [`cluster_by_louvain()`].
#' - `walktrap` see [`cluster_by_walktrap()`].
#' - `MCL` see [`cluster_by_MCL()`].
#' - `binary_cut` see [`binary_cut()`].
#' 
#' The additional argument in individual clustering functions can be set with the `control` argument
#' in `cluster_terms()`.
#'
#' @returns
#' A vector of numeric cluster labels.
#' @export
#' @rdname cluster_terms
cluster_terms = function(mat, method = "binary_cut", control = list(), verbose = se_opt$verbose) {
	
	if(nrow(mat) != ncol(mat)) {
		stop_wrap("The matrix should be square.")
	}

	if(verbose) message(qq("Cluster @{nrow(mat)} terms by '@{method}'..."), appendLF = FALSE)
	flush.console()

	fun = get_clustering_method(method, control = control)
	
	t1 = Sys.time()
	cl = fun(mat)
	t2 = Sys.time()

	t_diff = t2 - t1
	t_diff = format(t_diff)
	if(verbose) message(qq(" @{length(unique(cl))} clusters, used @{t_diff}."))

	if(length(unique(cl)) > 1) {
		if(method != "binary_cut") {
			#reorder the class labels
			class_mean = tapply(1:nrow(mat), cl, function(ind) {
				colMeans(mat[ind, , drop = FALSE])
			})
			class_mean = do.call(rbind, class_mean)
			ns = nrow(class_mean)

			class_dend = as.dendrogram(hclust(stats::dist(class_mean)))
			class_dend = reorder(class_dend, wts = rowSums(class_mean))
			map = structure(1:ns, names = order.dendrogram(class_dend))

			cl = map[as.character(cl)]
		}
	}

	attr(cl, "running_time") = t_diff

	return(cl)
}


#' @param max_k Maximal _k_ for k-means/PAM clustering. K-means/PAM clustering is applied from k = 2 to k = max_k.
#' @param ... Other arguments.
#'
#' @details
#' `cluster_by_kmeans()`: The best k for k-means clustering is determined according to the "elbow" or "knee" method on
#' the distribution of within-cluster sum of squares (WSS) on each k. All other arguments are passed
#' from `...` to [`stats::kmeans()`].
#'
#' @export
#' @rdname cluster_terms
cluster_by_kmeans = function(mat, max_k = max(2, min(round(nrow(mat)/5), 100)), ...) {
	
	cl = list()
	wss = NULL

	if(max_k <= 2) {
		stop_wrap("`max_k` should be larger than 2.")
	}

	for (i in 2:max_k) {
		suppressWarnings(km <- kmeans(mat, centers = i, iter.max = 50))
		cl[[i - 1]] = km$cluster
		wss[i - 1] = sum(km$withinss)
	}
	best_km = min(elbow_finder(2:max_k, wss)[1], knee_finder(2:max_k, wss)[1])

	cl[[best_km - 1]]
}


#' @details
#' `cluster_by_pam()`: PAM is applied by [`fpc::pamk()`] which can automatically select the best k.
#' All other arguments are passed from `...` to [`fpc::pamk()`].
#'
#' @export
#' @rdname cluster_terms
cluster_by_pam = function(mat, max_k = max(2, min(round(nrow(mat)/10), 100)), ...) {

	check_pkg("fpc", bioc = FALSE)

	best = fpc::pamk(mat, krange = 2:max_k, ...)
	best$pamobject$clustering
}

# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
elbow_finder = function(x_values, y_values) {
	# Max values to create line
	max_x_x = max(x_values)
	max_x_y = y_values[which.max(x_values)]
	max_y_y = max(y_values)
	max_y_x = x_values[which.max(y_values)]
	max_df = data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

	# Creating straight line between the max values
	fit = lm(max_df$y ~ max_df$x)

	# Distance from point to line
	distances = c()
	for(i in seq_along(x_values)) {
		distances = c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
	}

	# Max distance point
	x_max_dist = x_values[which.max(distances)]
	y_max_dist = y_values[which.max(distances)]

	return(c(x_max_dist, y_max_dist))
}

# https://raghavan.usc.edu//papers/kneedle-simplex11.pdf
knee_finder = function(x, y) {
	n = length(x)
	a = (y[n] - y[1])/(x[n] - x[1])
	b = y[1] - a*x[1]
	d = a*x - y
	x[which.max(d)]
}

#' @param minClusterSize Minimal number of objects in a cluster. Pass to [`dynamicTreeCut::cutreeDynamic()`].
#' 
#' @details
#' `cluster_by_dynamicTreeCut()`: All other arguments are passed from `...` to [`dynamicTreeCut::cutreeDynamic()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_dynamicTreeCut = function(mat, minClusterSize = 5, ...) {

	check_pkg("dynamicTreeCut", bioc = FALSE)

	cl = dynamicTreeCut::cutreeDynamic(hclust(stats::dist(mat)), distM = 1 - mat, minClusterSize = minClusterSize, verbose = 0, ...)
	unname(cl)
}


#' @details
#' `cluster_by_fast_greedy()`: All other arguments are passed from `...` to [`igraph::cluster_fast_greedy()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_fast_greedy = function(mat, ...) {
	check_pkg("igraph", bioc = FALSE)

	g = igraph::graph_from_adjacency_matrix(mat, mode = "upper", weighted = TRUE)
	imc = igraph::cluster_fast_greedy(g, ...)
	as.vector(igraph::membership(imc))
}

#' @details
#' `cluster_by_leading_eigen()`: All other arguments are passed from `...` to [`igraph::cluster_leading_eigen()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_leading_eigen = function(mat, ...) {
	check_pkg("igraph", bioc = FALSE)

	g = igraph::graph_from_adjacency_matrix(mat, mode = "upper", weighted = TRUE)
	imc = igraph::cluster_leading_eigen(g, ...)
	as.vector(igraph::membership(imc))
}

#' @details
#' `cluster_by_louvain()`: All other arguments are passed from `...` to [`igraph::cluster_louvain()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_louvain = function(mat, ...) {
	check_pkg("igraph", bioc = FALSE)

	g = igraph::graph_from_adjacency_matrix(mat, mode = "upper", weighted = TRUE)
	imc = igraph::cluster_louvain(g, ...)
	as.vector(igraph::membership(imc))
}

#' @details
#' `cluster_by_walktrap()`: All other arguments are passed from `...` to [`igraph::cluster_walktrap()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_walktrap = function(mat, ...) {
	check_pkg("igraph", bioc = FALSE)

	g = igraph::graph_from_adjacency_matrix(mat, mode = "upper", weighted = TRUE)
	imc = igraph::cluster_walktrap(g, ...)
	as.vector(igraph::membership(imc))
}

#' @param G Passed to the `G` argument in [`mclust::Mclust()`] which is the number of clusters.
#' @details
#' `cluster_by_mclust()`: All other arguments are passed from `...` to [`mclust::Mclust()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_mclust = function(mat, G = seq_len(max(2, min(round(nrow(mat)/5), 100))), ...) {
	
	check_pkg("mclust", bioc = FALSE)

	mclustBIC = mclust::mclustBIC

	fit = mclust::Mclust(mat, G = G, verbose = FALSE, control = mclust::emControl(itmax = c(1000, 1000)), ...)

	unname(fit$classification)
}

#' @param s Passed to the `s` argument in [`apcluster::apcluster()`].
#' @details
#' `cluster_by_apcluster()`: All other arguments are passed from `...` to [`apcluster::apcluster()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_apcluster = function(mat, s = apcluster::negDistMat(r = 2), ...) {
	
	check_pkg("apcluster", bioc = FALSE)

	x = apcluster::apcluster(s, mat, ...)
	cl = numeric(nrow(mat))
	for(i in seq_along(x@clusters)) {
		cl[x@clusters[[i]]] = i
	}
	cl
}


#' @param minPts Passed to the `minPts` argument in [`dbscan::hdbscan()`].
#' @details
#' `cluster_by_hdbscan()`: All other arguments are passed from `...` to [`dbscan::hdbscan()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_hdbscan = function(mat, minPts = 5, ...) {
	
	check_pkg("dbscan", bioc = FALSE)

	dbscan::hdbscan(mat, minPts = minPts, ...)$cluster
}

#' @param addLoops Passed to the `addLoops` argument in [`MCL::mcl()`].
#' @details
#' `cluster_by_MCL()`: All other arguments are passed from `...` to [`MCL::mcl()`].
#' 
#' @export
#' @rdname cluster_terms
cluster_by_MCL = function(mat, addLoops = TRUE, ...) {
	
	check_pkg("MCL", bioc = FALSE)

	MCL::mcl(mat, addLoops = addLoops, allow1 = TRUE, ...)$Cluster
}
