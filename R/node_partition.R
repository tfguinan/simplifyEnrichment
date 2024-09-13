

#' Partition the matrix
#'
#' @param mat The submatrix in the binary cut clustering process.
#' @param n_repeats Number of repeated runs of k-means clustering.
#'
#' @details
#' These functions can be set to the `partition_fun` argument in [`binary_cut()`].
#' 
#' `partition_by_kmeans()`: Since k-means clustering brings randomness, this function performs
#' k-means clustering several times (controlled by `n_repeats`) and uses the final consensus partitioning results.
#' 
#' @export
#' @import clue
#' @rdname node_partition
#' @return
#' All partitioning functions split the matrix into two groups and return a categorical vector of
#' labels of 1 and 2.
partition_by_kmeans = function(mat, n_repeats = 10) {
    partition_list = lapply(seq_len(n_repeats), function(i) {
        as.cl_hard_partition(kmeans(mat, 2))
    })
    partition_list = cl_ensemble(list = partition_list)
    partition_consensus = cl_consensus(partition_list)
    cl = as.vector(cl_class_ids(partition_consensus))
    if(length(unique(cl)) == 1) {
    	cl = partition_list[[1]]$.Data$cluster
    }
    cl
}

#' @details
#' `partition_by_pam()`: The clustering is performed by [`cluster::pam()`] with the `pamonce` argument set to 5.
#' 
#' @export
#' @importFrom cluster pam
#' @rdname node_partition
partition_by_pam = function(mat) {
    if(nrow(mat) > 10) {
        fit = pam(mat, 2, pamonce = 5)
    } else {
        fit = pam(mat, 2)
    }
    
    fit$clustering
}


#' @details
#' `partition_by_hclust()`: The "ward.D2" clusering method was used.
#'
#' @export
#' @rdname node_partition
partition_by_hclust = function(mat) {
    cutree(hclust(dist(mat), method = "ward.D2"), 2)
}

#' @details
#' `partition_by_kmeanspp()`: It uses the kmeanspp method from the **flexclust** package.
#' 
#' @export
#' @rdname node_partition
partition_by_kmeanspp = function(mat) {
    check_pkg("flexclust", bioc = FALSE)
    
    cl = flexclust::kcca(mat, k = 2, 
        family = flexclust::kccaFamily("kmeans"),
        control = list(initcent = "kmeanspp"))@cluster
    cl
}
